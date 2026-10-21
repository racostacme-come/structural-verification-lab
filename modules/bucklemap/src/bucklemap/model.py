"""Dense Hermite FEM for conservative, constant-force Euler column buckling."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from ._core import element_matrices


def _mesh(nodes):
    x = np.asarray(nodes, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.all(np.isfinite(x)):
        raise ValueError("nodes must be a finite 1D array of length >= 2")
    if not np.all(np.diff(x) > 0):
        raise ValueError("nodes must be strictly increasing")
    return x


def assemble(nodes, rigidity=1.0):
    """Return full (K, G); EI may be scalar or one positive value per element.

    Units: nodes [length], EI [force*length^2]; eigenvalues are forces.
    G is the slope-product integral for a unit uniform compressive force.
    """
    x = _mesh(nodes)
    ei = np.asarray(rigidity, dtype=float)
    if ei.ndim == 0:
        ei = np.full(len(x) - 1, float(ei))
    if ei.shape != (len(x) - 1,) or not np.all(np.isfinite(ei)) or np.any(ei <= 0):
        raise ValueError("EI must be positive finite scalar or one value per element")
    k = np.zeros((2 * len(x), 2 * len(x)))
    g = np.zeros_like(k)
    for i, (h, value) in enumerate(zip(np.diff(x), ei, strict=True)):
        ke, ge = element_matrices(float(h), float(value))
        sl = slice(2 * i, 2 * i + 4)
        k[sl, sl] += ke
        g[sl, sl] += ge
    if not np.all(np.isfinite(k)) or not np.all(np.isfinite(g)):
        raise ValueError("assembled matrices exceed floating-point range")
    return k, g


@dataclass(frozen=True)
class BucklingResult:
    """Sorted loads, G-normalized full modes, and reduced-space backward errors.

    `modes[:, j]` interleaves displacement/slope. Fixed DOFs are exactly zero.
    Backward errors are NOT bounds on continuum or physical error.
    """

    nodes: np.ndarray
    loads: np.ndarray
    modes: np.ndarray
    residuals: np.ndarray
    free_dofs: np.ndarray


def solve(nodes, rigidity=1.0, *, support="pinned", springs=(0.0, 0.0), modes=1):
    """Solve K q = P G q after exact constraint elimination.

    support: 'pinned' (w=0 both ends), 'clamped' (w=theta=0 both),
    'cantilever' (left w=theta=0), 'clamped-pinned' (left fixed, right w=0).
    springs: finite nonnegative (left, right) rotational stiffness [force*length].
    A nonzero spring at an already fixed rotation is rejected.
    Dense O(ndof^3) solve intended for small studies, not large-scale meshes.
    """
    x = _mesh(nodes)
    k, g = assemble(x, rigidity)
    end = 2 * len(x) - 2
    constraints = {
        "pinned": [0, end],
        "clamped": [0, 1, end, end + 1],
        "cantilever": [0, 1],
        "clamped-pinned": [0, 1, end],
    }
    if support not in constraints:
        raise ValueError("unknown support preset")
    fixed = constraints[support]
    spring = np.asarray(springs, dtype=float)
    if spring.shape != (2,) or not np.all(np.isfinite(spring)) or np.any(spring < 0):
        raise ValueError("springs must contain two finite nonnegative stiffnesses")
    for dof, value in zip((1, end + 1), spring, strict=True):
        if dof in fixed and value != 0:
            raise ValueError("a rotational spring cannot act at an already fixed rotation")
        k[dof, dof] += value
    free = np.setdiff1d(np.arange(len(k)), fixed)
    if (
        isinstance(modes, bool)
        or not isinstance(modes, Integral)
        or not 1 <= modes <= len(free)
    ):
        raise ValueError(f"modes must be an integer in [1, {len(free)}]")
    kr, gr = k[np.ix_(free, free)], g[np.ix_(free, free)]
    if not np.all(np.isfinite(kr)):
        raise ValueError("spring stiffness overflows assembled matrix")
    # Congruence scaling improves mixed displacement/slope conditioning.
    scale = 1 / np.sqrt(np.diag(kr))
    ks, gs = scale[:, None] * kr * scale, scale[:, None] * gr * scale
    loads, vectors = eigh(ks, gs, subset_by_index=(0, modes - 1), driver="gvx")
    vectors *= scale[:, None]
    if not np.all(np.isfinite(loads)) or np.any(loads <= 0):
        raise ArithmeticError("nonpositive or nonfinite load; check mesh and conditioning")
    full = np.zeros((len(k), modes))
    full[free] = vectors
    for j in range(modes):
        pivot = np.argmax(np.abs(full[::2, j])) * 2
        if full[pivot, j] < 0:
            full[:, j] *= -1
    residuals = np.empty(modes)
    for j, load in enumerate(loads):
        q = full[free, j]
        denominator = np.linalg.norm(kr, np.inf) + load * np.linalg.norm(gr, np.inf)
        residuals[j] = np.linalg.norm(kr @ q - load * gr @ q, np.inf) / (
            denominator * np.linalg.norm(q, np.inf)
        )
    return BucklingResult(x.copy(), loads, full, residuals, free)


def reconstruct(nodes, dofs, points):
    """Evaluate cubic Hermite displacement; no extrapolation or renormalization."""
    x = _mesh(nodes)
    q = np.asarray(dofs, dtype=float)
    p = np.asarray(points, dtype=float)
    if q.shape != (2 * len(x),) or not np.all(np.isfinite(q)):
        raise ValueError("dofs must be a finite vector with displacement/slope at every node")
    if not np.all(np.isfinite(p)) or np.any(p < x[0]) or np.any(p > x[-1]):
        raise ValueError("points must be finite and inside the mesh")
    i = np.clip(np.searchsorted(x, p, side="right") - 1, 0, len(x) - 2)
    h = x[i + 1] - x[i]
    t = (p - x[i]) / h
    return (
        (1 - 3 * t**2 + 2 * t**3) * q[2 * i]
        + h * (t - 2 * t**2 + t**3) * q[2 * i + 1]
        + (3 * t**2 - 2 * t**3) * q[2 * i + 2]
        + h * (-(t**2) + t**3) * q[2 * i + 3]
    )


def spring_reference(kappa):
    """Exact nondimensional first load P*L^2/EI for equal end springs.

    Uniform EI, w(0)=w(L)=0, kappa=k_theta*L/EI. The symmetric mode
    gives kappa*sin(z/2)+z*cos(z/2)=0, pi<=z<=2*pi; return z^2.
    This scalar continuum equation is independent of the FEM matrices.
    """
    value = float(kappa)
    if not np.isfinite(value) or value < 0:
        raise ValueError("kappa must be finite and nonnegative")
    if value == 0:
        return np.pi**2
    # Divide by kappa+1 without overflowing for very large finite kappa.
    a = 1 / (1 + value)
    b = value * a
    z = brentq(lambda z: b * np.sin(z / 2) + a * z * np.cos(z / 2), np.pi, 2 * np.pi)
    return z**2
