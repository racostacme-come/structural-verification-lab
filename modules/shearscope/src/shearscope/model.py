"""SI-unit beam models, sparse assembly and cantilever boundary conditions."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from . import _core


@dataclass(frozen=True)
class Beam:
    """Uniform rectangular section; lengths in m, modulus in Pa, density in kg/m³."""

    length: float = 1.0
    width: float = 0.02
    height: float = 0.01
    young: float = 70e9
    poisson: float = 0.3
    density: float = 2700.0
    kappa: float = 5.0 / 6.0

    def __post_init__(self):
        values = (self.length, self.width, self.height, self.young, self.density, self.kappa)
        if not all(np.isfinite(x) and x > 0 for x in values):
            raise ValueError(
                "dimensions, modulus, density and kappa must be finite and positive"
            )
        if not np.isfinite(self.poisson) or not -1 < self.poisson < 0.5:
            raise ValueError("poisson must lie strictly between -1 and 0.5")

    @property
    def area(self):
        return self.width * self.height

    @property
    def inertia(self):
        return self.width * self.height**3 / 12

    @property
    def ei(self):
        return self.young * self.inertia

    @property
    def shear(self):
        return self.kappa * self.young / (2 * (1 + self.poisson)) * self.area


def _count(value, name, upper):
    if isinstance(value, bool) or not isinstance(value, Integral) or not 1 <= value <= upper:
        raise ValueError(f"{name} must be an integer between 1 and {upper}")


def assemble(beam: Beam, elements: int = 32, integration: str = "reduced", nodes=None):
    """Return (coordinates, CSR stiffness, CSR consistent mass).

    Optional strictly increasing nodes must span [0, beam.length]; they override
    `elements`. DOFs are [w0, theta0, w1, theta1, ...]. No constraints are applied.
    """
    if integration not in ("full", "reduced"):
        raise ValueError("integration must be 'full' or 'reduced'")
    if nodes is None:
        _count(elements, "elements", 10000)
        x = np.linspace(0, beam.length, elements + 1)
    else:
        x = np.asarray(nodes, dtype=float)
        if x.ndim != 1 or not 2 <= x.size <= 10001 or not np.all(np.isfinite(x)):
            raise ValueError("nodes must be a finite vector of 2 to 10001 coordinates")
        if x[0] != 0 or x[-1] != beam.length or np.any(np.diff(x) <= 0):
            raise ValueError("nodes must increase strictly from 0 to beam.length")
    rows, cols, kvals, mvals = [], [], [], []
    for e, length in enumerate(np.diff(x)):
        k, m = _core.element(
            length,
            beam.ei,
            beam.shear,
            beam.density * beam.area,
            beam.density * beam.inertia,
            integration == "reduced",
        )
        dofs = np.arange(2 * e, 2 * e + 4)
        rows.extend(np.repeat(dofs, 4))
        cols.extend(np.tile(dofs, 4))
        kvals.extend(k.ravel())
        mvals.extend(m.ravel())
    shape = (2 * x.size, 2 * x.size)
    return (
        x.copy(),
        coo_matrix((kvals, (rows, cols)), shape=shape).tocsr(),
        coo_matrix((mvals, (rows, cols)), shape=shape).tocsr(),
    )


@dataclass(frozen=True)
class StaticResult:
    nodes: np.ndarray
    displacement: np.ndarray
    rotation: np.ndarray
    reactions: np.ndarray
    strain_energy: float
    relative_residual: float


def solve_static(
    beam: Beam,
    elements=32,
    integration="reduced",
    force=1.0,
    moment=0.0,
    distributed=0.0,
    nodes=None,
):
    """Clamp the left end; apply tip force/moment and uniform transverse load.

    Loads use N, N·m, N/m. Reactions are [root shear, root moment].
    The residual is ||Kff uf - ff||₂ / max(||ff||₂, 1 N in numeric units).
    """
    if not all(np.isfinite(v) for v in (force, moment, distributed)):
        raise ValueError("loads must be finite")
    x, k, _ = assemble(beam, elements, integration, nodes)
    f = np.zeros(k.shape[0])
    for e, length in enumerate(np.diff(x)):
        f[2 * e] += distributed * length / 2
        f[2 * e + 2] += distributed * length / 2
    f[-2:] += [force, moment]
    u = np.zeros_like(f)
    u[2:] = spsolve(k[2:, 2:], f[2:])
    if not np.all(np.isfinite(u)):
        raise ArithmeticError("nonfinite solution; check mesh and physical scaling")
    residual = k @ u - f
    return StaticResult(
        x,
        u[::2],
        u[1::2],
        residual[:2],
        float(0.5 * u @ (k @ u)),
        float(np.linalg.norm(residual[2:]) / max(np.linalg.norm(f[2:]), 1.0)),
    )


def exact_tip(beam: Beam, force=1.0, moment=0.0, distributed=0.0):
    """Closed-form Timoshenko tip displacement for the supported load family."""
    if not all(np.isfinite(v) for v in (force, moment, distributed)):
        raise ValueError("loads must be finite")
    length = beam.length
    return (
        force * (length**3 / (3 * beam.ei) + length / beam.shear)
        + moment * length**2 / (2 * beam.ei)
        + distributed * (length**4 / (8 * beam.ei) + length**2 / (2 * beam.shear))
    )


def modes(beam: Beam, elements=32, integration="reduced", count=3):
    """Return frequencies in Hz and full mass-normalized mode columns.

    Dense symmetric generalized eigensolve, limited to 500 elements. Mode signs
    are chosen so the largest transverse displacement is positive.
    """
    _count(elements, "elements", 500)
    _count(count, "count", 2 * elements)
    _, k, m = assemble(beam, elements, integration)
    values, vectors = eigh(
        k[2:, 2:].toarray(), m[2:, 2:].toarray(), subset_by_index=(0, count - 1), driver="gvx"
    )
    if np.any(values <= 0) or not np.all(np.isfinite(values)):
        raise ArithmeticError("nonpositive eigenvalues; check numerical conditioning")
    shapes = np.zeros((k.shape[0], count))
    shapes[2:, :] = vectors
    for j in range(count):
        peak = np.argmax(np.abs(shapes[::2, j]))
        if shapes[2 * peak, j] < 0:
            shapes[:, j] *= -1
    return np.sqrt(values) / (2 * np.pi), shapes
