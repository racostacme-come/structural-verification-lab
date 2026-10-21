"""P1 Galerkin solve, residual estimation, and bulk-marked bisection on [0, 1]."""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import solve_banded

FloatArray = NDArray[np.float64]
Field = Callable[[FloatArray], ArrayLike]


def _integer(value: int, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _positive(value: float, name: str) -> float:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def mesh(nodes: ArrayLike) -> FloatArray:
    """Validate and copy a strictly increasing mesh with endpoints exactly 0 and 1."""
    x = np.array(nodes, dtype=float, copy=True)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)):
        raise ValueError("nodes must be a finite one-dimensional array with >= 2 entries")
    if x[0] != 0 or x[-1] != 1 or np.any(np.diff(x) <= 0):
        raise ValueError("nodes must increase strictly from 0 to 1")
    return x


def _quadrature(x: FloatArray, order: int) -> tuple[FloatArray, FloatArray, FloatArray]:
    _integer(order, "quadrature order", 2)
    if order > 32:
        raise ValueError("quadrature order must be <= 32")
    z, w = np.polynomial.legendre.leggauss(order)
    t = (z + 1) / 2
    points = x[:-1, None] + np.diff(x)[:, None] * t
    weights = np.diff(x)[:, None] * w / 2
    return points, weights, t


def _values(field: Field, points: FloatArray) -> FloatArray:
    values = np.asarray(field(points), dtype=float)
    try:
        values = np.broadcast_to(values, points.shape)
    except ValueError as exc:
        raise ValueError("field must return a scalar or values matching input shape") from exc
    if not np.all(np.isfinite(values)):
        raise ValueError("field returned nonfinite values")
    return values


@dataclass(frozen=True)
class Solution:
    """Read-only nodal solution, cell flux, and squared residual contributions."""

    nodes: FloatArray
    values: FloatArray
    flux: FloatArray
    indicators_squared: FloatArray
    diffusion: float

    @property
    def estimator(self) -> float:
        return float(np.sqrt(np.sum(self.indicators_squared)))

    @property
    def elements(self) -> int:
        return self.nodes.size - 1

    def evaluate(self, points: ArrayLike) -> FloatArray:
        """Interpolate within the closed domain; extrapolation is rejected."""
        p = np.asarray(points, dtype=float)
        if not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
            raise ValueError("evaluation points must lie in [0, 1]")
        return np.asarray(np.interp(p, self.nodes, self.values))


def solve(
    nodes: ArrayLike,
    load: Field,
    *,
    diffusion: float = 1.0,
    boundary: tuple[float, float] = (0.0, 0.0),
    quadrature: int = 8,
) -> Solution:
    """Solve -a*u''=load with constant a>0 and prescribed endpoint values.

    Fields receive arrays of shape (elements, quadrature). Assembly and the
    tridiagonal direct solve use O(elements) storage and work for fixed quadrature.
    The estimator is an uncalibrated residual indicator, not a certified bound.
    """
    x = mesh(nodes)
    a = _positive(diffusion, "diffusion")
    bc = np.asarray(boundary, dtype=float)
    if bc.shape != (2,) or not np.all(np.isfinite(bc)):
        raise ValueError("boundary must contain two finite endpoint values")
    points, weights, t = _quadrature(x, quadrature)
    f = _values(load, points)
    h = np.diff(x)
    stiffness = a / h
    rhs = np.zeros(x.size)
    rhs[:-1] += np.sum(weights * f * (1 - t), axis=1)
    rhs[1:] += np.sum(weights * f * t, axis=1)
    u = np.zeros(x.size)
    u[[0, -1]] = bc
    if x.size > 2:
        interior = rhs[1:-1].copy()
        interior[0] += stiffness[0] * bc[0]
        interior[-1] += stiffness[-1] * bc[1]
        bands = np.zeros((3, interior.size))
        bands[1] = stiffness[:-1] + stiffness[1:]
        bands[0, 1:] = -stiffness[1:-1]
        bands[2, :-1] = -stiffness[1:-1]
        u[1:-1] = solve_banded((1, 1), bands, interior)
    flux = a * np.diff(u) / h
    indicators = h**2 / a * np.sum(weights * f**2, axis=1)
    jump_squared = np.diff(flux) ** 2
    indicators[:-1] += h[:-1] / (2 * a) * jump_squared
    indicators[1:] += h[1:] / (2 * a) * jump_squared
    for array in (x, u, flux, indicators):
        if not np.all(np.isfinite(array)):
            raise FloatingPointError("nonfinite solution or indicator; rescale the problem")
        array.flags.writeable = False
    return Solution(x, u, flux, indicators, a)


def mark(indicators_squared: ArrayLike, theta: float = 0.5) -> NDArray[np.int64]:
    """Smallest descending prefix containing theta of squared indicator mass.

    Equal contributions prefer the lower cell index. Zero total marks no cells.
    theta=1 marks every positive contribution, even across large dynamic ranges.
    """
    eta = np.asarray(indicators_squared, dtype=float)
    if eta.ndim != 1 or not np.all(np.isfinite(eta)) or np.any(eta < 0):
        raise ValueError("indicators must be a finite nonnegative vector")
    if not np.isfinite(theta) or not 0 < theta <= 1:
        raise ValueError("theta must lie in (0, 1]")
    order = np.argsort(-eta, kind="stable")
    order = order[eta[order] > 0]
    if order.size == 0 or theta == 1:
        return order
    scaled = eta[order] / eta[order[0]]
    cumulative = np.cumsum(scaled)
    count = np.searchsorted(cumulative, theta * cumulative[-1], side="left") + 1
    return order[:count]


def refine(nodes: ArrayLike, marked: ArrayLike) -> FloatArray:
    """Bisect selected distinct cells, preserving every existing node."""
    x = mesh(nodes)
    chosen = np.asarray(marked)
    if chosen.ndim != 1:
        raise ValueError("marked cells must be a one-dimensional integer array")
    if chosen.size == 0:
        return x
    if chosen.dtype.kind not in "iu" or np.any(chosen < 0) or np.any(chosen >= x.size - 1):
        raise ValueError("marked cell index out of range or not integer")
    chosen = np.unique(chosen)
    mids = x[chosen] + (x[chosen + 1] - x[chosen]) / 2
    if np.any(mids <= x[chosen]) or np.any(mids >= x[chosen + 1]):
        raise FloatingPointError("bisection reached floating-point resolution")
    return np.sort(np.concatenate((x, mids)))


def errors(solution: Solution, exact: Field, gradient: Field, *, quadrature: int = 16):
    """Return (L2 error, energy error), independently integrated from the solve."""
    points, weights, t = _quadrature(solution.nodes, quadrature)
    uh = solution.values[:-1, None] * (1 - t) + solution.values[1:, None] * t
    l2 = np.sum(weights * (uh - _values(exact, points)) ** 2)
    slope = solution.flux[:, None] / solution.diffusion
    energy = solution.diffusion * np.sum(weights * (slope - _values(gradient, points)) ** 2)
    return float(np.sqrt(l2)), float(np.sqrt(energy))


@dataclass(frozen=True)
class AdaptiveResult:
    """Accepted solutions and the explicit reason refinement stopped."""

    history: tuple[Solution, ...]
    stop_reason: str

    @property
    def solution(self) -> Solution:
        return self.history[-1]


def adapt(
    load: Field,
    *,
    initial_nodes: ArrayLike | None = None,
    diffusion: float = 1.0,
    boundary: tuple[float, float] = (0.0, 0.0),
    theta: float = 0.5,
    tolerance: float = 0.1,
    max_elements: int = 4096,
    max_steps: int = 40,
    quadrature: int = 8,
) -> AdaptiveResult:
    """Solve-estimate-mark-refine until tolerance or a resource limit.

    The full bulk-marked set is accepted or rejected: no partial marking is used
    to squeeze into the cell budget. Tolerance applies to the residual estimator.
    """
    _positive(tolerance, "tolerance")
    _integer(max_elements, "max_elements", 1)
    _integer(max_steps, "max_steps", 1)
    mark([], theta)  # Validate even if no refinement is needed.
    x = mesh(np.linspace(0, 1, 9) if initial_nodes is None else initial_nodes)
    if x.size - 1 > max_elements:
        raise ValueError("initial mesh exceeds max_elements")
    history = []
    for step in range(max_steps):
        current = solve(x, load, diffusion=diffusion, boundary=boundary, quadrature=quadrature)
        history.append(current)
        if current.estimator <= tolerance:
            return AdaptiveResult(tuple(history), "tolerance")
        if step == max_steps - 1:
            return AdaptiveResult(tuple(history), "max_steps")
        selected = mark(current.indicators_squared, theta)
        if current.elements + selected.size > max_elements:
            return AdaptiveResult(tuple(history), "max_elements")
        x = refine(x, selected)
    raise AssertionError("unreachable")
