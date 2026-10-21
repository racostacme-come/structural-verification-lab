"""Adaptive predictor/corrector orchestration; mechanics and Newton run in C++."""

from dataclasses import dataclass
from math import hypot, isfinite

from ._core import Truss, correct


@dataclass(frozen=True)
class Point:
    q: float
    load: float
    tangent: float
    energy: float
    step: float
    iterations: int
    retries: int
    residual: float
    constraint: float


@dataclass(frozen=True)
class Path:
    points: tuple[Point, ...]
    load_scale: float
    rejected_steps: int
    stop_reason: str


class ContinuationError(RuntimeError):
    """Continuation stopped; ``path`` retains only accepted equilibrium points."""

    def __init__(self, reason: str, path: Path):
        super().__init__(reason)
        self.path = path


def trace(
    model: Truss,
    *,
    q_stop: float = 2.2,
    step: float = 0.08,
    min_step: float = 1e-5,
    max_step: float = 0.12,
    load_scale: float | None = None,
    tolerance: float = 1e-11,
    max_iterations: int = 12,
    max_steps: int = 10000,
) -> Path:
    """Trace forward from the unloaded state until q >= q_stop.

    Uses a tangent-plane pseudo-arclength constraint in (q, load/load_scale).
    The final point can overshoot q_stop. Each rejected attempt halves the step.
    A correction > half the predictor length is rejected to limit branch jumps.
    """
    scale = model.ratio**2 if load_scale is None else load_scale
    values = (q_stop, step, min_step, max_step, scale, tolerance)
    if not all(isfinite(v) for v in values):
        raise ValueError("continuation parameters must be finite")
    if not 0 < q_stop <= 3 or not 0 < min_step <= step <= max_step <= 0.5:
        raise ValueError("require q_stop in (0, 3] and 0 < min_step <= step <= max_step <= 0.5")
    if scale <= 0 or not 1e-14 <= tolerance <= 1e-4:
        raise ValueError("load_scale must be positive; tolerance must be in [1e-14, 1e-4]")
    if any(type(v) is not int or v < 1 for v in (max_iterations, max_steps)):
        raise ValueError("iteration and step budgets must be positive integers")
    points = [Point(0.0, 0.0, model.tangent(0), 0.0, 0.0, 0, 0, 0.0, 0.0)]
    rejected = 0
    previous_tangent = (1.0, 0.0)

    def failure(reason):
        raise ContinuationError(reason, Path(tuple(points), scale, rejected, reason))

    for _ in range(max_steps):
        last = points[-1]
        slope = last.tangent / scale
        norm = hypot(1.0, slope)
        if not isfinite(norm):
            failure("ill_conditioned_scaling")
        tq, tl = 1.0 / norm, slope / norm
        if tq * previous_tangent[0] + tl * previous_tangent[1] < 0:
            tq, tl = -tq, -tl
        retries = 0
        while True:
            qp, lp = last.q + step * tq, last.load + scale * step * tl
            result = correct(model, qp, lp, tq, tl, scale, tolerance, max_iterations)
            correction = hypot(result.q - qp, (result.load - lp) / scale)
            forward = result.q > last.q
            if result.converged and correction <= 0.5 * step and forward:
                break
            rejected += 1
            retries += 1
            if step / 2 < min_step:
                failure("minimum_step_reached")
            step /= 2
        points.append(
            Point(
                result.q,
                result.load,
                model.tangent(result.q),
                model.energy(result.q),
                step,
                result.iterations,
                retries,
                result.residual,
                result.constraint,
            )
        )
        if result.q >= q_stop:
            return Path(tuple(points), scale, rejected, "target_reached")
        previous_tangent = (tq, tl)
        if result.iterations <= 3:
            step = min(max_step, 1.25 * step)
        elif result.iterations >= 7:
            step = max(min_step, 0.7 * step)
    failure("step_budget_exhausted")
