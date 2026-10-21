"""P1 axial bar assembly with committed-state constitutive evaluations."""

from dataclasses import dataclass

import numpy as np

from ._core import Material, State


class ConvergenceError(RuntimeError):
    """An increment failed; its trial states were not committed."""


def _readonly(values):
    array = np.array(values, dtype=float, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class Bar:
    """One material per element, positive lengths, constant positive area.

    Use any consistent unit system. The examples use mm, N, and MPa.
    Strictly positive total hardening is required for invertible FEM tangents.
    """

    lengths: np.ndarray
    materials: tuple[Material, ...]
    area: float = 1.0

    def __post_init__(self):
        lengths = _readonly(self.lengths)
        materials = tuple(self.materials)
        if (
            lengths.ndim != 1
            or len(lengths) == 0
            or not np.all(np.isfinite(lengths))
            or np.any(lengths <= 0)
            or not np.isfinite(np.sum(lengths))
            or len(materials) != len(lengths)
            or not np.isfinite(self.area)
            or self.area <= 0
        ):
            raise ValueError(
                "Require positive finite lengths/area and one material per element."
            )
        if any(
            not isinstance(m, Material) or m.isotropic + m.kinematic <= 0 for m in materials
        ):
            raise ValueError("FEM requires Material objects with positive total hardening.")
        object.__setattr__(self, "lengths", lengths)
        object.__setattr__(self, "materials", materials)


@dataclass(frozen=True)
class Step:
    """A converged increment; iterations counts linear Newton corrections."""

    displacement: np.ndarray
    states: tuple[State, ...]
    stress: np.ndarray
    iterations: int
    residual_history: tuple[float, ...]


def solve_step(
    bar: Bar,
    end_displacement: float,
    previous: Step | None = None,
    *,
    tangent: str = "consistent",
    max_iterations: int = 100,
    atol: float = 1e-10,
    rtol: float = 1e-10,
) -> Step:
    """Solve u(0)=0, u(L)=end_displacement, no body force.

    Every trial uses the same previous internal states. Failures leave previous
    unchanged. An elastic tangent is exposed only as an iteration benchmark.
    Tolerances are force units (atol) and dimensionless (rtol). Previous must
    come from this same bar and use the same units.
    """
    if (
        tangent not in {"consistent", "elastic"}
        or isinstance(max_iterations, bool)
        or not isinstance(max_iterations, int)
        or max_iterations < 0
        or not np.isfinite(end_displacement)
        or not np.isfinite(atol)
        or not np.isfinite(rtol)
        or atol <= 0
        or rtol < 0
    ):
        raise ValueError("Invalid Newton settings or end displacement.")
    n = len(bar.lengths)
    old = tuple(State() for _ in range(n)) if previous is None else previous.states
    u = np.zeros(n + 1) if previous is None else np.array(previous.displacement, copy=True)
    if len(old) != n or u.shape != (n + 1,) or not np.all(np.isfinite(u)) or u[0] != 0:
        raise ValueError("Previous step is incompatible with this bar.")
    x = np.concatenate(([0.0], np.cumsum(bar.lengths)))
    u += (end_displacement - u[-1]) * x / x[-1]
    u[-1] = end_displacement
    residuals = []
    for iteration in range(max_iterations + 1):
        strain = np.diff(u) / bar.lengths
        responses = [m.update(e, s) for m, e, s in zip(bar.materials, strain, old, strict=True)]
        stress = np.array([r.stress for r in responses])
        force = bar.area * stress
        residual = force[:-1] - force[1:]
        norm = float(np.max(np.abs(residual), initial=0))
        residuals.append(norm)
        if norm <= atol + rtol * max(1.0, float(np.max(np.abs(force)))):
            return Step(
                _readonly(u),
                tuple(r.state for r in responses),
                _readonly(stress),
                iteration,
                tuple(residuals),
            )
        if iteration == max_iterations:
            break
        moduli = np.array(
            [
                r.tangent if tangent == "consistent" else m.elastic
                for r, m in zip(responses, bar.materials, strict=True)
            ]
        )
        k = bar.area * moduli / bar.lengths
        matrix = np.diag(k[:-1] + k[1:])
        matrix += np.diag(-k[1:-1], 1) + np.diag(-k[1:-1], -1)
        try:
            correction = np.linalg.solve(matrix, -residual)
        except np.linalg.LinAlgError as exc:
            raise ConvergenceError(
                "Singular tangent; previous state remains unchanged."
            ) from exc
        # Backtracking prevents active-set cycles when a trial lands on a yield kink.
        for reduction in range(30):
            factor = 0.5**reduction
            candidate = u.copy()
            candidate[1:-1] += factor * correction
            trial_strain = np.diff(candidate) / bar.lengths
            trial_force = bar.area * np.array(
                [
                    m.update(e, s).stress
                    for m, e, s in zip(bar.materials, trial_strain, old, strict=True)
                ]
            )
            trial_norm = float(np.max(np.abs(np.diff(trial_force)), initial=0))
            threshold = atol + rtol * max(1.0, float(np.max(np.abs(trial_force))))
            if trial_norm <= (1 - 1e-4 * factor) * norm or trial_norm <= threshold:
                u = candidate
                break
        else:
            raise ConvergenceError("Line search failed; previous state remains unchanged.")
    raise ConvergenceError(
        f"No convergence after {max_iterations} corrections; residual={norm:g}."
    )
