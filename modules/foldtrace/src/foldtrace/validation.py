"""Reproducible analytical campaign; no imported datasets or course solutions."""

import csv
import json
from dataclasses import asdict
from math import expm1, log1p, sqrt
from pathlib import Path

import numpy as np

from . import Truss, trace


def exact_force(q, ratio):
    """Independent direct strain/projection expression (NumPy)."""
    q = np.asarray(q)
    length = np.sqrt(1 + ratio**2 * (1 - q) ** 2)
    return -(length - sqrt(1 + ratio**2)) * (1 - q) / length


def exact_folds(ratio):
    """Closed-form roots of d(load)/dq=0; stable for shallow geometry."""
    offset = sqrt(expm1(log1p(ratio**2) / 3) / ratio**2)
    return np.array([1 - offset, 1 + offset])


def locate_folds(model, path):
    """Bracket folds using accepted states, then bisect the C++ tangent."""
    folds = []
    for left, right in zip(path.points[:-1], path.points[1:], strict=True):
        if left.tangent == 0:
            folds.append(left.q)
        elif left.tangent * right.tangent < 0:
            lo, hi = left.q, right.q
            for _ in range(48):
                mid = (lo + hi) / 2
                if model.tangent(lo) * model.tangent(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            folds.append((lo + hi) / 2)
    if path.points[-1].tangent == 0:
        folds.append(path.points[-1].q)
    return np.array(folds)


def campaign(output: str | Path, ratio: float = 0.2) -> dict:
    """Write CSV, JSON, and PNG; raise if analytical acceptance checks fail."""
    model = Truss(ratio)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    path = trace(model)
    q = np.array([p.q for p in path.points])
    load = np.array([p.load for p in path.points])
    numerical_folds = locate_folds(model, path)
    analytical_folds = exact_folds(ratio)
    if numerical_folds.shape != (2,):
        raise RuntimeError("campaign did not bracket two folds")
    rows = []
    for step in [0.16, 0.08, 0.04, 0.02]:
        fine = trace(model, step=step, max_step=step)
        states = np.array([p.q for p in fine.points])
        loads = np.array([p.load for p in fine.points])
        # Mid-cell errors measure a useful continuous reconstruction of the path.
        mid = (states[:-1] + states[1:]) / 2
        error = np.max(np.abs((loads[:-1] + loads[1:]) / 2 - exact_force(mid, ratio)))
        rows.append({"step": step, "points": len(states), "max_midpoint_error": float(error)})
    orders = np.log2(
        np.array([row["max_midpoint_error"] for row in rows[:-1]])
        / np.array([row["max_midpoint_error"] for row in rows[1:]])
    )
    summary = {
        "ratio": ratio,
        "points": len(path.points),
        "stop_reason": path.stop_reason,
        "rejected_steps": path.rejected_steps,
        "load_scale": path.load_scale,
        "max_scaled_equilibrium_residual": float(max(abs(p.residual) for p in path.points)),
        "max_hyperplane_residual": float(max(abs(p.constraint) for p in path.points)),
        "max_independent_force_error": float(np.max(np.abs(load - exact_force(q, ratio)))),
        "analytical_fold_q": analytical_folds.tolist(),
        "numerical_fold_q": numerical_folds.tolist(),
        "max_fold_location_error": float(np.max(np.abs(numerical_folds - analytical_folds))),
        "fold_loads": exact_force(analytical_folds, ratio).tolist(),
        "interpolation_orders": orders.tolist(),
        "note": (
            "Static symmetric equilibrium only; unstable states are not a dynamic trajectory."
        ),
    }
    if summary["max_scaled_equilibrium_residual"] > 1e-10:
        raise RuntimeError("equilibrium acceptance criterion failed")
    if summary["max_fold_location_error"] > 1e-9:
        raise RuntimeError("analytical fold acceptance criterion failed")
    if ratio == 0.2 and orders[-1] < 1.8:
        raise RuntimeError("default benchmark reconstruction order below 1.8")
    with (output / "path.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(path.points[0])))
        writer.writeheader()
        writer.writerows(asdict(p) for p in path.points)
    with (output / "convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "validation.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    _plot(output, model, path, rows, analytical_folds)
    return summary


def _plot(output, model, path, rows, folds):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ratio = model.ratio
    grid = np.linspace(0, 2.3, 900)
    loads = exact_force(grid, ratio)
    unstable = (grid > folds[0]) & (grid < folds[1])
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    fig.suptitle(
        "FoldTrace | follow equilibrium beyond the peak load", fontsize=19, weight="bold"
    )
    ax = axes[0, 0]
    ax.plot(grid, loads, color="#156a76", lw=2, label="exact equilibrium")
    ax.plot(grid[unstable], loads[unstable], color="#d16037", lw=3, label="negative stiffness")
    ax.scatter(
        [p.q for p in path.points],
        [p.load for p in path.points],
        s=12,
        color="#162f45",
        label="continuation states",
        zorder=4,
    )
    ax.scatter(folds, exact_force(folds, ratio), marker="*", s=150, color="#c39722", zorder=5)
    ax.axhline(0, color="0.75", lw=0.7)
    ax.set(
        xlabel="Downward displacement q = w/h",
        ylabel="Load P / (2 EA h/L₀)",
        title="Both limit points crossed",
    )
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for displacement, color, label in [
        (0, "#156a76", "initial"),
        (1, "#d16037", "flat"),
        (2, "#162f45", "inverted"),
    ]:
        ax.plot(
            [-1, 0, 1], [0, ratio * (1 - displacement), 0], "o-", color=color, lw=2, label=label
        )
    ax.set(xlabel="x / a", ylabel="y / a", title=f"Two axial bars | h/a = {ratio:g}")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(fontsize=8, loc="upper right")
    ax = axes[1, 0]
    steps = np.array([row["step"] for row in rows])
    errors = np.array([row["max_midpoint_error"] for row in rows])
    ax.loglog(steps, errors, "o-", color="#156a76", label="measured")
    ax.loglog(
        steps, errors[-1] * (steps / steps[-1]) ** 2, "--", color="0.4", label="second order"
    )
    ax.set(
        xlabel="Maximum step in scaled path coordinates",
        ylabel="Maximum load interpolation error",
        title="Piecewise linear path reconstruction",
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2, which="both")
    ax = axes[1, 1]
    peak = float(exact_force(folds[0], ratio))
    applied = 0.6 * peak
    coefficient = 2 * ratio**2 / sqrt(1 + ratio**2)
    potential = np.array([model.energy(float(q)) for q in grid]) - coefficient * applied * grid
    ax.plot(grid, potential, color="#162f45", lw=2)
    ax.set(
        xlabel="Downward displacement q = w/h",
        ylabel="Potential / (EA a)",
        title="Energy landscape at 60% of the peak load",
    )
    ax.text(
        0.03,
        0.06,
        "Stability refers only to the constrained symmetric mode.\n"
        "No inertia, damping, or dynamic jump is simulated.",
        transform=ax.transAxes,
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    fig.savefig(output / "equilibrium_audit.png", dpi=170)
    plt.close(fig)
