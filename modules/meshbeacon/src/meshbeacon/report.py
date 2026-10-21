"""Deterministic analytical validation and an equal-cell-budget experiment."""

import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .core import adapt, errors, solve  # noqa: E402
from .problems import bump, sine  # noqa: E402


def _row(strategy, step, solution, problem):
    l2, energy = errors(solution, problem.exact, problem.gradient)
    return {
        "problem": problem.name,
        "strategy": strategy,
        "step": step,
        "elements": solution.elements,
        "free_dofs": solution.elements - 1,
        "l2_error": l2,
        "energy_error": energy,
        "estimator": solution.estimator,
        "effectivity": solution.estimator / energy,
    }


def _csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def benchmark(output: str | Path) -> dict:
    """Write the fixed, reproducible benchmark; raise if numerical checks fail.

    Existing files with the four documented output names are replaced. Other
    files are untouched. Wall-clock timing is intentionally not a quality metric.
    """
    smooth = sine()
    rows = []
    for step, n in enumerate((8, 16, 32, 64, 128, 256)):
        rows.append(_row("uniform", step, solve(np.linspace(0, 1, n + 1), smooth.load), smooth))
    smooth_rows = rows.copy()
    localized = bump()
    result = adapt(localized.load, tolerance=0.15)
    adaptive_rows, uniform_rows = [], []
    for step, solution in enumerate(result.history):
        adaptive_rows.append(_row("adaptive", step, solution, localized))
        uniform = solve(np.linspace(0, 1, solution.elements + 1), localized.load)
        uniform_rows.append(_row("uniform_equal_cells", step, uniform, localized))
    rows.extend(adaptive_rows)
    rows.extend(uniform_rows)
    final = result.solution
    higher_q = solve(final.nodes, localized.load, quadrature=16)
    energy24 = errors(final, localized.exact, localized.gradient, quadrature=24)[1]
    energy_ratio = adaptive_rows[-1]["energy_error"] / uniform_rows[-1]["energy_error"]
    summary = {
        "parameters": {
            "center": 0.37,
            "width": 0.04,
            "theta": 0.5,
            "tolerance": 0.15,
            "assembly_quadrature": 8,
            "error_quadrature": 16,
            "initial_elements": 8,
        },
        "stop_reason": result.stop_reason,
        "adaptive_solves": len(result.history),
        "final_elements": final.elements,
        "final_free_dofs": final.elements - 1,
        "adaptive_energy_error": adaptive_rows[-1]["energy_error"],
        "uniform_equal_cells_energy_error": uniform_rows[-1]["energy_error"],
        "adaptive_to_uniform_energy_ratio": energy_ratio,
        "final_estimator": final.estimator,
        "final_effectivity": adaptive_rows[-1]["effectivity"],
        "sine_finest_l2_order": float(np.log2(rows[4]["l2_error"] / rows[5]["l2_error"])),
        "sine_finest_energy_order": float(
            np.log2(rows[4]["energy_error"] / rows[5]["energy_error"])
        ),
        "quadrature_8_to_16_max_nodal_change": float(
            np.max(abs(final.values - higher_q.values))
        ),
        "error_quadrature_16_to_24_energy_change": abs(
            energy24 - adaptive_rows[-1]["energy_error"]
        ),
    }
    checks = {
        "tolerance_reached": result.stop_reason == "tolerance",
        "l2_order": abs(summary["sine_finest_l2_order"] - 2) < 0.01,
        "energy_order": abs(summary["sine_finest_energy_order"] - 1) < 0.01,
        "adaptive_advantage": energy_ratio < 0.4,
        "load_quadrature_stability": summary["quadrature_8_to_16_max_nodal_change"] < 1e-7,
        "error_quadrature_stability": summary["error_quadrature_16_to_24_energy_change"] < 1e-9,
    }
    summary["checks"] = checks
    if not all(checks.values()):
        raise RuntimeError(f"benchmark validation failed: {checks}")
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    _csv(destination / "convergence.csv", rows)
    _csv(
        destination / "adaptive_nodes.csv",
        [
            {"x": float(x), "u_h": float(u), "u_exact": float(localized.exact(x))}
            for x, u in zip(final.nodes, final.values, strict=True)
        ],
    )
    (destination / "validation.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    _plot(
        destination / "mesh_audit.png",
        final,
        localized,
        adaptive_rows,
        uniform_rows,
        smooth_rows,
    )
    return summary


def _plot(path, final, problem, adaptive, uniform, smooth):
    with plt.rc_context(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "#f7f9fc",
        }
    ):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        ax = axes[0, 0]
        fine = np.linspace(0, 1, 2001)
        ax.plot(fine, problem.exact(fine), color="#183b56", label="Exact bump")
        ax.plot(
            final.nodes, final.values, ".", color="#ce5b28", markersize=3, label="Adaptive P1"
        )
        ax.set(title="Localized displacement", xlabel="Position x", ylabel="Displacement u")
        ax.legend()
        ax = axes[0, 1]
        ax.semilogy(
            (final.nodes[:-1] + final.nodes[1:]) / 2,
            np.diff(final.nodes),
            ".-",
            color="#087e8b",
            markersize=3,
        )
        ax.axvline(0.37, color="#ce5b28", linestyle="--", linewidth=1)
        ax.set(
            title=f"Where {final.elements} elements were spent",
            xlabel="Cell midpoint",
            ylabel="Cell length h",
        )
        ax = axes[1, 0]
        for data, label, color in [
            (adaptive, "Adaptive", "#087e8b"),
            (uniform, "Uniform: same cell count", "#ce5b28"),
        ]:
            ax.loglog(
                [r["free_dofs"] for r in data],
                [r["energy_error"] for r in data],
                "o-",
                label=label,
                color=color,
                markersize=4,
            )
        ax.set(
            title="Fair budget comparison",
            xlabel="Free degrees of freedom",
            ylabel="Energy error",
        )
        ax.legend()
        ax = axes[1, 1]
        for key, label, color in [
            ("energy_error", "Energy error", "#183b56"),
            ("l2_error", "L2 error", "#ce5b28"),
        ]:
            ax.loglog(
                [r["elements"] for r in smooth],
                [r[key] for r in smooth],
                "o-",
                label=label,
                color=color,
                markersize=4,
            )
        ax.set(
            title="Smooth sine: analytical convergence",
            xlabel="Uniform elements",
            ylabel="Error",
        )
        ax.legend()
        for ax in axes.flat:
            ax.grid(alpha=0.18, which="both")
        fig.suptitle(
            "MeshBeacon | Spend the mesh where the error lives",
            fontsize=17,
            fontweight="bold",
            color="#183b56",
        )
        fig.savefig(path, dpi=160)
        plt.close(fig)
