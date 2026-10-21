"""Deterministic numerical experiments and original figures."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from .fem import assemble, critical_penalty, errors, solve, unit_square
from .problems import exact, gradient, source


def run():
    """Return JSON-serializable measurements; no timings or random data."""
    convergence = []
    thresholds = []
    for n in (4, 8, 16, 32):
        mesh = unit_square(n)
        system = assemble(mesh, source, exact)
        thresholds.append({"n": n, "gamma_star": critical_penalty(system)})
        for method in ("strong", "penalty", "nitsche"):
            solution = solve(system, method, 10)
            convergence.append(
                {
                    "n": n,
                    "dofs": len(mesh.points),
                    "method": method,
                    "gamma": 10,
                    **errors(mesh, solution.values, exact, gradient),
                    "backward_error": solution.backward_error,
                }
            )
    system = assemble(unit_square(8), source, exact)
    threshold = critical_penalty(system)
    sweep = []
    for gamma in (1, 1.4, 1.55, 1.57, 2, 3, 5, 10, 30, 100, 1000, 10000, 100000):
        for method in ("penalty", "nitsche"):
            eigenvalues = eigvalsh(system.matrix(method, gamma).toarray())
            stable = eigenvalues[0] > 0
            record = {
                "gamma": gamma,
                "method": method,
                "lambda_min": float(eigenvalues[0]),
                "condition_2": float(eigenvalues[-1] / eigenvalues[0]) if stable else None,
            }
            if stable:
                result = solve(system, method, gamma)
                record.update(errors(system.mesh, result.values, exact, gradient))
                record["backward_error"] = result.backward_error
            sweep.append(record)

    def affine(p):
        return 0.7 + 1.3 * p[:, 0] - 2.1 * p[:, 1]

    patch_system = assemble(unit_square(8), lambda p: np.zeros(len(p)), affine)
    patch = {}
    for method in ("strong", "penalty", "nitsche"):
        result = solve(patch_system, method, 10)
        patch[method] = float(np.max(abs(result.values - affine(patch_system.mesh.points))))
    low = assemble(unit_square(4), source, exact, order=6)
    high = assemble(unit_square(4), source, exact, order=10)
    value_low, value_high = solve(low).values, solve(high).values
    e8 = errors(low.mesh, value_low, exact, gradient, order=8)
    e12 = errors(low.mesh, value_low, exact, gradient, order=12)
    orders = {}
    for method in ("strong", "penalty", "nitsche"):
        last = [r for r in convergence if r["method"] == method][-2:]
        orders[method] = {
            key: float(np.log2(last[0][key] / last[1][key])) for key in ("l2", "h1")
        }
    return {
        "schema_version": 1,
        "convergence": convergence,
        "thresholds": thresholds,
        "sweep": sweep,
        "gamma_star_n8": threshold,
        "patch_max_error": patch,
        "orders": orders,
        "quadrature_check": {
            "load_max_difference": float(np.max(abs(low.load - high.load))),
            "solution_max_difference": float(np.max(abs(value_low - value_high))),
            "norm_max_difference": max(abs(e8[k] - e12[k]) for k in e8),
        },
    }


def write_results(output, data, figures=True):
    """Write measured JSON/CSV and, optionally, publication figures."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "metrics.json").write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )
    for key in ("convergence", "sweep", "thresholds"):
        rows = data[key]
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (output / f"{key}.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if figures:
        plot_results(output, data)


def plot_results(output, data):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, ScalarFormatter

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "savefig.dpi": 180,
        }
    )
    colors = {"strong": "#596773", "nitsche": "#087f8c", "penalty": "#cc6533"}
    labels = {"strong": "Strong", "nitsche": "Nitsche", "penalty": "Penalty"}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), layout="constrained")
    for method in colors:
        rows = [r for r in data["convergence"] if r["method"] == method]
        for ax, key, title in zip(
            axes, ("l2", "h1"), (r"Field $L^2$ error", r"Gradient $H^1$ error"), strict=True
        ):
            ax.loglog(
                [r["n"] for r in rows],
                [r[key] for r in rows],
                "o-",
                color=colors[method],
                label=labels[method],
                markersize=4,
            )
            ax.set(xlabel="Subdivisions per side, n", ylabel=title)
            ax.set_xticks([4, 8, 16, 32])
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.set_minor_locator(NullLocator())
    axes[0].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"convergence.{suffix}")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.0), layout="constrained")
    for method in ("nitsche", "penalty"):
        rows = [
            r for r in data["sweep"] if r["method"] == method and r["condition_2"] is not None
        ]
        for ax, key, label in zip(
            axes,
            ("boundary_l2", "l2", "condition_2"),
            (r"Boundary $L^2$ error", r"Field $L^2$ error", r"Matrix $\kappa_2$"),
            strict=True,
        ):
            ax.loglog(
                [r["gamma"] for r in rows],
                [r[key] for r in rows],
                "o-",
                label=labels[method],
                color=colors[method],
                markersize=3,
            )
            ax.axvline(data["gamma_star_n8"], color="#999999", ls=":", lw=1)
            ax.set(xlabel=r"Penalty parameter $\gamma$", ylabel=label)
            ax.set_xticks([1, 100, 10000])
    axes[0].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"parameter.{suffix}")
    plt.close(fig)
