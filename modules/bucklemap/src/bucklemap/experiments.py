"""Deterministic synthetic experiments; CSV/JSON contain unrounded measurements."""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from .model import reconstruct, solve, spring_reference


def _field_error(result, z):
    xi, wi = leggauss(8)
    x = result.nodes
    h = np.diff(x)
    points = x[:-1, None] + h[:, None] * (xi + 1) / 2
    values = reconstruct(x, result.modes[:, 0], points)
    values /= reconstruct(x, result.modes[:, 0], 0.5)
    exact = (np.cos(z * (points - 0.5)) - np.cos(z / 2)) / (1 - np.cos(z / 2))
    return float(np.sqrt(np.sum((values - exact) ** 2 * h[:, None] * wi / 2)))


def measurements():
    convergence = []
    for kappa in (0.0, 10.0):
        exact = spring_reference(kappa)
        for n in (4, 8, 16, 32):
            result = solve(np.linspace(0, 1, n + 1), springs=(kappa, kappa))
            convergence.append(
                {
                    "kappa": kappa,
                    "elements": n,
                    "load": float(result.loads[0]),
                    "exact": exact,
                    "relative_load_error": float(result.loads[0] / exact - 1),
                    "field_l2_error": _field_error(result, np.sqrt(exact)),
                    "backward_error": float(result.residuals[0]),
                }
            )
    restraint = []
    for kappa in np.r_[0.0, np.logspace(-2, 4, 25)]:
        result = solve(np.linspace(0, 1, 33), springs=(kappa, kappa))
        exact = spring_reference(kappa)
        # Feynman-Hellmann derivative: q^T K' q, because q^T G q = 1.
        sensitivity = float(result.modes[1, 0] ** 2 + result.modes[-1, 0] ** 2)
        z = np.sqrt(exact)
        dz = -np.sin(z / 2) / ((kappa / 2 + 1) * np.cos(z / 2) - z / 2 * np.sin(z / 2))
        restraint.append(
            {
                "kappa": float(kappa),
                "load": float(result.loads[0]),
                "exact": exact,
                "sensitivity": sensitivity,
                "exact_sensitivity": float(2 * z * dz),
            }
        )
    boundary = []
    references = {
        "pinned": np.pi**2,
        "clamped": 4 * np.pi**2,
        "cantilever": np.pi**2 / 4,
        "clamped-pinned": brentq(lambda z: np.tan(z) - z, 4.4, 4.6) ** 2,
    }
    for name, exact in references.items():
        result = solve(np.linspace(0, 1, 33), support=name)
        boundary.append(
            {
                "support": name,
                "load": float(result.loads[0]),
                "exact": float(exact),
                "relative_error": float(result.loads[0] / exact - 1),
            }
        )
    orders = {}
    for kappa in (0.0, 10.0):
        rows = [r for r in convergence if r["kappa"] == kappa]
        orders[str(kappa)] = float(
            np.log2(rows[-2]["relative_load_error"] / rows[-1]["relative_load_error"])
        )
    return {
        "schema_version": 1,
        "setup": {
            "length": 1.0,
            "rigidity": 1.0,
            "field_error_quadrature_points_per_cell": 8,
            "restraint_elements": 32,
        },
        "convergence": convergence,
        "restraint": restraint,
        "boundary": boundary,
        "final_orders": orders,
    }


def run_study(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    data = measurements()
    (output / "metrics.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    for key in ("convergence", "restraint", "boundary"):
        with (output / f"{key}.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[key][0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(data[key])
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 180,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout="constrained")
    for kappa, color in [(0.0, "#087e8b"), (10.0, "#be4b28")]:
        rows = [r for r in data["convergence"] if r["kappa"] == kappa]
        n = [r["elements"] for r in rows]
        axes[0].loglog(
            n,
            [r["relative_load_error"] for r in rows],
            "o-",
            color=color,
            label=rf"$\kappa={kappa:g}$: load error",
        )
        axes[0].loglog(
            n,
            [r["backward_error"] for r in rows],
            "x--",
            color=color,
            label=rf"$\kappa={kappa:g}$: backward error",
        )
        axes[1].loglog(
            n,
            [r["field_l2_error"] for r in rows],
            "o-",
            color=color,
            label=rf"$\kappa={kappa:g}$",
        )
    axes[0].set_ylabel("Relative load / normalized residual")
    axes[1].set_ylabel("Displacement L2 error (unit peak)")
    for ax in axes:
        ax.set_xlabel("Number of elements")
        ax.set_xticks([4, 8, 16, 32], ["4", "8", "16", "32"])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.savefig(output / "convergence.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout="constrained")
    rows = data["restraint"][1:]
    kap = [r["kappa"] for r in rows]
    axes[0].semilogx(
        kap,
        [r["exact"] / np.pi**2 for r in rows],
        "-",
        color="#087e8b",
        label="Continuum scalar root",
    )
    axes[0].semilogx(
        kap[::3],
        [r["load"] / np.pi**2 for r in rows[::3]],
        "o",
        color="#be4b28",
        markersize=4,
        label="32-element FEM",
    )
    axes[0].axhline(4, ls=":", color="0.4", label="Clamped limit")
    axes[0].set_ylabel(r"Load relative to pinned: $P/\pi^2$")
    axes[1].loglog(
        kap,
        [r["exact_sensitivity"] for r in rows],
        "-",
        color="#087e8b",
        label="Continuum derivative",
    )
    axes[1].loglog(
        kap[::3],
        [r["sensitivity"] for r in rows[::3]],
        "o",
        color="#be4b28",
        markersize=4,
        label="Modal sensitivity",
    )
    axes[1].set_ylabel(r"Sensitivity $dP/d\kappa$")
    for ax in axes:
        ax.set_xlabel(r"Equal rotational restraint $\kappa=k_\theta L/EI$")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.savefig(output / "restraint.png")
    plt.close(fig)
    print(f"Wrote reproducible metrics, CSV tables and two figures to {output}")
