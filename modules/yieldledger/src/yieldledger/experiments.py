"""Deterministic synthetic experiments; all numerical claims originate here."""

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

from ._core import Material  # noqa: E402
from .bar import Bar, solve_step  # noqa: E402
from .ledger import ledger  # noqa: E402


def cyclic_path(n: int) -> np.ndarray:
    """0 -> 0.008 -> -0.008 -> 0.008 -> 0; every reversal is sampled."""
    vertices = [0.0, 0.008, -0.008, 0.008, 0.0]
    return np.r_[
        0,
        np.concatenate(
            [
                np.linspace(a, b, n + 1)[1:]
                for a, b in zip(vertices[:-1], vertices[1:], strict=True)
            ]
        ),
    ]


def layered_bar(subdivisions: int = 1) -> Bar:
    a = Material(200000, 200, 2000, 4000)
    b = Material(100000, 400, 5000, 5000)
    n = subdivisions
    return Bar(np.r_[np.full(n, 40 / n), np.full(n, 60 / n)], (a,) * n + (b,) * n, 10)


def run(output: str | Path) -> dict:
    """Write CSV, JSON, PNG and PDF; raise if numerical acceptance gates fail."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    material = Material(200000, 250, 2000, 3000)
    material_data = ledger(material, cyclic_path(400))
    np.savetxt(
        output / "cycle.csv",
        np.column_stack(list(material_data.values())),
        delimiter=",",
        header=",".join(material_data),
        comments="",
        fmt="%.16e",
    )
    balance = material_data["work"] - sum(
        material_data[k] for k in ("stored", "dissipation", "defect")
    )
    refinement = []
    for n in (25, 50, 100, 200, 400, 800):
        data = ledger(material, cyclic_path(n))
        refinement.append([n, data["defect"][-1], data["dissipation"][-1]])
    refinement = np.array(refinement)
    np.savetxt(
        output / "refinement.csv",
        refinement,
        delimiter=",",
        header="increments_per_leg,work_defect_MPa,dissipation_MPa",
        comments="",
        fmt="%.16e",
    )
    monotonic = np.linspace(0, 0.02, 501)
    mono = ledger(material, monotonic)
    exact = np.where(
        monotonic <= 0.00125,
        200000 * monotonic,
        250 + 200000 * 5000 / 205000 * (monotonic - 0.00125),
    )
    old = material.update(0.004).state
    tangent_errors = []
    for strain in (-0.012, -0.0005, 0.0001, 0.008):
        h = 1e-8
        fd = (
            material.update(strain + h, old).stress - material.update(strain - h, old).stress
        ) / (2 * h)
        tangent_errors.append(abs(fd - material.update(strain, old).tangent) / material.elastic)
    mesh_rows = []
    exact_stress = (0.4 + 40 * 200 / 6000) / (40 / 200000 + 60 / 100000 + 40 / 6000)
    for n in (1, 4, 16):
        model = layered_bar(n)
        result = solve_step(model, 0.4)
        x = np.r_[0, np.cumsum(model.lengths)]
        eps_a = exact_stress / 200000 + (exact_stress - 200) / 6000
        exact_u = np.minimum(x, 40) * eps_a + np.maximum(x - 40, 0) * exact_stress / 100000
        mesh_rows.append(
            [
                2 * n,
                np.max(np.abs(result.stress - exact_stress)),
                np.max(np.abs(result.displacement - exact_u)),
                result.iterations,
            ]
        )
    mesh_rows = np.array(mesh_rows)
    np.savetxt(
        output / "mesh.csv",
        mesh_rows,
        delimiter=",",
        header="elements,max_stress_error_MPa,max_nodal_error_mm,corrections",
        comments="",
        fmt="%.16e",
    )
    consistent = solve_step(layered_bar(), 0.4)
    elastic = solve_step(layered_bar(), 0.4, tangent="elastic", max_iterations=1000)
    niter = len(elastic.residual_history)
    np.savetxt(
        output / "newton.csv",
        np.column_stack([np.arange(niter), elastic.residual_history]),
        delimiter=",",
        header="correction,elastic_tangent_residual_N",
        comments="",
        fmt="%.16e",
    )
    metrics = {
        "material": {"E_MPa": 200000, "yield_MPa": 250, "Hi_MPa": 2000, "Hk_MPa": 3000},
        "cycle": {
            "increments": len(material_data["strain"]) - 1,
            "stored_MPa": float(material_data["stored"][-1]),
            "dissipation_MPa": float(material_data["dissipation"][-1]),
            "endpoint_work_MPa": float(material_data["work"][-1]),
            "work_defect_MPa": float(material_data["defect"][-1]),
            "max_balance_error_MPa": float(np.max(np.abs(balance))),
        },
        "monotonic_max_stress_error_MPa": float(np.max(np.abs(mono["stress"] - exact))),
        "tangent_max_error_normalized_by_E": float(max(tangent_errors)),
        "defect_refinement_orders": np.log2(refinement[:-1, 1] / refinement[1:, 1]).tolist(),
        "bar": {
            "exact_stress_MPa": exact_stress,
            "max_stress_error_MPa": float(np.max(mesh_rows[:, 1])),
            "max_nodal_error_mm": float(np.max(mesh_rows[:, 2])),
            "consistent_corrections": consistent.iterations,
            "elastic_corrections": elastic.iterations,
            "consistent_residuals_N": list(consistent.residual_history),
            "elastic_final_residual_N": elastic.residual_history[-1],
        },
    }
    gates = [
        metrics["monotonic_max_stress_error_MPa"] < 1e-8,
        metrics["tangent_max_error_normalized_by_E"] < 1e-8,
        metrics["cycle"]["max_balance_error_MPa"] < 1e-10,
        metrics["bar"]["max_stress_error_MPa"] < 1e-7,
        metrics["bar"]["max_nodal_error_mm"] < 1e-10,
        min(metrics["defect_refinement_orders"][-2:]) > 0.95,
        elastic.iterations > 10 * consistent.iterations,
    ]
    if not all(gates):
        raise RuntimeError(f"An experiment acceptance gate failed: {metrics}")
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    _plot(output, material_data, refinement, consistent, elastic)
    return metrics


def _plot(output, data, refinement, consistent, elastic):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.3), layout="constrained")
    ax = axes[0, 0]
    ax.plot(100 * data["strain"], data["stress"], color="#176b87", lw=1.8)
    ax.set(
        xlabel="Strain (%)", ylabel="Stress (MPa)", title="(a) Mixed-hardening cyclic response"
    )
    ax = axes[0, 1]
    for key, label, color in [
        ("stored", "Stored energy", "#176b87"),
        ("dissipation", "Plastic dissipation", "#b65c26"),
        ("work", "Endpoint work", "#333333"),
    ]:
        ax.plot(np.linspace(0, 4, len(data[key])), data[key], label=label, color=color)
    ax.set(
        xlabel="Loading leg coordinate",
        ylabel="Energy density (MPa)",
        title="(b) Cumulative energy ledger",
    )
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 0]
    ax.loglog(
        refinement[:, 0], refinement[:, 1], "o-", color="#b65c26", label="Measured defect"
    )
    ax.loglog(
        refinement[:, 0],
        refinement[-1, 1] * refinement[-1, 0] / refinement[:, 0],
        "--",
        color="gray",
        label="Slope -1",
    )
    ax.set(
        xlabel="Increments per loading leg",
        ylabel="Endpoint-work excess (MPa)",
        title="(c) Work quadrature refinement",
    )
    ax.set_xticks(refinement[:, 0], [str(int(n)) for n in refinement[:, 0]])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5], ["0.02", "0.05", "0.1", "0.2", "0.5"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 1]
    ax.semilogy(elastic.residual_history, label="Elastic tangent", color="#b65c26")
    ax.semilogy(
        np.maximum(consistent.residual_history, 1e-13),
        "o-",
        label="Consistent tangent",
        color="#176b87",
    )
    ax.set(
        xlabel="Newton corrections",
        ylabel="Free-force residual (N)",
        title="(d) Two-layer bar, U = 0.4 mm",
    )
    ax.legend(frameon=False, fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=0.18)
    fig.savefig(output / "overview.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(output / "overview.png", dpi=180)
    plt.close(fig)
