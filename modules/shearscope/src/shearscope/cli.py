"""Command line analysis and deterministic validation campaign."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .model import Beam, exact_tip, modes, solve_static


def campaign(output: Path):
    """Generate original numerical data and a three-panel validation figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator

    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for slenderness in (5, 10, 20, 50, 100, 200, 500, 1000):
        beam = Beam(height=1 / slenderness)
        for n in (4, 8, 16, 32, 64):
            for integration in ("full", "reduced"):
                result = solve_static(beam, n, integration)
                rows.append(
                    {
                        "slenderness": slenderness,
                        "elements": n,
                        "integration": integration,
                        "normalized_tip": result.displacement[-1] / exact_tip(beam),
                        "relative_residual": result.relative_residual,
                    }
                )
    with (output / "locking.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    beam = Beam(height=0.005)
    roots = np.array([1.875104068711961, 4.694091132974174, 7.854757438237613])
    reference = roots**2 * np.sqrt(beam.ei / (beam.density * beam.area)) / (2 * np.pi)
    modal_rows = []
    for n in (8, 16, 32, 64, 80):
        for integration in ("full", "reduced"):
            frequencies, _ = modes(beam, n, integration, count=3)
            for j, frequency in enumerate(frequencies):
                modal_rows.append(
                    {
                        "elements": n,
                        "integration": integration,
                        "mode": j + 1,
                        "frequency_hz": frequency,
                        "eb_limit_hz": reference[j],
                        "relative_to_eb": frequency / reference[j],
                    }
                )
    with (output / "modes.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(modal_rows[0]))
        writer.writeheader()
        writer.writerows(modal_rows)

    colors = {"full": "#c45a30", "reduced": "#007d87"}
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.9), layout="constrained")
    fig.set_facecolor("#fafaf7")
    fig.suptitle("ShearScope  |  When a beam mesh lies", fontsize=21, fontweight="bold")
    for integration in ("full", "reduced"):
        selected = [r for r in rows if r["elements"] == 8 and r["integration"] == integration]
        axes[0].semilogx(
            [r["slenderness"] for r in selected],
            [r["normalized_tip"] for r in selected],
            "o-",
            color=colors[integration],
            label=integration.capitalize(),
        )
        selected = [
            r for r in rows if r["slenderness"] == 100 and r["integration"] == integration
        ]
        axes[1].loglog(
            [r["elements"] for r in selected],
            [abs(1 - r["normalized_tip"]) for r in selected],
            "o-",
            color=colors[integration],
            label=integration.capitalize(),
        )
        selected = [r for r in modal_rows if r["mode"] == 1 and r["integration"] == integration]
        axes[2].semilogy(
            [r["elements"] for r in selected],
            [r["relative_to_eb"] for r in selected],
            "o-",
            color=colors[integration],
            label=integration.capitalize(),
        )
    axes[0].axhline(1, color="#555555", linestyle=":", label="Exact Timoshenko")
    axes[0].set(
        xlabel="Slenderness L/h",
        ylabel="Computed / exact tip deflection",
        title="01  Static compliance · 8 elements",
        ylim=(-0.04, 1.08),
    )
    axes[1].set(
        xlabel="Number of elements",
        ylabel="Relative tip error",
        title="02  Mesh convergence · L/h = 100",
    )
    axes[1].set_xticks([4, 8, 16, 32, 64], labels=["4", "8", "16", "32", "64"])
    axes[1].xaxis.set_minor_locator(NullLocator())
    axes[2].axhline(1, color="#555555", linestyle=":", label="Euler–Bernoulli limit")
    axes[2].set(
        xlabel="Number of elements",
        ylabel="First frequency / EB limit",
        title="03  Modal bias · L/h = 200",
    )
    for ax in axes:
        ax.grid(alpha=0.2, which="both")
        ax.legend(fontsize=8, frameon=False)
    fig.savefig(output / "shearscope.png", dpi=180, metadata={"Software": "ShearScope"})
    plt.close(fig)

    full = next(
        r["normalized_tip"]
        for r in rows
        if r["slenderness"] == 100 and r["elements"] == 8 and r["integration"] == "full"
    )
    reduced = next(
        r["normalized_tip"]
        for r in rows
        if r["slenderness"] == 100 and r["elements"] == 8 and r["integration"] == "reduced"
    )
    summary = {
        "static_cases": len(rows),
        "modal_cases": len(modal_rows) // 3,
        "L_over_h_100_n8_full_tip_ratio": full,
        "L_over_h_100_n8_reduced_tip_ratio": reduced,
        "max_static_relative_residual": max(r["relative_residual"] for r in rows),
        "note": "Modal reference is the slender Euler-Bernoulli limit, not exact Timoshenko.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    analysis = sub.add_parser("analyze", help="Solve a rectangular cantilever in SI units")
    analysis.add_argument("--elements", type=int, default=32)
    analysis.add_argument("--length", type=float, default=1.0)
    analysis.add_argument("--height", type=float, default=0.01)
    analysis.add_argument("--force", type=float, default=1.0)
    analysis.add_argument("--integration", choices=("full", "reduced"), default="reduced")
    study = sub.add_parser("campaign", help="Write validation CSV, JSON and PNG files")
    study.add_argument("--output", type=Path, default=Path("out"))
    args = parser.parse_args(argv)
    try:
        if args.command == "campaign":
            result = campaign(args.output)
        else:
            beam = Beam(length=args.length, height=args.height)
            solution = solve_static(beam, args.elements, args.integration, args.force)
            result = {
                "tip_m": float(solution.displacement[-1]),
                "exact_tip_m": exact_tip(beam, args.force),
                "root_reactions_N_Nm": solution.reactions.tolist(),
                "relative_residual": solution.relative_residual,
            }
    except (ValueError, ArithmeticError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
