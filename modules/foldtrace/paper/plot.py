"""Regenerate the manuscript figure from committed numerical CSV outputs."""

import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parent.parent


def read(name):
    with (ROOT / "results" / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def values(rows, key):
    if not rows:
        raise ValueError(f"No recorded rows selected for {key}")
    return np.array([float(row[key]) for row in rows])


plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(6.35, 2.35), layout="constrained")
rows = read("path.csv")
axes[0].plot(values(rows, "q"), values(rows, "load"), "o-", markersize=3)
axes[0].axhline(0, color="0.7", linewidth=0.7)
axes[0].set(title="Accepted equilibrium path", xlabel="Displacement q", ylabel="Load lambda")
rows = read("convergence.csv")
axes[1].loglog(values(rows, "step"), values(rows, "max_midpoint_error"), "o-")
axes[1].set(
    title="Path interpolation", xlabel="Predictor step", ylabel="Maximum midpoint error"
)

for ax in axes:
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    if ax.get_xscale() == "log":
        points = np.unique(ax.lines[0].get_xdata())
        if 2 <= len(points) <= 6:
            ax.set_xticks(points, labels=[f"{x:.3g}" for x in points])
    ax.grid(alpha=0.2, which="both")
fig.savefig(ROOT / "paper" / "figure.pdf")
plt.close(fig)
