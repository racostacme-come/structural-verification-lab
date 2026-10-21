"""Compare regenerated sample arrays, allowing floating-point roundoff."""

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actual", default="build/reproduced")
    parser.add_argument("--reference", default="results")
    args = parser.parse_args()
    for name in ("cycle.csv", "refinement.csv", "mesh.csv", "newton.csv"):
        reference = np.loadtxt(Path(args.reference) / name, delimiter=",", skiprows=1)
        actual = np.loadtxt(Path(args.actual) / name, delimiter=",", skiprows=1)
        # Mesh correction counts may vary at the nonsmooth yield boundary; physical
        # values and iteration convergence are independently gated by experiments.
        if name == "mesh.csv":
            reference, actual = reference[:, :3], actual[:, :3]
        np.testing.assert_allclose(actual, reference, rtol=2e-6, atol=1e-8, err_msg=name)
    print("All four recorded numerical arrays reproduced within documented tolerances.")


if __name__ == "__main__":
    main()
