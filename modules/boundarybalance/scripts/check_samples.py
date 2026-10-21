"""Recompute numerical outputs and compare with committed measurements."""

import json
from pathlib import Path

import numpy as np

from boundarybalance.experiments import run


def compare(actual, expected, location="root"):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), location
        for key in expected:
            compare(actual[key], expected[key], f"{location}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), location
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, f"{location}[{i}]")
    elif isinstance(expected, float):
        # Roundoff-level diagnostics vary by BLAS and sparse factorization.
        np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=3e-11, err_msg=location)
    else:
        assert actual == expected, location


if __name__ == "__main__":
    path = Path(__file__).resolve().parents[1] / "examples/results/metrics.json"
    compare(run(), json.loads(path.read_text()))
    print("All recorded numerical measurements reproduced within stated tolerances.")
