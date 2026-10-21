"""Recompute numerical records; tolerate documented BLAS roundoff, not changed physics."""

import csv
import json
from pathlib import Path

import numpy as np

from bucklemap.experiments import measurements


def main():
    folder = Path(__file__).resolve().parents[1] / "examples/results"
    saved = json.loads((folder / "metrics.json").read_text())
    fresh = measurements()
    assert saved["setup"] == fresh["setup"]
    for group in ("convergence", "restraint", "boundary"):
        assert len(saved[group]) == len(fresh[group])
        for old, new in zip(saved[group], fresh[group], strict=True):
            assert old.keys() == new.keys()
            for key in old:
                if isinstance(old[key], str):
                    assert old[key] == new[key]
                elif key == "backward_error":
                    assert old[key] < 1e-12 and new[key] < 1e-12
                else:
                    np.testing.assert_allclose(
                        new[key], old[key], rtol=3e-8, atol=3e-10, err_msg=f"{group}: {key}"
                    )
        with (folder / f"{group}.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == len(saved[group])
        for csv_row, json_row in zip(rows, saved[group], strict=True):
            for key, value in json_row.items():
                assert csv_row[key] == str(value)
    for key in saved["final_orders"]:
        assert abs(saved["final_orders"][key] - fresh["final_orders"][key]) < 0.02
    print("All numerical records reproduced; CSV/JSON agree; backward errors below 1e-12.")


if __name__ == "__main__":
    main()
