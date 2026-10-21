"""Verify that the delivered PDF's source and numerical inputs have not changed."""

import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
report = json.loads((root / "paper/validation.json").read_text())
for name, expected in report["inputs_sha256"].items():
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
assert (
    hashlib.sha256((root / "paper/paper.pdf").read_bytes()).hexdigest() == report["pdf_sha256"]
)
print(f"Verified {len(report['inputs_sha256'])} manuscript input hashes and PDF hash.")
