import csv
import json

from meshbeacon.cli import main
from meshbeacon.report import benchmark


def test_reproducible_report(tmp_path):
    a = benchmark(tmp_path / "first")
    b = benchmark(tmp_path / "second")
    assert a == b
    assert all(a["checks"].values())
    for name in ("convergence.csv", "adaptive_nodes.csv", "validation.json"):
        assert (tmp_path / "first" / name).read_bytes() == (
            tmp_path / "second" / name
        ).read_bytes()
    assert (tmp_path / "first" / "mesh_audit.png").read_bytes().startswith(b"\x89PNG")
    rows = list(csv.DictReader((tmp_path / "first" / "convergence.csv").open()))
    assert len(rows) == 6 + 2 * a["adaptive_solves"]
    assert json.loads((tmp_path / "first" / "validation.json").read_text()) == a


def test_cli_failure(tmp_path, capsys):
    path = tmp_path / "file"
    path.write_text("existing file")
    assert main(["--output", str(path)]) == 1
    assert "meshbeacon:" in capsys.readouterr().err
    assert path.read_text() == "existing file"
