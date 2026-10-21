import json
import subprocess
import sys

import pytest

from boundarybalance.cli import main


def test_module_cli(tmp_path):
    path = tmp_path / "result.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "boundarybalance.cli",
            "solve",
            "--n",
            "4",
            "--output",
            str(path),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    data = json.loads(result.stdout)
    assert data == json.loads(path.read_text())
    assert data["l2"] < 0.1
    assert data["backward_error"] < 1e-14


@pytest.mark.parametrize("arguments", [["--n", "1"], ["--gamma", "1"], ["--gamma", "nan"]])
def test_cli_rejects_invalid_input(arguments, capsys):
    with pytest.raises(SystemExit) as error:
        main(["solve", *arguments])
    assert error.value.code == 2
    assert "error:" in capsys.readouterr().err
