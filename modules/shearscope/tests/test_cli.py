import json
import subprocess
import sys

import pytest

from shearscope.cli import campaign, main


def test_cli_analysis(capsys):
    main(["analyze", "--force", "-2"])
    result = json.loads(capsys.readouterr().out)
    assert result["tip_m"] < 0
    assert result["root_reactions_N_Nm"] == pytest.approx([2, 2], abs=1e-7)


def test_cli_validation(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["analyze", "--height", "0"])
    assert exc.value.code == 2
    assert "positive" in capsys.readouterr().err


def test_module_entrypoint():
    run = subprocess.run(
        [sys.executable, "-m", "shearscope.cli", "analyze"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(run.stdout)["tip_m"] > 0


def test_campaign_outputs(tmp_path):
    summary = campaign(tmp_path)
    assert summary["static_cases"] == 80
    assert summary["modal_cases"] == 10
    assert summary["L_over_h_100_n8_reduced_tip_ratio"] > 0.996
    assert summary["max_static_relative_residual"] < 1e-4
    assert (tmp_path / "shearscope.png").stat().st_size > 10000
    assert len((tmp_path / "locking.csv").read_text().splitlines()) == 81
