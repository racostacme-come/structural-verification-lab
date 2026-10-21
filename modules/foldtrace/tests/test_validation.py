import json

import numpy as np
import pytest

from foldtrace import Truss, trace
from foldtrace.cli import main
from foldtrace.validation import campaign, exact_folds, locate_folds


@pytest.mark.parametrize("ratio", [0.01, 0.2, 1, 10])
def test_bracketed_folds_agree_with_analytical_solution(ratio):
    model = Truss(ratio)
    np.testing.assert_allclose(
        locate_folds(model, trace(model)), exact_folds(ratio), atol=1e-10
    )


def test_campaign_outputs_and_order(tmp_path):
    summary = campaign(tmp_path)
    assert summary["interpolation_orders"][-1] > 1.8
    assert summary["max_independent_force_error"] < 1e-11
    assert summary["max_fold_location_error"] < 1e-10
    assert summary == json.loads((tmp_path / "validation.json").read_text())
    assert (tmp_path / "equilibrium_audit.png").read_bytes().startswith(b"\x89PNG")
    assert (tmp_path / "path.csv").read_text().splitlines()[0].startswith("q,load,tangent")


def test_cli_outputs(tmp_path, capsys):
    assert main(["--output", str(tmp_path)]) == 0
    assert json.loads(capsys.readouterr().out)["stop_reason"] == "target_reached"


def test_cli_invalid_ratio(tmp_path, capsys):
    with pytest.raises(SystemExit) as caught:
        main(["--ratio", "0", "--output", str(tmp_path)])
    assert caught.value.code == 2
    assert "height/span ratio" in capsys.readouterr().err
