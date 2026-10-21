import math

import numpy as np
import pytest

from foldtrace import ContinuationError, Truss, trace


@pytest.mark.parametrize("ratio", [0.01, 0.1, 0.2, 0.5, 1, 3, 10])
def test_path_crosses_both_folds(ratio):
    model = Truss(ratio)
    path = trace(model)
    assert path.stop_reason == "target_reached"
    assert path.points[-1].q >= 2.2
    assert np.all(np.diff([p.q for p in path.points]) > 0)
    signs = np.sign([p.tangent for p in path.points])
    assert np.count_nonzero(np.diff(signs)) == 2
    assert max(abs(p.residual) for p in path.points) <= 1e-11
    assert max(abs(p.constraint) for p in path.points) <= 1e-11
    for point in path.points:
        assert abs(model.force(point.q) - point.load) / path.load_scale <= 1e-11


def test_rejections_reduce_step_and_preserve_equilibrium():
    path = trace(Truss(10), step=0.5, max_step=0.5, load_scale=0.01)
    assert path.rejected_steps > 0
    assert sum(p.retries for p in path.points) == path.rejected_steps
    assert any(p.step < 0.5 for p in path.points[1:])


def test_budget_failure_retains_only_accepted_points():
    with pytest.raises(ContinuationError, match="step_budget_exhausted") as caught:
        trace(Truss(), max_steps=1)
    assert len(caught.value.path.points) == 2
    assert caught.value.path.points[-1].residual < 1e-11


def test_minimum_step_failure():
    with pytest.raises(ContinuationError, match="minimum_step_reached") as caught:
        trace(Truss(), step=0.5, min_step=0.5, max_step=0.5, max_iterations=1)
    assert len(caught.value.path.points) == 1
    assert caught.value.path.rejected_steps == 1


@pytest.mark.parametrize(
    "kwargs",
    [
        {"q_stop": 0},
        {"q_stop": 3.1},
        {"step": math.nan},
        {"step": 0},
        {"max_step": 0.6},
        {"min_step": 0.2},
        {"load_scale": 0},
        {"load_scale": math.inf},
        {"tolerance": 1e-16},
        {"max_steps": 0},
        {"max_steps": True},
        {"max_iterations": 1.2},
    ],
)
def test_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        trace(Truss(), **kwargs)


def test_determinism_and_immutable_history():
    one = trace(Truss())
    assert one == trace(Truss())
    with pytest.raises(AttributeError):
        one.points[-1].q = 1
