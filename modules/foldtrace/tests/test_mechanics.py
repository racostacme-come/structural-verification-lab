import math

import numpy as np
import pytest
from foldtrace._core import correct

from foldtrace import Truss


@pytest.mark.parametrize("ratio", [0.01, 0.1, 0.2, 1.0, 10.0])
def test_geometry_force_tangent_and_energy(ratio):
    model = Truss(ratio)
    l0 = math.sqrt(1 + ratio**2)
    for q in np.linspace(-0.2, 2.2, 23):
        length = math.sqrt(1 + ratio**2 * (1 - q) ** 2)
        # Independent engineering-strain axial force projected onto vertical axis.
        expected = -(length - l0) * (1 - q) / length
        assert model.force(q) == pytest.approx(expected, abs=2e-14)
        eps = 1e-6
        numerical = (model.force(q + eps) - model.force(q - eps)) / (2 * eps)
        assert model.tangent(q) == pytest.approx(numerical, rel=2e-6, abs=2e-9)
        energy_derivative = (model.energy(q + eps) - model.energy(q - eps)) / (2 * eps)
        assert energy_derivative == pytest.approx(
            2 * ratio**2 / l0 * model.force(q), rel=2e-6, abs=2e-8
        )


@pytest.mark.parametrize("ratio", [0.01, 0.2, 1, 10])
def test_exact_folds_and_symmetry(ratio):
    model = Truss(ratio)
    offset = math.sqrt(math.expm1(math.log1p(ratio**2) / 3) / ratio**2)
    for q in [1 - offset, 1 + offset]:
        assert abs(model.tangent(q)) < 2e-14
    for q in [0.0, 0.13, 0.8, 1.0]:
        assert model.force(q) == pytest.approx(-model.force(2 - q), abs=1e-15)
        assert model.energy(q) == pytest.approx(model.energy(2 - q), abs=1e-15)
    assert model.force(0) == model.force(1) == model.force(2) == 0
    assert model.energy(0) == model.energy(2) == 0
    assert model.energy(1) > 0


@pytest.mark.parametrize("ratio", [0, -1, 0.009, 10.1, math.nan, math.inf])
def test_invalid_geometry(ratio):
    with pytest.raises(ValueError):
        Truss(ratio)


@pytest.mark.parametrize("q", [math.nan, math.inf, -10.1, 10.1])
def test_invalid_state(q):
    for method in [Truss().force, Truss().tangent, Truss().energy]:
        with pytest.raises(ValueError):
            method(q)


def test_corrector_at_fold_is_regular():
    model = Truss()
    q = 1 - math.sqrt(math.expm1(math.log1p(0.04) / 3) / 0.04)
    result = correct(model, q, model.force(q) + 0.001, 1.0, 0.0, 0.04, 1e-12, 5)
    assert result.converged
    assert result.q == q
    assert result.load == pytest.approx(model.force(q), abs=1e-14)


@pytest.mark.parametrize("field,value", [(2, 2.0), (4, 0.0), (5, -1.0), (6, 0)])
def test_corrector_rejects_bad_settings(field, value):
    args = [0.1, 0.001, 1.0, 0.0, 0.04, 1e-12, 10]
    args[field] = value
    with pytest.raises(ValueError):
        correct(Truss(), *args)
