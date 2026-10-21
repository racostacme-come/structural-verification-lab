import numpy as np
import pytest

from yieldledger import Material, State, ledger


@pytest.fixture
def material():
    return Material(200000, 250, 2000, 3000)


def test_monotonic_bilinear_analytical_solution(material):
    strains = np.linspace(0, 0.02, 501)
    result = ledger(material, strains)
    ey = 250 / 200000
    et = 200000 * 5000 / 205000
    exact = np.where(strains <= ey, 200000 * strains, 250 + et * (strains - ey))
    np.testing.assert_allclose(result["stress"], exact, rtol=2e-14, atol=2e-12)


@pytest.mark.parametrize("strain", [-0.012, -0.0005, 0.0001, 0.008])
def test_consistent_tangent_by_finite_difference(material, strain):
    old = material.update(0.004).state
    h = 1e-8
    fd = (material.update(strain + h, old).stress - material.update(strain - h, old).stress) / (
        2 * h
    )
    assert material.update(strain, old).tangent == pytest.approx(fd, rel=2e-8)


def test_reverse_yield_and_elastic_unloading(material):
    peak = material.update(0.006)
    state = peak.state
    reverse_yield = material.kinematic * state.plastic_strain - (
        material.yield_stress + material.isotropic * state.accumulated
    )
    reverse_strain = state.plastic_strain + reverse_yield / material.elastic
    inside = material.update(reverse_strain + 1e-7, state)
    outside = material.update(reverse_strain - 1e-7, state)
    assert inside.plastic_increment == 0
    assert inside.tangent == material.elastic
    assert outside.plastic_increment > 0
    assert outside.state.plastic_strain < state.plastic_strain
    assert outside.state.accumulated > state.accumulated


def test_random_cyclic_complementarity_and_energy(material):
    rng = np.random.default_rng(42)
    path = np.r_[0, np.cumsum(rng.normal(0, 0.0004, 1000))]
    result = ledger(material, path)
    f = np.abs(result["stress"] - material.kinematic * result["plastic_strain"]) - (
        material.yield_stress + material.isotropic * result["alpha"]
    )
    da = np.diff(result["alpha"], prepend=0)
    assert np.max(f) < 1e-9
    assert np.min(da) >= 0
    assert np.max(np.abs(f * da)) < 1e-10
    np.testing.assert_allclose(
        result["work"], result["stored"] + result["dissipation"] + result["defect"], atol=2e-12
    )
    assert np.all(np.diff(result["defect"]) >= 0)


def test_perfect_plasticity_and_yield_kink():
    m = Material(1000, 10)
    at_yield = m.update(0.01)
    assert at_yield.tangent == 1000  # documented elastic-sided choice at f_trial == 0
    assert m.update(0.1).stress == pytest.approx(10)
    assert m.update(0.1).tangent == 0


def test_endpoint_work_defect_refines(material):
    defects = []
    for n in (100, 200, 400):
        path = np.r_[np.linspace(0, 0.008, n + 1), np.linspace(0.008, -0.008, 2 * n + 1)[1:]]
        data = ledger(material, path)
        defects.append(data["defect"][-1])
    assert 1.95 < defects[0] / defects[1] < 2.05
    assert 1.95 < defects[1] / defects[2] < 2.05


@pytest.mark.parametrize(
    "args", [(0, 1, 0, 0), (1, -1, 0, 0), (1, 1, -1, 0), (np.inf, 1, 0, 0), (1, 1, 0, np.nan)]
)
def test_bad_material_rejected(args):
    with pytest.raises(ValueError):
        Material(*args)


@pytest.mark.parametrize("state", [State(1, 0), State(0, -1), State(np.nan, 1)])
def test_bad_state_rejected(material, state):
    with pytest.raises(ValueError):
        material.update(0, state)


@pytest.mark.parametrize("path", [[], [1], [0, np.nan], [[0, 1]]])
def test_bad_path_rejected(material, path):
    with pytest.raises(ValueError):
        ledger(material, path)


def test_overflow_rejected(material):
    with pytest.raises(OverflowError):
        material.update(1e308)


def test_trial_does_not_mutate_committed_state(material):
    old = material.update(0.005).state
    snapshot = old.plastic_strain, old.accumulated
    material.update(-0.01, old)
    assert (old.plastic_strain, old.accumulated) == snapshot
