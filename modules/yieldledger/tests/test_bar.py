import numpy as np
import pytest

from yieldledger import Bar, ConvergenceError, Material, solve_step


def bar_with_subdivisions(n):
    a, b = Material(200000, 200, 2000, 4000), Material(100000, 400, 5000, 5000)
    return Bar(np.r_[np.full(n, 40 / n), np.full(n, 60 / n)], (a,) * n + (b,) * n, 10)


@pytest.mark.parametrize("n", [1, 4, 16])
def test_two_layer_continuum_solution(n):
    bar = bar_with_subdivisions(n)
    u = 0.4
    stress = (u + 40 * 200 / 6000) / (40 / 200000 + 60 / 100000 + 40 / 6000)
    result = solve_step(bar, u)
    x = np.r_[0, np.cumsum(bar.lengths)]
    eps_a = stress / 200000 + (stress - 200) / 6000
    exact_u = np.minimum(x, 40) * eps_a + np.maximum(x - 40, 0) * stress / 100000
    np.testing.assert_allclose(result.stress, stress, atol=2e-8, rtol=0)
    np.testing.assert_allclose(result.displacement, exact_u, atol=1e-11, rtol=0)


def test_consistent_and_elastic_tangent_converge_to_same_solution():
    bar = bar_with_subdivisions(1)
    exact = solve_step(bar, 0.4)
    slow = solve_step(bar, 0.4, tangent="elastic", max_iterations=1000)
    np.testing.assert_allclose(slow.displacement, exact.displacement, atol=1e-10)
    assert slow.iterations > 10 * exact.iterations


def test_failure_leaves_old_state_intact():
    bar = bar_with_subdivisions(1)
    old = solve_step(bar, 0.3)
    u = old.displacement.copy()
    state = [(s.plastic_strain, s.accumulated) for s in old.states]
    with pytest.raises(ConvergenceError):
        solve_step(bar, -0.6, previous=old, max_iterations=0)
    np.testing.assert_array_equal(old.displacement, u)
    assert [(s.plastic_strain, s.accumulated) for s in old.states] == state


def test_cyclic_bar_equilibrium_and_global_work_identity():
    bar = bar_with_subdivisions(4)
    history = np.r_[np.linspace(0, 0.6, 61), np.linspace(0.6, -0.6, 121)[1:]]
    old = solve_step(bar, 0)
    for value in history[1:]:
        new = solve_step(bar, value, old)
        old_e = np.diff(old.displacement) / bar.lengths
        new_e = np.diff(new.displacement) / bar.lengths
        stored = dissipation = defect = 0.0
        for m, length, eo, en, so, sn in zip(
            bar.materials, bar.lengths, old_e, new_e, old.states, new.states, strict=True
        ):
            dep = sn.plastic_strain - so.plastic_strain
            da = sn.accumulated - so.accumulated
            volume = bar.area * length
            stored += volume * (m.energy(en, sn) - m.energy(eo, so))
            dissipation += volume * m.yield_stress * da
            defect += (
                volume
                * (
                    m.elastic * (en - eo - dep) ** 2
                    + m.isotropic * da**2
                    + m.kinematic * dep**2
                )
                / 2
            )
        work = new.stress[-1] * bar.area * (value - old.displacement[-1])
        assert work == pytest.approx(stored + dissipation + defect, abs=1e-8)
        assert np.ptp(new.stress) < 1e-7
        old = new


def test_one_element_has_no_free_dofs():
    m = Material(200000, 200, 1000)
    result = solve_step(Bar([10], (m,)), 0.03)
    assert result.iterations == 0
    assert result.stress[0] == pytest.approx(m.update(0.003).stress)


@pytest.mark.parametrize("lengths,area", [([], 1), ([0], 1), ([1], 0), ([np.nan], 1)])
def test_invalid_geometry(lengths, area):
    with pytest.raises(ValueError):
        Bar(lengths, (Material(1000, 10, 100),), area)


def test_perfect_plastic_fem_rejected():
    with pytest.raises(ValueError):
        Bar([1], (Material(1000, 10),))


@pytest.mark.parametrize(
    "settings", [{"tangent": "bad"}, {"max_iterations": -1}, {"atol": 0}, {"rtol": np.nan}]
)
def test_invalid_solver_settings(settings):
    with pytest.raises(ValueError):
        solve_step(bar_with_subdivisions(1), 0.1, **settings)
