import numpy as np
import pytest
from scipy.optimize import brentq

from bucklemap import assemble, solve, spring_reference


@pytest.mark.parametrize(
    "support,exact",
    [
        ("pinned", np.pi**2),
        ("clamped", 4 * np.pi**2),
        ("cantilever", np.pi**2 / 4),
        ("clamped-pinned", brentq(lambda z: np.tan(z) - z, 4.4, 4.6) ** 2),
    ],
)
def test_euler_boundary_conditions(support, exact):
    result = solve(np.linspace(0, 1, 33), support=support)
    assert 0 < result.loads[0] / exact - 1 < 3e-6
    assert result.residuals[0] < 1e-12


def test_higher_modes_and_geometric_orthogonality():
    x = np.linspace(0, 1, 33)
    result = solve(x, modes=4)
    np.testing.assert_allclose(result.loads, np.pi**2 * np.arange(1, 5) ** 2, rtol=4e-5)
    _, g = assemble(x)
    np.testing.assert_allclose(result.modes.T @ g @ result.modes, np.eye(4), atol=2e-12)
    assert np.all(result.modes[[0, -2]] == 0)


def test_nested_ritz_upper_bounds_and_fourth_order():
    loads = np.array([solve(np.linspace(0, 1, n + 1)).loads[0] for n in (4, 8, 16, 32)])
    errors = loads / np.pi**2 - 1
    assert np.all(errors > 0) and np.all(np.diff(loads) < 0)
    assert 3.9 < np.log2(errors[-2] / errors[-1]) < 4.1


@pytest.mark.parametrize("kappa", [0, 0.1, 1, 10, 100, 10000])
def test_spring_continuum_equation(kappa):
    value = solve(np.linspace(0, 1, 33), springs=(kappa, kappa)).loads[0]
    assert abs(value / spring_reference(kappa) - 1) < 3e-6


def test_dimensional_scaling():
    base = solve(np.linspace(0, 1, 17), springs=(3, 5), modes=3)
    length, rigidity = 2.7, 19.0
    scaled = solve(
        np.linspace(5, 5 + length, 17),
        rigidity,
        springs=(3 * rigidity / length, 5 * rigidity / length),
        modes=3,
    )
    np.testing.assert_allclose(scaled.loads, base.loads * rigidity / length**2, rtol=2e-11)


def test_nonuniform_heterogeneous_refinement_and_spring_monotonicity():
    x = np.array([0, 0.1, 0.35, 0.5, 0.8, 1])
    ei = np.array([1, 1, 3, 3, 2])
    fine = np.sort(np.r_[x, (x[:-1] + x[1:]) / 2])
    coarse_load = solve(x, ei).loads[0]
    fine_load = solve(fine, np.repeat(ei, 2)).loads[0]
    assert fine_load <= coarse_load
    restrained = solve(fine, np.repeat(ei, 2), springs=(1, 7)).loads[0]
    assert restrained > fine_load


@pytest.mark.parametrize(
    "nodes,rigidity",
    [
        ([0, 0, 1], 1),
        ([0], 1),
        ([1, 0], 1),
        ([0, np.nan], 1),
        ([[0, 1]], 1),
        ([0, 1], -1),
        ([0, 1], [1, 2]),
        ([0, 1], np.inf),
    ],
)
def test_reject_bad_mesh_and_rigidity(nodes, rigidity):
    with pytest.raises(ValueError):
        solve(nodes, rigidity)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"modes": 0},
        {"modes": 1.2},
        {"modes": True},
        {"modes": 99},
        {"support": "free"},
        {"springs": [-1, 0]},
        {"springs": [0, np.inf]},
        {"springs": [0]},
        {"support": "clamped", "springs": [1, 1]},
    ],
)
def test_reject_bad_options(kwargs):
    with pytest.raises(ValueError):
        solve([0, 0.5, 1], **kwargs)


def test_clamped_one_element_has_no_free_modes():
    with pytest.raises(ValueError):
        solve([0, 1], support="clamped")


@pytest.mark.parametrize("kappa", [-1, np.nan, np.inf])
def test_invalid_reference(kappa):
    with pytest.raises(ValueError):
        spring_reference(kappa)
