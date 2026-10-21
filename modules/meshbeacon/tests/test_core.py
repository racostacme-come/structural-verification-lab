import numpy as np
import pytest

from meshbeacon.core import adapt, errors, mark, mesh, refine, solve
from meshbeacon.problems import bump, sine


def test_patch_nonuniform_nonzero_boundary():
    x = [0, 0.07, 0.3, 0.9, 1]
    s = solve(x, lambda x: 0, diffusion=3, boundary=(2, 5))
    np.testing.assert_allclose(s.values, 2 + 3 * np.array(x), atol=1e-14)
    np.testing.assert_allclose(s.flux, 9, atol=1e-13)
    assert s.estimator < 1e-13
    with pytest.raises(ValueError):
        s.values[1] = 2


def test_quadratic_nodal_exactness_and_energy_error():
    x = np.array([0, 0.1, 0.4, 0.7, 1])
    s = solve(x, lambda x: 2, diffusion=1)
    np.testing.assert_allclose(s.values, x * (1 - x), atol=1e-14)
    l2, energy = errors(s, lambda x: x * (1 - x), lambda x: 1 - 2 * x)
    assert energy == pytest.approx(np.sqrt(np.sum(np.diff(x) ** 3) / 3))
    assert l2 == pytest.approx(np.sqrt(np.sum(np.diff(x) ** 5) / 30))


def test_hand_computed_indicator_and_single_element():
    s = solve([0, 0.5, 1], lambda x: 2)
    np.testing.assert_allclose(s.indicators_squared, [0.75, 0.75])
    s = solve([0, 1], lambda x: 2, boundary=(1, 3))
    np.testing.assert_allclose(s.values, [1, 3])
    assert s.estimator == pytest.approx(2)


def test_diffusion_scaling():
    p = sine(7)
    s = solve(np.linspace(0, 1, 17), p.load, diffusion=7)
    ref = solve(s.nodes, sine().load)
    np.testing.assert_allclose(s.values, ref.values)
    assert s.estimator == pytest.approx(np.sqrt(7) * ref.estimator)


def test_sine_convergence():
    p = sine()
    e = [
        errors(solve(np.linspace(0, 1, n + 1), p.load), p.exact, p.gradient)
        for n in (16, 32, 64)
    ]
    orders = np.log2(np.array(e[:-1]) / np.array(e[1:]))
    np.testing.assert_allclose(orders, [[2, 1], [2, 1]], atol=0.003)


def test_bulk_marking_minimal_and_deterministic():
    np.testing.assert_array_equal(mark([1, 4, 4, 1], 0.5), [1, 2])
    np.testing.assert_array_equal(mark([1e300, 1e300, 1], 0.5), [0])
    np.testing.assert_array_equal(mark([1e300, 1, 0], 1), [0, 1])
    assert mark([0, 0]).size == 0
    rng = np.random.default_rng(1729)
    for _ in range(20):
        eta = rng.random(21)
        selected = mark(eta, 0.6)
        assert eta[selected].sum() >= 0.6 * eta.sum()
        assert eta[selected[:-1]].sum() < 0.6 * eta.sum()


def test_refinement_nested_and_duplicates():
    refined = refine([0, 0.2, 1], [1, 1])
    np.testing.assert_allclose(refined, [0, 0.2, 0.6, 1], rtol=0, atol=2e-16)
    np.testing.assert_array_equal(refined[[0, 1, 3]], [0, 0.2, 1])
    np.testing.assert_array_equal(refine([0, 1], []), [0, 1])
    with pytest.raises(FloatingPointError):
        refine([0, np.nextafter(0.0, 1.0), 1], [0])


def test_adaptive_accuracy_and_localization():
    p = bump()
    result = adapt(p.load, tolerance=0.3)
    assert result.stop_reason == "tolerance"
    s = result.solution
    uniform = solve(np.linspace(0, 1, s.elements + 1), p.load)
    energy = errors(s, p.exact, p.gradient)[1]
    uniform_energy = errors(uniform, p.exact, p.gradient)[1]
    assert energy < 0.4 * uniform_energy
    mids = (s.nodes[1:] + s.nodes[:-1]) / 2
    assert np.count_nonzero(abs(mids - 0.37) < 0.12) > s.elements * 0.7
    previous = np.inf
    for solution in result.history:
        current = errors(solution, p.exact, p.gradient)[1]
        assert current <= previous * (1 + 1e-8)
        previous = current


def test_stopping_reasons_and_zero_load():
    assert adapt(lambda x: 0).stop_reason == "tolerance"
    assert adapt(lambda x: 1, max_steps=1).stop_reason == "max_steps"
    r = adapt(lambda x: 1, max_elements=8, tolerance=1e-10)
    assert r.stop_reason == "max_elements"
    assert r.solution.elements == 8


def test_interpolation():
    s = solve([0, 0.5, 1], lambda x: 0, boundary=(1, 3))
    np.testing.assert_allclose(s.evaluate([[0, 0.25], [0.75, 1]]), [[1, 1.5], [2.5, 3]])
    assert s.evaluate(0.5) == pytest.approx(2)
    with pytest.raises(ValueError):
        s.evaluate([1.1])


@pytest.mark.parametrize("nodes", [[0], [[0, 1]], [0, 0.5, 0.5, 1], [0, np.nan, 1], [0.1, 1]])
def test_invalid_mesh(nodes):
    with pytest.raises(ValueError):
        mesh(nodes)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"diffusion": 0},
        {"diffusion": np.inf},
        {"boundary": (1,)},
        {"boundary": (0, np.nan)},
        {"quadrature": 1},
        {"quadrature": 33},
        {"quadrature": 3.5},
        {"quadrature": True},
    ],
)
def test_invalid_solve(kwargs):
    with pytest.raises(ValueError):
        solve([0, 0.5, 1], lambda x: 1, **kwargs)


@pytest.mark.parametrize("field", [lambda x: np.nan, lambda x: np.ones(3)])
def test_bad_field(field):
    with pytest.raises(ValueError):
        solve([0, 0.5, 1], field)


@pytest.mark.parametrize(
    "eta,theta",
    [([-1], 0.5), ([np.inf], 0.5), ([[1]], 0.5), ([1], 0), ([1], np.nan), ([1], 1.1)],
)
def test_bad_mark(eta, theta):
    with pytest.raises(ValueError):
        mark(eta, theta)


@pytest.mark.parametrize("chosen", [[-1], [1], [0.5], [[0]], [True]])
def test_bad_refine(chosen):
    with pytest.raises(ValueError):
        refine([0, 1], chosen)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tolerance": 0},
        {"theta": 0},
        {"max_elements": 4},
        {"max_steps": 0},
        {"max_elements": 8.5},
    ],
)
def test_bad_adapt(kwargs):
    with pytest.raises(ValueError):
        adapt(lambda x: 1, **kwargs)


@pytest.mark.parametrize("kwargs", [{"width": 0}, {"center": 1}, {"center": np.nan}])
def test_bad_bump(kwargs):
    with pytest.raises(ValueError):
        bump(**kwargs)
