import numpy as np
import pytest

from shearscope import Beam, _core, assemble, exact_tip, modes, solve_static


@pytest.mark.parametrize("integration", ["full", "reduced"])
def test_element_against_independent_polynomial_integration(integration):
    length, ei, shear = 0.7, 13.0, 41.0
    k, m = _core.element(length, ei, shear, 2.0, 0.4, integration == "reduced")
    # Exact polynomial antiderivative: full shear matrix minus the midpoint term.
    s = np.array(
        [
            [1 / length, 0.5, -1 / length, 0.5],
            [0.5, length / 3, -0.5, length / 6],
            [-1 / length, -0.5, 1 / length, -0.5],
            [0.5, length / 6, -0.5, length / 3],
        ]
    )
    if integration == "reduced":
        s[np.ix_([1, 3], [1, 3])] = length / 4
    bending = np.zeros((4, 4))
    bending[np.ix_([1, 3], [1, 3])] = ei / length * np.array([[1, -1], [-1, 1]])
    np.testing.assert_allclose(k, bending + shear * s, atol=1e-13)
    assert np.linalg.eigvalsh(m).min() > 0
    np.testing.assert_allclose(k @ [0, 1, length, 1], 0, atol=1e-13)
    assert np.count_nonzero(np.linalg.eigvalsh(k) > 1e-9) == 2


@pytest.mark.parametrize(
    "loads", [{"force": 2.0}, {"force": 0, "moment": 0.3}, {"force": 0, "distributed": 3.0}]
)
def test_analytical_tip_and_global_equilibrium(loads):
    beam = Beam(height=0.05)
    result = solve_static(beam, elements=64, **loads)
    assert result.displacement[-1] == pytest.approx(exact_tip(beam, **loads), rel=2e-4)
    force, moment, q = (
        loads.get("force", 1),
        loads.get("moment", 0),
        loads.get("distributed", 0),
    )
    np.testing.assert_allclose(
        result.reactions,
        [-force - q * beam.length, -moment - force * beam.length - q * beam.length**2 / 2],
        atol=1e-7,
    )
    assert result.relative_residual < 1e-7


def test_quadratic_convergence_and_locking():
    beam = Beam(height=0.01)
    errors = [
        abs(solve_static(beam, n).displacement[-1] / exact_tip(beam) - 1)
        for n in (4, 8, 16, 32)
    ]
    np.testing.assert_allclose(np.array(errors[:-1]) / errors[1:], 4, rtol=0.002)
    full = solve_static(beam, 8, "full").displacement[-1] / exact_tip(beam)
    reduced = solve_static(beam, 8).displacement[-1] / exact_tip(beam)
    assert full < 0.03
    assert reduced > 0.996


def test_nonuniform_mesh_and_energy():
    beam = Beam()
    x = beam.length * np.linspace(0, 1, 65) ** 1.3
    result = solve_static(beam, nodes=x, force=2.0)
    assert result.displacement[-1] == pytest.approx(exact_tip(beam, 2), rel=2e-4)
    assert result.strain_energy == pytest.approx(result.displacement[-1], rel=1e-8)


def test_modes_euler_bernoulli_limit_and_orthogonality():
    beam = Beam(height=0.005)
    frequencies, shapes = modes(beam, 80, count=3)
    roots = np.array([1.875104068711961, 4.694091132974174, 7.854757438237613])
    reference = (
        roots**2 / (2 * np.pi * beam.length**2) * np.sqrt(beam.ei / (beam.density * beam.area))
    )
    np.testing.assert_allclose(frequencies, reference, rtol=0.002)
    _, k, m = assemble(beam, 80)
    np.testing.assert_allclose(shapes.T @ m @ shapes, np.eye(3), atol=1e-10)
    residual = (k @ shapes - (m @ shapes) * (2 * np.pi * frequencies) ** 2)[2:]
    assert np.linalg.norm(residual) / np.linalg.norm((k @ shapes)[2:]) < 1e-5


@pytest.mark.parametrize(
    "kwargs",
    [{"height": 0}, {"length": -1}, {"young": np.nan}, {"poisson": 0.5}, {"density": np.inf}],
)
def test_invalid_beam(kwargs):
    with pytest.raises(ValueError):
        Beam(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"elements": 0},
        {"elements": True},
        {"elements": 1.5},
        {"integration": "typo"},
        {"nodes": [0, 1, 0.5]},
        {"nodes": [0, np.nan, 1]},
        {"nodes": [[0, 1]]},
    ],
)
def test_invalid_mesh(kwargs):
    with pytest.raises(ValueError):
        assemble(Beam(), **kwargs)


def test_zero_negative_and_invalid_loads():
    beam = Beam()
    assert solve_static(beam, force=0).strain_energy == 0
    assert solve_static(beam, force=-1).displacement[-1] < 0
    with pytest.raises(ValueError):
        solve_static(beam, force=np.inf)
    with pytest.raises(ValueError):
        modes(beam, count=0)
    with pytest.raises(ValueError):
        _core.element(0, 1, 1, 1, 1)
