"""Required native integration checks; missing builds fail instead of skipping."""
import json
from pathlib import Path
import numpy as np
import pytest
from thermomechanics import BarMesh, Material, assemble, solve
from thermomechanics.native import NativeAssembler
from thermomechanics.native_experiment import experiment


@pytest.fixture(scope='module')
def native():
    return NativeAssembler()


@pytest.mark.parametrize('dt', [-30., 0., 60.])
@pytest.mark.parametrize('right', [None, 0.])
def test_native_uniform_patch(native, dt, right):
    m = Material(80e9, 12e-6)
    x = np.array([0., .13, .55, 1.7])
    mesh = BarMesh(x, (m,)*3, (3e-4,)*3)
    result = solve(mesh, [dt]*3, right=right, assembler=native.assemble)
    force = 0. if right is None else -m.young*3e-4*m.expansion*dt
    np.testing.assert_allclose(result.displacement, m.expansion*dt*x if right is None else x*0, atol=1e-16)
    np.testing.assert_allclose(result.element_force, force, rtol=1e-12, atol=1e-8)
    np.testing.assert_allclose(result.residual, [-force, 0, 0, force], atol=1e-8)


def test_native_energy_gradient_and_noncontiguous_input(native):
    mesh = BarMesh((0., .2, 1.), (Material(5., .01), Material(20., -.02)), (2., 3.))
    dt = np.array([15., -999., 30., -999.])[::2]
    k, f = native.assemble(mesh, dt)
    np.testing.assert_allclose(k, k.T)
    np.testing.assert_allclose(k@np.ones(3), 0., atol=1e-13)
    u = np.array([.4, -.7, .2])
    h = np.diff(mesh.nodes)
    ea = np.array([10., 60.])
    thermal_strain = np.array([.01, -.02])*dt
    # Differentiation of stored strain energy independently checks load signs.
    def energy(v):
        return .5*np.sum(ea*h*(np.diff(v)/h-thermal_strain)**2)
    eps = 1e-5
    gradient = [(energy(u+eps*v)-energy(u-eps*v))/(2*eps) for v in np.eye(3)]
    np.testing.assert_allclose(k@u-f, gradient, rtol=1e-9, atol=1e-8)
    np.testing.assert_allclose(f.sum(), 0., atol=1e-13)
    np.testing.assert_array_equal(dt, [15., 30.])


@pytest.mark.parametrize('dt', [[], [1., 2.], [float('nan')]])
def test_native_rejects_invalid_temperature(native, dt):
    with pytest.raises(ValueError):
        native.assemble(BarMesh((0., 1.), (Material(1., 0.),), (1.,)), dt)


def test_native_arithmetic_overflow_is_reported(native):
    mesh = BarMesh((0., 1.), (Material(1e308, 1.),), (2.,))
    with pytest.raises(FloatingPointError):
        native.assemble(mesh, [1.])


def test_explicit_missing_library_never_falls_back(tmp_path):
    with pytest.raises(FileNotFoundError):
        NativeAssembler(tmp_path/'missing-library')


def test_native_seeded_experiment_and_versioned_sources():
    actual = experiment()
    expected = json.loads((Path(__file__).resolve().parents[1]/'results/native-thermoelastic.json').read_text())
    assert actual['source_sha256'] == expected['source_sha256']
    assert len(actual['cases']) == len(expected['cases']) == 36
    assert all(error < 1e-11 for error in actual['maximum_scaled_errors'].values())
    # Roundoff metrics need bounded agreement, not byte equality across compilers.
    for a, b in zip(actual['cases'], expected['cases'], strict=True):
        for key in ('elements', 'sample', 'right_displacement_m'):
            assert a[key] == b[key]
        for key in actual['maximum_scaled_errors']:
            assert a[key] == pytest.approx(b[key], abs=1e-11)
