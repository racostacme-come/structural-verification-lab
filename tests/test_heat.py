"""Independent resistance, balance and analytical eigenmode checks."""
import numpy as np
import pytest
from thermomechanics.heat import HeatMesh, fluxes, operator, step, steady


def mesh(contact=0.):
    return HeatMesh((0., .2, .55, 1.), (2., 8., 8.), (3., 4., 4.), (contact, 0.), .7)


@pytest.mark.parametrize('resistance', [0., 1e-10, .3, 1e5])
def test_series_resistance_and_interface_jump(resistance):
    m = mesh(resistance)
    t = steady(m, left=12., right=2.)
    q = .7*10/(.2/2+.8/8+resistance)
    np.testing.assert_allclose(fluxes(m, t, left=12., right=2.), q, rtol=2e-9)
    # Extrapolate each adjacent cell centre to the contact face.
    hot_face = t[0]-q*.1/(.7*2)
    cold_face = t[1]+q*.175/(.7*8)
    assert hot_face-cold_face == pytest.approx(q*resistance/.7, abs=2e-13)


@pytest.mark.parametrize('dt', [1e-6, .1, 100.])
@pytest.mark.parametrize('boundaries', [(None, None), (10., None), (10., 2.)])
def test_local_and_global_energy_balance(dt, boundaries):
    m = mesh(.4)
    initial = np.array([3., 9., 5.])
    source = np.array([2., -1., 4.])
    left, right = boundaries
    t = step(m, initial, dt, left=left, right=right, source=source)
    q = fluxes(m, t, left=left, right=right)
    power = source*.7*np.diff(m.faces)
    np.testing.assert_allclose(m.mass*(t-initial), dt*(q[:-1]-q[1:]+power), rtol=2e-9, atol=1e-12)
    assert np.dot(m.mass, t-initial) == pytest.approx(dt*(q[0]-q[-1]+sum(power)), rel=2e-9, abs=1e-12)


@pytest.mark.parametrize('dt', [1e-6, .1, 100.])
def test_maximum_principle_and_insulated_conservation(dt):
    m = mesh(.4)
    initial = np.array([1., 9., 3.])
    t = step(m, initial, dt)
    assert 1.-1e-12 <= min(t) <= max(t) <= 9.+1e-12
    assert np.dot(m.mass, t) == pytest.approx(np.dot(m.mass, initial), rel=2e-12)
    np.testing.assert_allclose(step(m, np.ones(3)*7, dt), 7, atol=1e-12)


def test_single_cell_lumped_solution():
    m = HeatMesh((0., 2.), (3.,), (4.,), (), .5)
    # M=4 J/K; each boundary conductance=1.5 W/K; source total=2 W.
    assert step(m, [5.], .2, left=10., right=0., source=[2.])[0] == pytest.approx((20+ .2*17)/4.6)


def test_cosine_discrete_mode_and_first_order_time_refinement():
    n = 12
    m = HeatMesh(tuple(np.linspace(0, 1, n+1)), (1.,)*n, (1.,)*n, (0.,)*(n-1))
    mode = np.cos(np.pi*(np.arange(n)+.5)/n)
    lam = 4*n*n*np.sin(np.pi/(2*n))**2
    h, _ = operator(m)
    np.testing.assert_allclose(h@mode, m.mass*lam*mode, atol=2e-14)
    errors = []
    for count in [10, 20, 40, 80]:
        dt = .1/count
        t = mode.copy()
        for _ in range(count):
            t = step(m, t, dt)
        np.testing.assert_allclose(t, mode*(1+dt*lam)**(-count), atol=2e-14)
        errors.append(max(abs(t-mode*np.exp(-lam*.1))))
    assert all(.9 < np.log2(a/b) < 1.1 for a,b in zip(errors,errors[1:]))


@pytest.mark.parametrize('kwargs', [dict(conductivity=(0.,)), dict(capacity=(-1.,)), dict(contact=(1.,)), dict(area=0.), dict(faces=(0., 0.))])
def test_reject_invalid_mesh(kwargs):
    args = dict(faces=(0., 1.), conductivity=(1.,), capacity=(1.,), contact=())
    args.update(kwargs)
    with pytest.raises(ValueError):
        HeatMesh(**args)


@pytest.mark.parametrize('kwargs', [dict(dt=0.), dict(dt=np.inf), dict(left=np.nan), dict(source=[1.])])
def test_reject_invalid_step(kwargs):
    args = dict(dt=.1)
    args.update(kwargs)
    with pytest.raises(ValueError):
        step(mesh(), [1., 2., 3.], **args)


def test_reject_underdetermined_steady_state():
    with pytest.raises(ValueError):
        steady(mesh())


def test_spatial_spectrum_is_the_assembled_operator():
    # Prevent the spatial study becoming an unconnected analytical calculation.
    for n in [8, 16, 32, 64]:
        m = HeatMesh(tuple(np.linspace(0, 1, n+1)), (1.,)*n, (1.,)*n, (0.,)*(n-1))
        mode = np.cos(np.pi*(np.arange(n)+.5)/n)*np.sinc(1/(2*n))
        h, _ = operator(m)
        lam = 4*n*n*np.sin(np.pi/(2*n))**2
        np.testing.assert_allclose(h@mode, m.mass*lam*mode, atol=1e-13)


def test_versioned_heat_experiment():
    import json
    from pathlib import Path
    from thermomechanics.heat_experiment import experiment
    actual = experiment()
    expected = json.loads((Path(__file__).resolve().parents[1]/'results/thermal-contact.json').read_text())
    assert actual['source_sha256'] == expected['source_sha256']
    for row in actual['contact']:
        assert row['max_relative_rate_error'] < 2e-9
        assert abs(row['transient_energy_residual_J']) < 1e-12
    for key, order in [('temporal_against_semidiscrete', 1.), ('spatial_exact_semidiscrete_against_continuum_cell_average', 2.)]:
        for a,b in zip(actual[key], expected[key], strict=True):
            assert a['error_K'] == pytest.approx(b['error_K'], rel=1e-9)
        assert all(abs(row['observed_order']-order) < .1 for row in actual[key][1:])
