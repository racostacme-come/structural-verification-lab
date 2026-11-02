"""Analytical-oracle, convergence and bounded contact verification gates."""
import json
from pathlib import Path
import numpy as np
import pytest
from thermomechanics.heat import HeatMesh
from thermomechanics.coupled_experiment import AREA, INTERFACE, MATERIALS
from thermomechanics.refinement_experiment import cosine_reference, spectral_temperature, experiment


@pytest.fixture(scope='module')
def evidence():
    return experiment()


@pytest.mark.parametrize('amplitude', [-40., 0., 40.])
def test_cosine_oracle_against_independent_gauss_integrals(amplitude):
    # Gaussian integration is independent of the reference sine antiderivative.
    z, weights = np.polynomial.legendre.leggauss(32)
    def thermal(left, right, material):
        x = left+(right-left)*(z+1)/2
        return material.expansion*amplitude*(right-left)/2*np.dot(weights, np.cos(np.pi*x))
    m0, m1 = MATERIALS
    compliance = INTERFACE/(m0.young*AREA)+(1-INTERFACE)/(m1.young*AREA)
    expected_force = -(thermal(0., INTERFACE, m0)+thermal(INTERFACE, 1., m1))/compliance
    x = np.array([0., .13, INTERFACE, .57, 1.])
    force, displacement = cosine_reference(x, amplitude)
    assert force == pytest.approx(expected_force, abs=1e-11)
    expected_u = []
    for position in x:
        split = min(position, INTERFACE)
        expected_u.append(force*(split/(m0.young*AREA)+(position-split)/(m1.young*AREA))
                          +thermal(0., split, m0)+thermal(split, position, m1))
    np.testing.assert_allclose(displacement, expected_u, atol=1e-18)
    np.testing.assert_allclose(displacement[[0, -1]], 0., atol=1e-18)


def test_spectral_reference_unequal_mass_two_cell_closed_form():
    mesh = HeatMesh((0., .25, 1.), (2., 5.), (3., 7.), (.2,), .4)
    initial = np.array([10., -2.])
    mass = mesh.mass
    mean = np.dot(mass, initial)/sum(mass)
    rate = mesh.conductance[1]*sum(1/mass)
    for time in (0., .1, 10.):
        expected = mean+(initial-mean)*np.exp(-rate*time)
        actual = spectral_temperature(mesh, initial, time)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-13)
        assert np.dot(mass, actual) == pytest.approx(np.dot(mass, initial), abs=1e-12)


@pytest.mark.parametrize('quantity', ['temperature_error_K', 'force_error_N', 'nodal_error_m'])
def test_backward_euler_coupled_order(evidence, quantity):
    rows = evidence['temporal_against_semidiscrete']
    assert [row['steps'] for row in rows] == [10, 20, 40, 80]
    assert all(.95 < row['orders'][quantity] < 1.02 for row in rows[1:])


@pytest.mark.parametrize('quantity', ['temperature_error_K', 'force_error_N', 'nodal_error_m'])
def test_smooth_continuum_coupled_order(evidence, quantity):
    rows = evidence['spatial_exact_semidiscrete_against_continuum']
    assert [row['cells'] for row in rows] == [10, 20, 40, 80]
    assert all(1.9 < row['orders'][quantity] < 2.1 for row in rows[1:])
    # Operator spectral evolution independently checked against its known cosine mode.
    assert max(row['spectral_mode_residual_K'] for row in rows) < 1e-9


@pytest.mark.parametrize('resistance', [0., 1e-10, .02, 2., 1e5])
def test_contact_property_limits(evidence, resistance):
    rows = [row for row in evidence['contact_property_sweep'] if row['resistance_m2_K_per_W'] == resistance]
    assert len(rows) == 9
    for row in rows:
        # High-R boundary fluxes subtract near-equal temperatures: absolute
        # temperature/force remain accurate while relative tiny-flux error grows.
        assert row['relative_rate_error'] < (2e-7 if resistance == 1e5 else 1e-10)
        assert row['scaled_temperature_error'] < 1e-11
        assert row['relative_force_error'] < 1e-11
        assert row['heat_condition_number'] < 2e4
        assert 0. <= row['reference_jump_K'] < 80.
    if resistance == 1e5:
        assert min(row['reference_jump_K'] for row in rows) > 79.99


def test_reproduce_all_versioned_values_with_dimensioned_tolerances(evidence):
    saved = json.loads((Path(__file__).resolve().parents[1]/'results/coupled-refinement.json').read_text())
    for key in ('source_sha256', 'parameters', 'scope', 'schema_version'):
        assert evidence[key] == saved[key]
    # Absolute floors are expressed in each output's units. Tiny solver residuals
    # are not required to have a portable relative mantissa; separate gates above
    # enforce their numerical meaning.
    floors = {'nodal_error_m': 1e-14, 'spectral_mode_residual_K': 1e-9,
              'relative_rate_error': 2e-8, 'scaled_temperature_error': 1e-12,
              'relative_force_error': 1e-12}
    for group in ('temporal_against_semidiscrete', 'spatial_exact_semidiscrete_against_continuum', 'contact_property_sweep'):
        assert len(evidence[group]) == len(saved[group])
        for actual, expected in zip(evidence[group], saved[group], strict=True):
            assert actual.keys() == expected.keys()
            for key in actual:
                if key == 'orders':
                    assert actual[key] == pytest.approx(expected[key], rel=1e-6, abs=1e-6)
                else:
                    assert actual[key] == pytest.approx(expected[key], rel=1e-7, abs=floors.get(key, 1e-9))
