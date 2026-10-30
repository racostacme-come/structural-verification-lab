"""Independent integral checks and failure contracts for nonmatching transfer."""
import json
from pathlib import Path
import numpy as np
import pytest
from thermomechanics.model import BarMesh, Material
from thermomechanics.heat import HeatMesh
from thermomechanics.transfer import overlap_temperatures
from thermomechanics.coupled import advance
from thermomechanics.native import NativeAssembler
from thermomechanics.coupled_experiment import experiment


def meshes():
    heat = HeatMesh((0., .2, .4, .7, 1.), (1.,)*4, (2.,)*4, (0., .1, 0.))
    bar = BarMesh((0., .1, .4, .6, 1.),
                  (Material(20., 2e-5),)*2+(Material(50., -1e-5),)*2, (1.,)*4)
    return heat, bar


def test_hand_integrated_discontinuous_field_and_signed_expansion():
    heat, bar = meshes()
    temperature = np.array([310., 330., 280., 340.])
    result = overlap_temperatures(heat, temperature, bar, reference_temperature=300.)
    np.testing.assert_allclose(result, [10., 70./3., -20., 25.], atol=1e-13)
    alpha = np.array([m.expansion for m in bar.materials])
    assert np.dot(alpha*np.diff(bar.nodes), result) == pytest.approx(1e-4, abs=1e-18)
    np.testing.assert_array_equal(temperature, [310., 330., 280., 340.])


@pytest.mark.parametrize('seed', range(5))
def test_random_partition_conservation_and_convexity(seed):
    rng = np.random.default_rng(seed)
    faces = np.r_[0., np.sort(rng.uniform(.01, .99, 19)), 1.]
    nodes = np.r_[0., np.sort(rng.uniform(.01, .99, 11)), 1.]
    heat = HeatMesh(faces, (1.,)*20, (1.,)*20, (0.,)*19)
    bar = BarMesh(nodes, (Material(1., 0.),)*12, (1.,)*12)
    field = rng.normal(size=20)
    mapped = overlap_temperatures(heat, field, bar)
    assert np.dot(np.diff(nodes), mapped) == pytest.approx(np.dot(np.diff(faces), field), abs=2e-15)
    roundoff = 8*np.finfo(float).eps*max(1., max(abs(field)))
    assert min(mapped) >= min(field)-roundoff and max(mapped) <= max(field)+roundoff
    np.testing.assert_allclose(overlap_temperatures(heat, np.full(20, 7.), bar), 7., atol=1e-13)


def test_identity_map_preserves_contact_jump():
    heat, _ = meshes()
    bar = BarMesh(heat.faces, (Material(1., 1.),)*4, (1.,)*4)
    np.testing.assert_allclose(overlap_temperatures(heat, [1., 3., -20., 40.], bar), [1., 3., -20., 40.])


@pytest.mark.parametrize('field', [[], [1., 2.], [1., 2., 3., np.nan], [1., 2., 3., np.inf]])
def test_invalid_field(field):
    heat, bar = meshes()
    with pytest.raises(ValueError, match='one finite value'):
        overlap_temperatures(heat, field, bar)


def test_invalid_reference():
    heat, bar = meshes()
    with pytest.raises(ValueError, match='reference temperature'):
        overlap_temperatures(heat, [1., 2., 3., 4.], bar, reference_temperature=np.nan)


@pytest.mark.parametrize('nodes', [(0., .9), (-.1, 1.), (0., 1.1)])
def test_domain_mismatch_never_extrapolates(nodes):
    heat, _ = meshes()
    bar = BarMesh(nodes, (Material(1., 1.),), (1.,))
    with pytest.raises(ValueError, match='identical endpoints'):
        overlap_temperatures(heat, [1., 2., 3., 4.], bar)


@pytest.mark.parametrize('change', ['material', 'area'])
def test_unresolved_property_jump(change):
    heat, _ = meshes()
    materials = (Material(1., 1.), Material(2. if change == 'material' else 1., 1.))
    bar = BarMesh((0., .5, 1.), materials, (1., 2. if change == 'area' else 1.))
    with pytest.raises(ValueError, match='property interfaces'):
        overlap_temperatures(heat, [1., 2., 3., 4.], bar)


def test_coupling_rejects_area_mismatch():
    heat, bar = meshes()
    bar = BarMesh(bar.nodes, bar.materials, (2.,)*4)
    with pytest.raises(ValueError, match='same constant area'):
        advance(heat, bar, [1., 2., 3., 4.], 1.)


@pytest.mark.parametrize('native', [False, True])
@pytest.mark.parametrize('right', [None, 0., .02])
def test_insulated_uniform_coupled_patch(native, right):
    heat = HeatMesh((0., .3, 1.), (2., 3.), (4., 5.), (.7,), 2.)
    material = Material(100., .001)
    bar = BarMesh((0., .1, .4, .8, 1.), (material,)*4, (2.,)*4)
    assembler = NativeAssembler().assemble if native else None
    field = np.array([320., 320.])
    state = advance(heat, bar, field, 100., reference_temperature=300.,
                    displacement_right=right, assembler=assembler)
    np.testing.assert_allclose(state.temperature, 320., atol=1e-12)
    np.testing.assert_allclose(state.mean_temperature_increment, 20., atol=1e-12)
    expected_force = 0. if right is None else 200.*(right-.02)
    expected_u = np.array(bar.nodes)*(.02 if right is None else right)
    np.testing.assert_allclose(state.mechanics.displacement, expected_u, atol=1e-14)
    np.testing.assert_allclose(state.mechanics.element_force, expected_force, atol=1e-12)
    np.testing.assert_array_equal(field, [320., 320.])


def test_coupled_experiment_reproduces_versioned_evidence():
    actual = experiment()
    expected = json.loads((Path(__file__).resolve().parents[1]/'results/coupled-transfer.json').read_text())
    assert actual['source_sha256'] == expected['source_sha256']
    assert actual['parameters'] == expected['parameters']
    assert len(actual['transfer_only']) == 4
    assert max(r['scaled_integral_error'] for r in actual['transfer_only']) < 5e-12
    assert max(r['overlap_relative_force_error'] for r in actual['transfer_only']) < 1e-11
    assert all(v < 1e-11 for v in actual['maximum_scaled_errors'].values())
    assert len(actual['coupled_histories']) == len(expected['coupled_histories']) == 3
    for a, b in zip(actual['transfer_only'], expected['transfer_only'], strict=True):
        assert a.keys() == b.keys()
        for key in a:
            assert a[key] == pytest.approx(b[key], rel=1e-10, abs=1e-12)
    for history, saved in zip(actual['coupled_histories'], expected['coupled_histories'], strict=True):
        assert history['contact_resistance_m2_K_per_W'] == saved['contact_resistance_m2_K_per_W']
        assert len(history['steps']) == len(saved['steps']) == 20
        for a, b in zip(history['steps'], saved['steps'], strict=True):
            for key in a:
                np.testing.assert_allclose(a[key], b[key], rtol=1e-10, atol=1e-11)
