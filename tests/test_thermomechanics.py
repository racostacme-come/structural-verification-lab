"""Analytical patches, signed reactions and independent continuum checks."""
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from thermomechanics import BarMesh, Material, Layer, assemble, solve, continuum_solution
from thermomechanics.__main__ import experiment


@pytest.mark.parametrize("dt", [-35., 0., 60.])
@pytest.mark.parametrize("right", [None, 0.])
def test_uniform_patch(dt, right):
    m = Material(80e9, 12e-6)
    x = np.array([0., .13, .55, 1.7])
    mesh = BarMesh(x, (m,)*3, (3e-4,)*3)
    result = solve(mesh, [dt]*3, right=right)
    force = 0. if right is None else -m.young*3e-4*m.expansion*dt
    expected_u = m.expansion*dt*x if right is None else x*0
    np.testing.assert_allclose(result.displacement, expected_u, atol=1e-16, rtol=1e-12)
    np.testing.assert_allclose(result.element_force, force, atol=1e-8, rtol=1e-12)
    np.testing.assert_allclose(result.residual, [-force,0,0,force], atol=1e-8)


def test_two_material_hand_calculation_and_nonzero_boundary():
    # Deliberately dimensionless-looking SI values permit hand arithmetic:
    # compliance = 2/(10*2)+3/(20*1)=.25 m/N; thermal extension = .8 m.
    a, b = Material(10., .1), Material(20., .2)
    mesh = BarMesh((0.,2.,5.), (a,b), (2.,1.))
    result = solve(mesh, (1.,1.), left=.3, right=.6)
    np.testing.assert_allclose(result.element_force, -2., atol=1e-14)
    np.testing.assert_allclose(result.displacement, [.3,.3,.6], atol=1e-14)
    np.testing.assert_allclose(result.residual, [2.,0.,-2.], atol=1e-14)
    ref = continuum_solution([Layer(2.,2.,a,1.,1.),Layer(3.,1.,b,1.,1.)],left=.3,right=.6)
    assert ref.force == pytest.approx(-2.)
    np.testing.assert_allclose(ref.displacement([0,2,5]), result.displacement, atol=1e-14)


def test_affine_reference_is_quadratic_and_reactions_have_sign():
    m = Material(100., .01)
    ref = continuum_solution([Layer(2., 3., m, 0., 10.)])
    assert ref.force == pytest.approx(-15.)
    assert ref.reactions == pytest.approx((15., -15.))
    np.testing.assert_allclose(ref.displacement([0.,1.,2.]), [0.,-.025,0.], atol=1e-16)


def test_stiffness_nullspace_symmetry_and_positive_energy():
    mesh = BarMesh((0.,.2,1.), (Material(5.,.01),Material(20.,-.02)), (2.,3.))
    k,f = assemble(mesh, (15.,30.))
    np.testing.assert_allclose(k,k.T)
    np.testing.assert_allclose(k@np.ones(3), 0., atol=1e-14)
    assert f.sum() == pytest.approx(0.)
    u = np.array([.4,-.7,.2])
    assert u@k@u == pytest.approx(5*2/.2*1.1**2+20*3/.8*.9**2)


def test_random_layered_affine_fields_match_independent_integrals():
    rng = np.random.default_rng(20261027)
    for _ in range(12):
        lengths = rng.uniform(.1,.7,5)
        materials = tuple(Material(e,a) for e,a in zip(rng.uniform(20e9,200e9,5),rng.uniform(-5e-6,25e-6,5)))
        areas = rng.uniform(1e-4,5e-4,5)
        tl,tr = rng.uniform(-30,100,(2,5))
        layers = [Layer(h,a,m,l,r) for h,a,m,l,r in zip(lengths,areas,materials,tl,tr)]
        x = np.r_[0.,np.cumsum(lengths)]
        # Avoid endpoint roundoff differing from math.fsum's oracle length.
        x[-1] = math.fsum(layer.length for layer in layers)
        mesh = BarMesh(x,materials,areas)
        left,right = rng.uniform(-1e-4,1e-4,2)
        ref = continuum_solution(layers,left=left,right=right)
        result = solve(mesh,(tl+tr)/2,left=left,right=right)
        np.testing.assert_allclose(result.element_force,ref.force,rtol=1e-11,atol=1e-8)
        np.testing.assert_allclose(result.displacement,ref.displacement(x),rtol=1e-11,atol=1e-15)


def test_one_element_with_both_nodes_constrained():
    mesh = BarMesh((0.,1.),(Material(10.,.01),),(2.,))
    result = solve(mesh,(10.,),left=.2,right=.4)
    np.testing.assert_allclose(result.element_force,[2.])
    np.testing.assert_allclose(result.residual,[-2.,2.])


@pytest.mark.parametrize("nodes,areas", [((0.,0.),(1.,)),((0.,float('nan')),(1.,)),((0.,1.),(-1.,)),((0.,1.),())])
def test_bad_mesh_rejected(nodes,areas):
    with pytest.raises(ValueError):
        BarMesh(nodes,(Material(1.,0.),),areas)


@pytest.mark.parametrize("young,alpha", [(0.,0.),(-1.,0.),(float('inf'),0.),(1.,float('nan'))])
def test_bad_material_rejected(young,alpha):
    with pytest.raises(ValueError):
        Material(young,alpha)


@pytest.mark.parametrize("dt", [[],[1.,2.],[float('nan')]])
def test_bad_temperature_rejected(dt):
    with pytest.raises(ValueError):
        solve(BarMesh((0.,1.),(Material(1.,0.),),(1.,)),dt)


def test_experiment_separates_nodal_exactness_from_continuum_error():
    data = experiment()
    for case in data['cases']:
        assert case['scaled_force_error'] < 1e-10
        assert case['scaled_nodal_displacement_error'] < 1e-10
    levels = data['affine_temperature_refinement']
    for row in levels:
        # Exact midpoint interpolation error: alpha * temperature slope * h^2 / 8.
        assert row['midpoint_displacement_error_m'] == pytest.approx(23e-6*80/(8*row['elements']**2),rel=1e-9)
        assert row['max_nodal_error_m'] < 1e-15
    assert all(abs(row['observed_order']-2.) < 1e-8 for row in levels[1:])


def test_versioned_experiment_reproduces():
    expected = json.loads((Path(__file__).resolve().parents[1]/'results/thermoelastic-baseline.json').read_text())
    actual = experiment()
    assert actual['source_sha256'] == expected['source_sha256']
    for a,b in zip(actual['cases'],expected['cases'],strict=True):
        assert a['name'] == b['name']
        for key in ('reference_force_N','element_force_N','displacement_m','support_force_N'):
            np.testing.assert_allclose(a[key],b[key],rtol=1e-11,atol=1e-12)
    for a,b in zip(actual['affine_temperature_refinement'],expected['affine_temperature_refinement'],strict=True):
        assert a['elements'] == b['elements']
        assert a['midpoint_displacement_error_m'] == pytest.approx(b['midpoint_displacement_error_m'],rel=1e-9)


def test_cli_writes_json(tmp_path):
    path = tmp_path/'baseline.json'
    subprocess.run([sys.executable,'-m','thermomechanics','--output',str(path)],check=True)
    assert len(json.loads(path.read_text())['cases']) == 4
