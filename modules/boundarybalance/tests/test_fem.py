import numpy as np
import pytest
from scipy.linalg import eigvalsh

from boundarybalance.fem import (
    assemble,
    critical_penalty,
    errors,
    solve,
    triangle_rule,
    unit_square,
)
from boundarybalance.problems import exact, gradient, source


def affine(p):
    return 0.7 + 1.3 * p[:, 0] - 2.1 * p[:, 1]


def zero(p):
    return np.zeros(len(p))


@pytest.mark.parametrize("n", [2, 3, 7])
def test_geometry_and_operators(n):
    mesh = unit_square(n)
    assert len(mesh.triangles) == 2 * n**2
    assert len(mesh.boundary) == 4 * n
    system = assemble(mesh, zero, affine)
    k, c, b = system.stiffness, system.flux, system.boundary_mass
    one = np.ones(len(mesh.points))
    x, y = mesh.points.T
    assert np.max(abs(k @ one)) < 2e-14
    assert one @ b @ one == pytest.approx(4 * n)
    assert x @ k @ x == pytest.approx(1)
    assert y @ k @ y == pytest.approx(1)
    assert x @ k @ y == pytest.approx(0, abs=1e-14)
    # Divergence theorem independently checks outward normals and index convention.
    assert x @ c @ x == pytest.approx(1)
    assert y @ c @ y == pytest.approx(1)
    np.testing.assert_allclose(k.toarray(), k.T.toarray(), atol=1e-14)
    np.testing.assert_allclose(system.matrix().toarray(), system.matrix().T.toarray())


@pytest.mark.parametrize("method", ["strong", "nitsche"])
@pytest.mark.parametrize("n", [2, 5, 9])
def test_affine_patch(method, n):
    system = assemble(unit_square(n), zero, affine)
    result = solve(system, method)
    np.testing.assert_allclose(result.values, affine(system.mesh.points), atol=2e-13)
    assert result.backward_error < 1e-14


def test_penalty_has_consistency_bias_even_on_affine_patch():
    system = assemble(unit_square(6), zero, affine)
    weak = solve(system, "penalty", 10)
    stiff = solve(system, "penalty", 1000)
    target = affine(system.mesh.points)
    assert np.max(abs(weak.values - target)) > 0.01
    assert np.linalg.norm(stiff.values - target) < np.linalg.norm(weak.values - target) / 80


@pytest.mark.parametrize("n", [2, 4, 8])
def test_schur_threshold_matches_full_matrix_inertia(n):
    system = assemble(unit_square(n), source, exact)
    threshold = critical_penalty(system)
    assert 1 < threshold < 10
    assert eigvalsh(system.matrix(gamma=threshold * 0.999).toarray())[0] < 0
    assert abs(eigvalsh(system.matrix(gamma=threshold).toarray())[0]) < 2e-13
    assert eigvalsh(system.matrix(gamma=threshold * 1.001).toarray())[0] > 0
    with pytest.raises(ValueError, match="SPD threshold"):
        solve(system, gamma=threshold * 0.999)


@pytest.mark.parametrize("method", ["strong", "nitsche"])
def test_continuum_convergence(method):
    records = []
    for n in (8, 16, 32):
        mesh = unit_square(n)
        result = solve(assemble(mesh, source, exact), method)
        records.append(errors(mesh, result.values, exact, gradient))
        assert result.backward_error < 2e-14
    order_l2 = np.log2(records[-2]["l2"] / records[-1]["l2"])
    order_h1 = np.log2(records[-2]["h1"] / records[-1]["h1"])
    assert 1.85 < order_l2 < 2.15
    assert 0.9 < order_h1 < 1.1


def test_quadrature_refinement_and_affine_error_norm():
    mesh = unit_square(4)
    systems = [assemble(mesh, source, exact, order=q) for q in (6, 10)]
    np.testing.assert_allclose(systems[0].load, systems[1].load, atol=2e-13)
    err = errors(mesh, affine(mesh.points), affine, lambda p: np.tile([1.3, -2.1], (len(p), 1)))
    assert max(err.values()) < 5e-15
    basis, weights = triangle_rule()
    assert weights.sum() == pytest.approx(0.5)
    np.testing.assert_allclose(weights @ basis, np.full(3, 1 / 6))


@pytest.mark.parametrize("n", [True, 1, -3, 2.5, float("nan")])
def test_invalid_mesh(n):
    with pytest.raises(ValueError):
        unit_square(n)


@pytest.mark.parametrize("gamma", [0, -1, float("nan"), float("inf")])
def test_invalid_penalty(gamma):
    with pytest.raises(ValueError):
        solve(assemble(unit_square(2), source, exact), gamma=gamma)


@pytest.mark.parametrize("callback", [lambda p: 0, lambda p: np.full(len(p), np.nan)])
def test_invalid_callback(callback):
    with pytest.raises(ValueError, match="callback"):
        assemble(unit_square(2), callback, exact)


def test_zero_problem_and_bad_method():
    system = assemble(unit_square(2), zero, zero)
    for method in ("strong", "penalty", "nitsche"):
        solution = solve(system, method)
        assert np.array_equal(solution.values, np.zeros(9))
        assert solution.backward_error == 0
    with pytest.raises(ValueError, match="method"):
        solve(system, "unknown")


def test_constant_source_load_integral():
    system = assemble(unit_square(5), lambda p: np.full(len(p), 3.0), zero)
    assert np.sum(system.load) == pytest.approx(3)
