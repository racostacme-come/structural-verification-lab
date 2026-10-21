import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from bucklemap import element_matrices, reconstruct


@pytest.mark.parametrize("h,ei", [(0.17, 3.4), (2.3, 0.07), (1, 1)])
def test_native_matrices_against_independent_polynomial_energy(h, ei):
    # Fit physical cubics to nodal data independently of the native formulas.
    interpolation = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [1, h, h**2, h**3], [0, 1, 2 * h, 3 * h**2]]
    )
    polynomials = np.linalg.solve(interpolation, np.eye(4))
    z, weight = leggauss(5)
    x, w = h * (z + 1) / 2, h * weight / 2
    first = np.column_stack([0 * x, 0 * x + 1, 2 * x, 3 * x * x]) @ polynomials
    second = np.column_stack([0 * x, 0 * x, 0 * x + 2, 6 * x]) @ polynomials
    k, g = element_matrices(h, ei)
    np.testing.assert_allclose(k, ei * second.T @ (w[:, None] * second), rtol=2e-14)
    np.testing.assert_allclose(g, first.T @ (w[:, None] * first), rtol=2e-14)


@pytest.mark.parametrize(
    "h,ei", [(0, 1), (-1, 1), (1, 0), (1, -1), (float("nan"), 1), (1, float("inf"))]
)
def test_invalid_native_parameters(h, ei):
    with pytest.raises(ValueError):
        element_matrices(h, ei)


def test_cubic_reconstruction_on_nonuniform_mesh():
    nodes = np.array([-1, -0.3, 0.2, 2.0])
    q = np.empty(2 * len(nodes))
    q[::2] = 2 - nodes + 3 * nodes**2 - 0.7 * nodes**3
    q[1::2] = -1 + 6 * nodes - 2.1 * nodes**2
    points = np.linspace(-1, 2, 71)
    np.testing.assert_allclose(
        reconstruct(nodes, q, points), 2 - points + 3 * points**2 - 0.7 * points**3, atol=3e-14
    )


def test_no_extrapolation():
    with pytest.raises(ValueError):
        reconstruct([0, 1], [0] * 4, [1.01])


def test_native_overflow():
    with pytest.raises(OverflowError):
        element_matrices(1e-200, 1.0)
