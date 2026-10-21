import numpy as np
import pytest

from bucklemap import solve, spring_reference


@pytest.mark.parametrize("kappa", [0.1, 1.0, 10.0, 100.0])
def test_modal_derivative_against_perturbed_eigenvalue(kappa):
    x = np.linspace(0, 1, 17)
    q = solve(x, springs=(kappa, kappa)).modes[:, 0]
    step = kappa * 1e-3
    plus = solve(x, springs=(kappa + step, kappa + step)).loads[0]
    minus = solve(x, springs=(kappa - step, kappa - step)).loads[0]
    finite_difference = (plus - minus) / (2 * step)
    np.testing.assert_allclose(q[1] ** 2 + q[-1] ** 2, finite_difference, rtol=3e-6)


def test_restraint_limits():
    assert spring_reference(0) == np.pi**2
    np.testing.assert_allclose(spring_reference(1e12), 4 * np.pi**2, rtol=1e-10)
    values = [spring_reference(k) for k in (0, 0.1, 1, 10, 100)]
    assert np.all(np.diff(values) > 0)
