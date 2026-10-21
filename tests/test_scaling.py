"""Independent dimensional and continuum checks across the two beam APIs."""
from dataclasses import replace

import numpy as np
import pytest

from bucklemap.model import solve
from shearscope.model import Beam, solve_static


def outputs(beam, elements):
    load = solve(np.linspace(0, beam.length, elements + 1), beam.ei,
                 support="cantilever").loads[0]
    compliance = solve_static(beam, elements, "reduced").displacement[-1]
    return load, compliance


@pytest.mark.parametrize("factor", [0.1, 3.0, 10.0])
def test_young_modulus_scaling(factor):
    beam = Beam()
    load, compliance = outputs(beam, 16)
    scaled_load, scaled_compliance = outputs(replace(beam, young=beam.young*factor), 16)
    assert scaled_load / load == pytest.approx(factor, rel=2e-7)
    assert scaled_compliance / compliance == pytest.approx(1/factor, rel=2e-8)


@pytest.mark.parametrize("scale", [0.25, 4.0])
def test_geometrically_similar_beams(scale):
    beam = Beam()
    load, compliance = outputs(beam, 16)
    scaled = replace(beam, length=scale*beam.length, width=scale*beam.width,
                     height=scale*beam.height)
    scaled_load, scaled_compliance = outputs(scaled, 16)
    # EI scales as length^4, Euler load as length^2, compliance as 1/length.
    assert scaled_load / load == pytest.approx(scale**2, rel=2e-7)
    assert scaled_compliance / compliance == pytest.approx(1/scale, rel=2e-8)


@pytest.mark.parametrize("slenderness", [20.0, 100.0])
def test_continuum_limits(slenderness):
    beam = Beam(height=1/slenderness)
    exact_load = np.pi**2 * beam.ei / (4*beam.length**2)
    exact_compliance = beam.length**3/(3*beam.ei) + beam.length/beam.shear
    errors = []
    for n in (8, 16, 32):
        load, compliance = outputs(beam, n)
        errors.append((abs(load/exact_load-1), abs(compliance/exact_compliance-1)))
    errors = np.array(errors)
    # Hermite buckling converges at order four; reduced static bending at order two.
    assert np.all(errors[:-1, 0] / errors[1:, 0] > 10)
    assert np.all(errors[:-1, 1] / errors[1:, 1] > 3.8)
    assert errors[-1, 0] < 1e-6
    assert errors[-1, 1] < 3e-4
