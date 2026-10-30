"""One-way backward-Euler heat -> conservative transfer -> quasi-static bar."""
from dataclasses import dataclass
import numpy as np
from .heat import step
from .mechanics import MechanicalSolution, solve
from .transfer import overlap_temperatures, validate_meshes


@dataclass
class CoupledState:
    temperature: np.ndarray
    mean_temperature_increment: np.ndarray
    mechanics: MechanicalSolution


def advance(heat_mesh, bar_mesh, temperature, dt, *, reference_temperature=0.,
            thermal_left=None, thermal_right=None, source=None,
            displacement_left=0., displacement_right=0., assembler=None):
    """Advance heat one step and solve mechanics without mechanical feedback.

    Temperatures and reference must use the same scale. Boundary heat values
    and sources are evaluated at the end of this step by the caller. Both
    models use the same constant cross-section; variable areas are unsupported.
    An explicit native assembler can replace mechanics assembly, never the solve.
    """
    validate_meshes(heat_mesh, bar_mesh)
    if any(a != heat_mesh.area for a in bar_mesh.areas):
        raise ValueError('coupled heat and mechanics require the same constant area')
    new_temperature = step(heat_mesh, temperature, dt, left=thermal_left,
                           right=thermal_right, source=source)
    means = overlap_temperatures(heat_mesh, new_temperature, bar_mesh,
                                 reference_temperature=reference_temperature)
    mechanics = solve(bar_mesh, means, left=displacement_left,
                      right=displacement_right, assembler=assembler)
    return CoupledState(new_temperature, means, mechanics)
