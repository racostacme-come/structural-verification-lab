"""Dense Python P1 baseline; intended for small verification problems."""
from dataclasses import dataclass
import numpy as np
from .model import BarMesh, element_temperatures


def assemble(mesh: BarMesh, delta_temperature):
    """Return K (N/m), thermal load f (N); K u = f at free nodes."""
    dt = element_temperatures(mesh, delta_temperature)
    n = len(mesh.nodes)
    stiffness = np.zeros((n, n))
    thermal_load = np.zeros(n)
    for e, (material, area) in enumerate(zip(mesh.materials, mesh.areas)):
        rigidity = material.young * area
        length = mesh.nodes[e+1] - mesh.nodes[e]
        stiffness[e:e+2, e:e+2] += rigidity/length * np.array([[1., -1.], [-1., 1.]])
        thermal_load[e:e+2] += rigidity * material.expansion * dt[e] * np.array([-1., 1.])
    return stiffness, thermal_load


@dataclass
class MechanicalSolution:
    displacement: np.ndarray
    residual: np.ndarray  # K u - f: support forces at constrained nodes
    element_force: np.ndarray  # tension positive


def solve(mesh: BarMesh, delta_temperature, *, left=0.0, right=0.0, assembler=None):
    """Fix end displacement(s); optionally use NativeAssembler.assemble.

    Only assembly is replaced. Boundary elimination and solution use NumPy.
    """
    if not np.isfinite(left) or (right is not None and not np.isfinite(right)):
        raise ValueError("prescribed displacements must be finite")
    dt = element_temperatures(mesh, delta_temperature)
    stiffness, load = (assemble if assembler is None else assembler)(mesh, dt)
    n = len(mesh.nodes)
    constrained = np.array([0] if right is None else [0, n-1])
    free = np.setdiff1d(np.arange(n), constrained)
    u = np.zeros(n)
    u[0] = left
    if right is not None:
        u[-1] = right
    if free.size:
        u[free] = np.linalg.solve(stiffness[np.ix_(free, free)],
                                  load[free] - stiffness[np.ix_(free, constrained)] @ u[constrained])
    force = np.array([m.young*a*((u[e+1]-u[e])/(mesh.nodes[e+1]-mesh.nodes[e]) - m.expansion*dt[e])
                      for e, (m,a) in enumerate(zip(mesh.materials, mesh.areas))])
    return MechanicalSolution(u, stiffness @ u - load, force)
