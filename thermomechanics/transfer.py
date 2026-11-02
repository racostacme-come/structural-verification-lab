"""Conservative integration of a piecewise-constant thermal field, in SI units.

Both grids cover exactly the same interval. Mechanical property jumps must be
thermal faces; interior grid points otherwise need not match. No interpolation
across contact jumps is performed. Dense overlaps target small verification runs.
"""
import numpy as np
from .heat import _field


def validate_meshes(heat_mesh, bar_mesh):
    if (heat_mesh.faces[0], heat_mesh.faces[-1]) != (bar_mesh.nodes[0], bar_mesh.nodes[-1]):
        raise ValueError('thermal and mechanical domains must have identical endpoints')
    for i in range(1, len(bar_mesh.materials)):
        changed = (bar_mesh.materials[i] != bar_mesh.materials[i-1]
                   or bar_mesh.areas[i] != bar_mesh.areas[i-1])
        if changed and bar_mesh.nodes[i] not in heat_mesh.faces:
            raise ValueError('mechanical property interfaces must be thermal faces')


def overlap_temperatures(heat_mesh, temperature, bar_mesh, *, reference_temperature=0.):
    """Return mechanical element means of T-Tref (K).

    For overlap length w[e,i], mean[e] = sum_i w[e,i]*(T[i]-Tref)/h[e].
    Since alpha is constant per mechanical element and its interfaces are heat
    faces, sum_e alpha[e]*h[e]*mean[e] preserves the source eigenstrain integral.
    This preserves the supplied cell field, not an unknown continuum solution.
    """
    validate_meshes(heat_mesh, bar_mesh)
    t = _field(heat_mesh, temperature)
    if not np.isfinite(reference_temperature):
        raise ValueError('reference temperature must be finite')
    x, y = np.asarray(bar_mesh.nodes), np.asarray(heat_mesh.faces)
    overlap = np.maximum(0., np.minimum(x[1:, None], y[None, 1:])
                         - np.maximum(x[:-1, None], y[None, :-1]))
    return overlap @ (t-reference_temperature) / np.diff(x)
