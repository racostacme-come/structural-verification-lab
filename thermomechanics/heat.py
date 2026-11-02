"""Original cell-centred 1D conduction, constant area, SI units.

Material interfaces must be cell faces. Contact resistance is areal (m² K/W).
None at a boundary means insulation; a finite number is prescribed temperature.
Dense backward Euler is intended for small verification problems.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class HeatMesh:
    faces: tuple[float, ...]
    conductivity: tuple[float, ...]  # W/(m K)
    capacity: tuple[float, ...]  # volumetric rho*c, J/(m³ K)
    contact: tuple[float, ...]  # one nonnegative resistance per interior face
    area: float = 1.0  # m², constant

    def __post_init__(self):
        for name in ('faces', 'conductivity', 'capacity', 'contact'):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        x = np.asarray(self.faces, dtype=float)
        if x.ndim != 1 or len(x) < 2 or not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
            raise ValueError('faces must be finite and strictly increasing')
        n = len(x)-1
        for name, size, positive in [('conductivity', n, True), ('capacity', n, True), ('contact', n-1, False)]:
            a = np.asarray(getattr(self, name), dtype=float)
            if a.shape != (size,) or not np.all(np.isfinite(a)) or np.any(a <= 0 if positive else a < 0):
                raise ValueError(f'invalid {name}')
        if not np.isfinite(self.area) or self.area <= 0:
            raise ValueError('area must be finite and positive')

    @property
    def mass(self):
        return self.area*np.diff(self.faces)*np.asarray(self.capacity)

    @property
    def conductance(self):
        half = np.diff(self.faces)/(2*np.asarray(self.conductivity))
        return self.area/np.r_[half[0], half[:-1]+np.asarray(self.contact)+half[1:], half[-1]]


def _field(mesh, value):
    a = np.asarray(value, dtype=float)
    if a.shape != (len(mesh.capacity),) or not np.all(np.isfinite(a)):
        raise ValueError('one finite value per cell required')
    return a


def operator(mesh, *, left=None, right=None):
    """Return H and b in M dT/dt + H T = b (W)."""
    n = len(mesh.capacity)
    h, b = np.zeros((n, n)), np.zeros(n)
    g = mesh.conductance
    for i in range(n-1):
        h[i:i+2, i:i+2] += g[i+1]*np.array([[1., -1.], [-1., 1.]])
    for i, v, conductance in [(0, left, g[0]), (n-1, right, g[-1])]:
        if v is not None:
            if not np.isfinite(v):
                raise ValueError('boundary temperature must be finite or None')
            h[i, i] += conductance
            b[i] += conductance*v
    return h, b


def fluxes(mesh, temperature, *, left=None, right=None):
    """Face heat rates (W), positive toward increasing x, including boundaries."""
    t = _field(mesh, temperature)
    operator(mesh, left=left, right=right)  # validate boundary contract
    g = mesh.conductance
    return np.r_[0. if left is None else g[0]*(left-t[0]),
                 g[1:-1]*(t[:-1]-t[1:]),
                 0. if right is None else g[-1]*(t[-1]-right)]


def step(mesh, temperature, dt, *, left=None, right=None, source=None):
    """Backward Euler with end-of-step boundary values and cell source W/m³."""
    t = _field(mesh, temperature)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('dt must be finite and positive')
    power = np.zeros_like(t) if source is None else _field(mesh, source)*mesh.area*np.diff(mesh.faces)
    h, b = operator(mesh, left=left, right=right)
    mass_rate = mesh.mass/dt
    return np.linalg.solve(np.diag(mass_rate)+h, mass_rate*t+b+power)


def steady(mesh, *, left=None, right=None, source=None):
    """Steady solve with at least one prescribed temperature (unique solution)."""
    if left is None and right is None:
        raise ValueError('steady temperature needs at least one prescribed boundary')
    h, b = operator(mesh, left=left, right=right)
    if source is not None:
        b += _field(mesh, source)*mesh.area*np.diff(mesh.faces)
    return np.linalg.solve(h, b)
