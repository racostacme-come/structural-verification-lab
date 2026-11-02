"""Continuum oracle from axial equilibrium and integrated compatibility.

This module does not call the FE assembly/solver or consume an FE mesh.
Each layer has constant E,A,alpha and an affine prescribed temperature.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import Material


@dataclass(frozen=True)
class Layer:
    length: float
    area: float
    material: Material
    temperature_left: float  # increment relative to stress-free reference, K
    temperature_right: float

    def __post_init__(self):
        if not isinstance(self.material, Material):
            raise ValueError("material must be a Material")
        if not all(math.isfinite(v) for v in (self.length, self.area, self.temperature_left, self.temperature_right)):
            raise ValueError("layer values must be finite")
        if self.length <= 0 or self.area <= 0:
            raise ValueError("length and area must be positive")


@dataclass(frozen=True)
class ContinuumSolution:
    layers: tuple[Layer, ...]
    force: float  # tension positive, N
    left_displacement: float

    @property
    def reactions(self):
        return (-self.force, self.force)

    def displacement(self, positions):
        x = np.asarray(positions, dtype=float)
        length = math.fsum(layer.length for layer in self.layers)
        if not np.all(np.isfinite(x)) or np.any(x < 0) or np.any(x > length):
            raise ValueError("reference coordinates must lie in [0, total length]")
        u = np.full_like(x, self.left_displacement)
        start = 0.0
        for layer in self.layers:
            s = np.clip(x - start, 0.0, layer.length)
            thermal_integral = layer.material.expansion * (
                layer.temperature_left * s
                + (layer.temperature_right - layer.temperature_left) * s*s / (2*layer.length))
            u += self.force * s / (layer.material.young * layer.area) + thermal_integral
            start += layer.length
        return u


def continuum_solution(layers, *, left=0.0, right=0.0):
    """Fixed end displacements; right=None means a traction-free right end.

    With no distributed axial load, N is constant and
    N = (u_R-u_L - integral(alpha*deltaT dx)) / integral(1/(EA) dx).
    """
    layers = tuple(layers)
    if not layers or any(not isinstance(layer, Layer) for layer in layers):
        raise ValueError("at least one Layer is required")
    if not math.isfinite(left) or (right is not None and not math.isfinite(right)):
        raise ValueError("prescribed displacements must be finite")
    thermal = math.fsum(layer.material.expansion * layer.length *
                        (layer.temperature_left + layer.temperature_right)/2 for layer in layers)
    compliance = math.fsum(layer.length/(layer.material.young*layer.area) for layer in layers)
    force = 0.0 if right is None else (right-left-thermal)/compliance
    return ContinuumSolution(layers, force, left)
