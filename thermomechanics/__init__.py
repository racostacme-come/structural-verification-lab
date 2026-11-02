"""Original small-strain, one-dimensional thermoelastic verification tools."""
from .model import BarMesh, Material
from .reference import Layer, continuum_solution
from .mechanics import assemble, solve

__all__ = ["BarMesh", "Material", "Layer", "continuum_solution", "assemble", "solve"]
