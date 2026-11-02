"""SI contracts. No imports from the preserved component snapshots."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Material:
    young: float  # Pa
    expansion: float  # 1/K; signed (zero and negative are admissible)

    def __post_init__(self):
        if not np.isfinite(self.young) or self.young <= 0:
            raise ValueError("young must be finite and positive")
        if not np.isfinite(self.expansion):
            raise ValueError("expansion must be finite")


@dataclass(frozen=True)
class BarMesh:
    """P1 mesh; every material/area interface must coincide with a node.

    Immutable tuples prevent caller mutation after validation. Areas are m^2,
    coordinates m. One material and one area are required for each element.
    """
    nodes: tuple[float, ...]
    materials: tuple[Material, ...]
    areas: tuple[float, ...]

    def __post_init__(self):
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(self, "materials", tuple(self.materials))
        object.__setattr__(self, "areas", tuple(self.areas))
        x = np.asarray(self.nodes, dtype=float)
        a = np.asarray(self.areas, dtype=float)
        if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
            raise ValueError("nodes must be finite and strictly increasing")
        if len(self.materials) != x.size - 1 or any(not isinstance(m, Material) for m in self.materials):
            raise ValueError("one Material is required per element")
        if a.shape != (x.size - 1,) or not np.all(np.isfinite(a)) or np.any(a <= 0):
            raise ValueError("one finite positive area is required per element")


def element_temperatures(mesh, delta_temperature):
    """Mean T-Tref in K for each element, not nodal temperatures."""
    values = np.asarray(delta_temperature, dtype=float)
    if values.shape != (len(mesh.materials),) or not np.all(np.isfinite(values)):
        raise ValueError("one finite mean temperature increment is required per element")
    return values
