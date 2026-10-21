"""Static and modal analysis of straight, uniform Timoshenko cantilevers."""

from .model import Beam, assemble, exact_tip, modes, solve_static

__all__ = ["Beam", "assemble", "exact_tip", "modes", "solve_static"]
__version__ = "0.1.0"
