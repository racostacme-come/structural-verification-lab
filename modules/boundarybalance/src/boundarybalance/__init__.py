"""Weak Dirichlet finite elements and numerical verification."""

__version__ = "0.1.0"

from .fem import assemble, critical_penalty, errors, solve, unit_square

__all__ = ["assemble", "critical_penalty", "errors", "solve", "unit_square"]
