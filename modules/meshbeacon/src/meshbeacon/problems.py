"""Original manufactured solutions; no course datasets are used."""

from dataclasses import dataclass

import numpy as np

from .core import Field, _positive


@dataclass(frozen=True)
class Problem:
    name: str
    exact: Field
    gradient: Field
    load: Field


def sine(diffusion: float = 1.0) -> Problem:
    a = _positive(diffusion, "diffusion")
    return Problem(
        "sine",
        lambda x: np.sin(np.pi * x),
        lambda x: np.pi * np.cos(np.pi * x),
        lambda x: a * np.pi**2 * np.sin(np.pi * x),
    )


def bump(center: float = 0.37, width: float = 0.04, diffusion: float = 1.0) -> Problem:
    """Gaussian displacement with a linear correction for homogeneous boundaries."""
    a = _positive(diffusion, "diffusion")
    _positive(width, "width")
    if not np.isfinite(center) or not 0 < center < 1:
        raise ValueError("center must lie in (0, 1)")
    left, right = np.exp(-((np.array([0 - center, 1 - center]) / width) ** 2))

    def exact(x):
        return np.exp(-(((x - center) / width) ** 2)) - ((1 - x) * left + x * right)

    def gradient(x):
        z = (x - center) / width
        return -2 * z / width * np.exp(-(z**2)) + left - right

    def load(x):
        z = (x - center) / width
        return a * (2 - 4 * z**2) / width**2 * np.exp(-(z**2))

    return Problem("bump", exact, gradient, load)
