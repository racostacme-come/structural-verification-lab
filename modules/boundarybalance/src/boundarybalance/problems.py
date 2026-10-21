"""Original smooth manufactured fields; no private data or course solutions."""

import numpy as np


def exact(points):
    x, y = points.T
    return 1 + x + 2 * y + np.sin(np.pi * x) * np.sin(np.pi * y)


def gradient(points):
    x, y = points.T
    return np.column_stack(
        (
            1 + np.pi * np.cos(np.pi * x) * np.sin(np.pi * y),
            2 + np.pi * np.sin(np.pi * x) * np.cos(np.pi * y),
        )
    )


def source(points):
    x, y = points.T
    return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
