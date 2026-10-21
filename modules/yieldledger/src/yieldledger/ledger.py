"""Separate stored energy, irreversible dissipation, and endpoint-work excess."""

import numpy as np

from ._core import Material, State


def ledger(material: Material, strains) -> dict[str, np.ndarray]:
    """Integrate a strain path starting at zero from a virgin state.

    Each entry contains cumulative right-endpoint work, stored energy,
    yield_stress * alpha (irreversible dissipation), and the exactly derived
    nonnegative work-quadrature defect. Strain reversals must be sampled.
    """
    strains = np.asarray(strains, dtype=float)
    if (
        strains.ndim != 1
        or not len(strains)
        or not np.all(np.isfinite(strains))
        or strains[0] != 0
    ):
        raise ValueError("Path must be a finite nonempty 1D array starting at zero.")
    names = (
        "strain",
        "stress",
        "plastic_strain",
        "alpha",
        "work",
        "stored",
        "dissipation",
        "defect",
    )
    data = {name: np.zeros(len(strains)) for name in names}
    data["strain"] = strains.copy()
    old = State()
    old_strain = work = defect = 0.0
    for i, strain in enumerate(strains):
        response = material.update(strain, old)
        new = response.state
        dep = new.plastic_strain - old.plastic_strain
        da = new.accumulated - old.accumulated
        dee = strain - old_strain - dep
        work += response.stress * (strain - old_strain)
        defect += (
            material.elastic * dee**2 + material.isotropic * da**2 + material.kinematic * dep**2
        ) / 2
        for name, value in zip(
            names[1:],
            (
                response.stress,
                new.plastic_strain,
                new.accumulated,
                work,
                material.energy(strain, new),
                material.yield_stress * new.accumulated,
                defect,
            ),
            strict=True,
        ):
            data[name][i] = value
        old, old_strain = new, strain
    return data
