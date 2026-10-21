"""A common SI cantilever connects Timoshenko compliance to Euler buckling."""
import numpy as np
from shearscope.model import Beam, solve_static, exact_tip
from bucklemap.model import solve

def experiment():
    beam = Beam(length=1.0, width=0.02, height=0.01)
    reference = np.pi**2 * beam.ei / (4 * beam.length**2)
    rows = []
    for n in (4, 8, 16, 32):
        load = float(solve(np.linspace(0, beam.length, n+1), beam.ei, support="cantilever").loads[0])
        reduced = solve_static(beam, n, "reduced").displacement[-1]
        full = solve_static(beam, n, "full").displacement[-1]
        # Subtract shear flexibility: Euler buckling and Timoshenko bending
        # refer to the same EI, L but not the same shear-deformable theory.
        bending = reduced - beam.length / beam.shear
        rows.append(dict(elements=n, buckling_relative_error=abs(load/reference-1),
                         reduced_compliance_relative_error=abs(reduced/exact_tip(beam)-1),
                         full_compliance_relative_error=abs(full/exact_tip(beam)-1),
                         dimensionless_product=float(load*bending/beam.length)))
    return {"experiment": "common-cantilever", "units": "SI", "rows": rows,
            "product_reference": float(np.pi**2/12),
            "identity": "P_Euler * (C_Timoshenko - L/(kGA)) / L = pi^2/12",
            "scope": "Linear elastic slender cantilever; no nonlinear stability or physical validation."}
