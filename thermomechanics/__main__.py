"""Run original prescribed-temperature experiments: python -m thermomechanics."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from . import BarMesh, Material, Layer, continuum_solution, solve


def experiment():
    soft = Material(70e9, 23e-6)
    stiff = Material(190e9, 11e-6)
    cases = []
    for name, materials, areas, dt, right in [
        ("uniform_restrained", (soft, soft), (2e-4, 2e-4), (40., 40.), 0.),
        ("uniform_free", (soft, soft), (2e-4, 2e-4), (40., 40.), None),
        ("two_material_restrained", (soft, stiff), (2e-4, 1.5e-4), (60., 25.), 0.),
        ("two_material_imposed", (soft, stiff), (2e-4, 1.5e-4), (60., 25.), 1e-4),
    ]:
        mesh = BarMesh((0., .4, 1.), materials, areas)
        layers = [Layer(h,a,m,t,t) for h,a,m,t in zip(np.diff(mesh.nodes),areas,materials,dt)]
        reference = continuum_solution(layers, right=right)
        result = solve(mesh, dt, right=right)
        force_scale = max(1., abs(reference.force), *(m.young*a*abs(m.expansion*t) for m,a,t in zip(materials,areas,dt)))
        displacement_scale = max(1e-12, max(abs(reference.displacement(mesh.nodes))),
                                 sum(abs(m.expansion*t)*h for m,t,h in zip(materials,dt,np.diff(mesh.nodes))))
        cases.append({"name": name, "nodes_m": list(mesh.nodes),
                      "young_Pa": [m.young for m in materials], "expansion_per_K": [m.expansion for m in materials],
                      "areas_m2": list(areas), "mean_temperature_increment_K": list(dt),
                      "right_displacement_m": right, "reference_force_N": reference.force,
                      "element_force_N": result.element_force.tolist(),
                      "displacement_m": result.displacement.tolist(),
                      "support_force_N": result.residual[[0,-1]].tolist(),
                      "force_scale_N": force_scale, "displacement_scale_m": displacement_scale,
                      "scaled_force_error": float(np.max(abs(result.element_force-reference.force))/force_scale),
                      "scaled_nodal_displacement_error": float(np.max(abs(result.displacement-reference.displacement(mesh.nodes)))/displacement_scale)})
    # Smooth affine heating: exact element averages yield exact nodal values,
    # but the P1 displacement between nodes still has O(h^2) continuum error.
    reference = continuum_solution([Layer(1., 2e-4, soft, 0., 80.)])
    refinement = []
    for n in (4, 8, 16, 32):
        x = np.linspace(0, 1, n+1)
        mid = (x[:-1]+x[1:])/2
        mesh = BarMesh(x, (soft,)*n, (2e-4,)*n)
        result = solve(mesh, 80*mid)
        error = float(np.max(abs((result.displacement[:-1]+result.displacement[1:])/2-reference.displacement(mid))))
        refinement.append({"elements": n, "midpoint_displacement_error_m": error,
                           "max_nodal_error_m": float(np.max(abs(result.displacement-reference.displacement(x)))),
                           "observed_order": None if not refinement else float(np.log2(refinement[-1]["midpoint_displacement_error_m"]/error))})
    root = Path(__file__).resolve().parent
    return {"schema_version": 1, "scope": "Prescribed temperature; no heat solver or mesh transfer yet; synthetic verification only.",
            "source_sha256": {p.name: hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in sorted(root.glob('*.py'))},
            "cases": cases, "affine_temperature_refinement": refinement}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="JSON destination (stdout if omitted)")
    args = parser.parse_args()
    data = json.dumps(experiment(), indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(data, encoding="utf-8", newline="\n")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
