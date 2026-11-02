"""Seeded layered-bar verification of native assembly, not physical validation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from . import BarMesh, Material, Layer, continuum_solution, assemble, solve
from .native import NativeAssembler


def experiment(library=None):
    backend = NativeAssembler(library)
    rng = np.random.default_rng(20261029)
    rows = []
    for n in (1, 2, 7, 24):
        for sample in range(3):
            lengths = rng.uniform(.02, .3, n)
            areas = rng.uniform(1e-4, 8e-4, n)
            materials = tuple(Material(e, a) for e, a in
                              zip(rng.uniform(20e9, 210e9, n), rng.uniform(-8e-6, 25e-6, n)))
            tl, tr = rng.uniform(-40, 120, (2, n))
            layers = [Layer(h, a, m, l, r) for h, a, m, l, r in zip(lengths, areas, materials, tl, tr)]
            x = np.r_[0., np.cumsum(lengths)]
            x[-1] = math.fsum(lengths)
            mesh = BarMesh(x, materials, areas)
            dt = (tl+tr)/2
            kp, fp = assemble(mesh, dt)
            kn, fn = backend.assemble(mesh, dt)
            for right in (None, 0., 2e-4):
                left = -1e-4
                ref = continuum_solution(layers, left=left, right=right)
                result = solve(mesh, dt, left=left, right=right, assembler=backend.assemble)
                force_scale = max(1., abs(ref.force), float(np.max(np.abs(fp))))
                displacement_scale = max(1e-12, float(np.max(np.abs(ref.displacement(x)))))
                free = result.residual[1:] if right is None else result.residual[1:-1]
                rows.append(dict(elements=n, sample=sample, right_displacement_m=right,
                                 scaled_stiffness_difference=float(np.max(abs(kn-kp))/np.max(abs(kp))),
                                 scaled_load_difference=float(np.max(abs(fn-fp)))/force_scale,
                                 scaled_force_error=float(np.max(abs(result.element_force-ref.force)))/force_scale,
                                 scaled_nodal_error=float(np.max(abs(result.displacement-ref.displacement(x))))/displacement_scale,
                                 scaled_free_residual=float(np.max(abs(free), initial=0.))/force_scale,
                                 scaled_support_error=float(np.max(abs(result.residual[[0,-1]]-ref.reactions)))/force_scale))
    root = Path(__file__).resolve().parents[1]
    sources = sorted((root/'thermomechanics').glob('*.py')) + sorted((root/'native').glob('*'))
    return dict(schema_version=1, seed=20261029,
                scope='C++ assembly consistency and independent prescribed-temperature continuum checks; shared NumPy solve; no heat/transfer error or physical validation',
                source_sha256={p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in sources if p.is_file()},
                cases=rows, maximum_scaled_errors={key: max(row[key] for row in rows) for key in rows[0] if key.startswith('scaled_')})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path)
    parser.add_argument('--output', type=Path, default=Path('results/native-thermoelastic.json'))
    args = parser.parse_args()
    data = experiment(args.library)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(data['maximum_scaled_errors'], indent=2))
