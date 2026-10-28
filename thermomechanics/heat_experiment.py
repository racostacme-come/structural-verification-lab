"""Synthetic contact sweep and separate temporal/spatial verification studies."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .heat import HeatMesh, fluxes, step, steady


def experiment():
    contact_rows = []
    for resistance in [0., 1e-10, .3, 1e5]:
        m = HeatMesh((0., .2, .55, 1.), (2., 8., 8.), (3., 4., 4.), (resistance, 0.), .7)
        t = steady(m, left=12., right=2.)
        exact = .7*10/(.2/2+.8/8+resistance)
        q = fluxes(m, t, left=12., right=2.)
        transient = step(m, [3., 9., 5.], .1, left=12., right=2.)
        rates = fluxes(m, transient, left=12., right=2.)
        balance = float(np.dot(m.mass, transient-np.array([3., 9., 5.]))-.1*(rates[0]-rates[-1]))
        contact_rows.append(dict(resistance_m2_K_per_W=resistance, reference_rate_W=exact,
                                 max_relative_rate_error=float(max(abs(q-exact))/abs(exact)),
                                 transient_energy_residual_J=balance))
    temporal, spatial = [], []
    n = 12
    m = HeatMesh(tuple(np.linspace(0, 1, n+1)), (1.,)*n, (1.,)*n, (0.,)*(n-1))
    mode = np.cos(np.pi*(np.arange(n)+.5)/n)
    lam = 4*n*n*np.sin(np.pi/(2*n))**2
    for count in [10, 20, 40, 80]:
        dt = .1/count
        t = mode.copy()
        for _ in range(count):
            t = step(m, t, dt)
        temporal.append(dict(steps=count, dt_s=dt,
                             error_K=float(max(abs(t-mode*np.exp(-lam*.1))))))
    # Exact semidiscrete evolution removes time error; compare actual cell averages.
    for n in [8, 16, 32, 64]:
        mode = np.cos(np.pi*(np.arange(n)+.5)/n)*np.sinc(1/(2*n))
        lam = 4*n*n*np.sin(np.pi/(2*n))**2
        spatial.append(dict(cells=n, error_K=float(max(abs(mode*(np.exp(-lam*.1)-np.exp(-np.pi**2*.1)))))))
    for rows in [temporal, spatial]:
        for a,b in zip(rows,rows[1:]):
            b['observed_order'] = float(np.log2(a['error_K']/b['error_K']))
    root = Path(__file__).resolve().parents[1]
    names = ['thermomechanics/heat.py', 'thermomechanics/heat_experiment.py']
    return dict(scope='Synthetic verification only; no coupled transfer or physical validation',
                source_sha256={name: hashlib.sha256((root/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names},
                contact=contact_rows, temporal_against_semidiscrete=temporal,
                spatial_exact_semidiscrete_against_continuum_cell_average=spatial)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/thermal-contact.json'))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(experiment(), indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
