"""Separate coupled temporal, continuum spatial and contact/property studies.

All inputs are original synthetic SI verification cases. Native mechanics
assembly is exercised; dense heat evolution/linear solves remain in NumPy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .heat import HeatMesh, operator, steady, fluxes
from .model import BarMesh, Material
from .mechanics import solve
from .native import NativeAssembler
from .coupled import advance
from .transfer import overlap_temperatures
from .reference import Layer, continuum_solution
from .coupled_experiment import AREA, INTERFACE, MATERIALS, grids, source_reference


AMPLITUDE = 40.  # K; temperature increments relative to a stress-free state
END_TIME = .1  # s; synthetic thermal diffusivity = 1 m2/s


def cosine_reference(positions, amplitude=AMPLITUDE):
    """Exact continuum mechanics for deltaT=amplitude*cos(pi*x), x in [0,1].

    Integrates alpha*T and compliance analytically on the two materials.
    Does not use a heat matrix, transfer, FE assembly or numerical quadrature.
    """
    x = np.asarray(positions, dtype=float)
    if not np.isfinite(amplitude) or not np.all(np.isfinite(x)) or np.any((x < 0) | (x > 1)):
        raise ValueError('finite amplitude and coordinates in [0, 1] required')
    m0, m1 = MATERIALS
    compliance = INTERFACE/(m0.young*AREA)+(1-INTERFACE)/(m1.young*AREA)
    integral = amplitude*(m0.expansion-m1.expansion)*np.sin(np.pi*INTERFACE)/np.pi
    force = -integral/compliance
    split = np.minimum(x, INTERFACE)
    thermal = amplitude/np.pi*(m0.expansion*np.sin(np.pi*split)
                              + m1.expansion*(np.sin(np.pi*x)-np.sin(np.pi*split)))
    displacement = force*(split/(m0.young*AREA)+(x-split)/(m1.young*AREA))+thermal
    return float(force), displacement


def mode_mesh(cells):
    # Multiples of ten resolve the mechanical material interface at x=.4.
    if not isinstance(cells, int) or cells < 10 or cells % 10:
        raise ValueError('cells must be a positive multiple of ten')
    faces = np.linspace(0., 1., cells+1)
    heat = HeatMesh(faces, (1.,)*cells, (1.,)*cells, (0.,)*(cells-1), AREA)
    initial = AMPLITUDE*np.diff(np.sin(np.pi*faces))/(np.pi*np.diff(faces))
    return heat, initial


def spectral_temperature(mesh, initial, time):
    """Exact semidiscrete insulated evolution via a symmetric eigensolve.

    Independent of backward Euler; shares only the spatial heat operator.
    Mass similarity transform also supports nonuniform positive capacities.
    """
    h, _ = operator(mesh)
    root_mass = np.sqrt(mesh.mass)
    eigenvalues, vectors = np.linalg.eigh(h/root_mass[:, None]/root_mass[None, :])
    return (vectors @ (np.exp(-eigenvalues*time)*(vectors.T @ (root_mass*initial))))/root_mass


def _orders(rows):
    for previous, current in zip(rows, rows[1:]):
        current['orders'] = {key: float(np.log2(previous[key]/current[key]))
                             for key in ('temperature_error_K', 'force_error_N', 'nodal_error_m')}


def refinements(backend):
    _, bar = grids(3, 5)
    heat, initial = mode_mesh(20)
    exact_discrete = spectral_temperature(heat, initial, END_TIME)
    # source_reference expects absolute Kelvin, whereas this study evolves increments.
    from .coupled_experiment import TREF
    discrete_mechanics = source_reference(heat, exact_discrete+TREF)
    temporal = []
    for steps in (10, 20, 40, 80):
        temperature = initial.copy()
        dt = END_TIME/steps
        for _ in range(steps):
            state = advance(heat, bar, temperature, dt, assembler=backend.assemble)
            temperature = state.temperature
        temporal.append(dict(steps=steps, dt_s=dt,
                             temperature_error_K=float(max(abs(temperature-exact_discrete))),
                             force_error_N=float(max(abs(state.mechanics.element_force-discrete_mechanics.force))),
                             nodal_error_m=float(max(abs(state.mechanics.displacement-discrete_mechanics.displacement(bar.nodes))))))
    spatial = []
    exact_force, exact_u = cosine_reference(bar.nodes, AMPLITUDE*np.exp(-np.pi**2*END_TIME))
    for cells in (10, 20, 40, 80):
        heat, initial = mode_mesh(cells)
        temperature = spectral_temperature(heat, initial, END_TIME)
        continuum_temperature = initial*np.exp(-np.pi**2*END_TIME)
        means = overlap_temperatures(heat, temperature, bar)
        result = solve(bar, means, assembler=backend.assemble)
        eigenvalue = 4*cells*cells*np.sin(np.pi/(2*cells))**2
        spatial.append(dict(cells=cells, h_m=1/cells,
                            temperature_error_K=float(max(abs(temperature-continuum_temperature))),
                            force_error_N=float(max(abs(result.element_force-exact_force))),
                            nodal_error_m=float(max(abs(result.displacement-exact_u))),
                            spectral_mode_residual_K=float(max(abs(temperature-initial*np.exp(-eigenvalue*END_TIME))))))
    _orders(temporal)
    _orders(spatial)
    return temporal, spatial


def contact_sweep(backend):
    """Exact series-resistance heat field and affine-layer continuum mechanics.

    Fixed 7+11 heat cells, 3+5 mechanical elements. Contact changes temperature
    jump only: mechanical bonding remains perfect. k and E contrasts vary
    independently; no empirical material identification is implied.
    """
    rows = []
    for k_ratio in (.1, 1., 10.):
        for e_ratio in (.1, 1., 10.):
            for resistance in (0., 1e-10, .02, 2., 1e5):
                original_heat, original_bar = grids(3, 5, resistance)
                k0, k1 = 13., 13.*k_ratio
                materials = (MATERIALS[0], Material(MATERIALS[0].young*e_ratio, MATERIALS[1].expansion))
                heat = HeatMesh(original_heat.faces, (k0,)*7+(k1,)*11,
                                original_heat.capacity, original_heat.contact, AREA)
                bar = BarMesh(original_bar.nodes, (materials[0],)*3+(materials[1],)*5, (AREA,)*8)
                # Heat rates are W; resistance and flux density use area-normalized units.
                flux_density = 80./(INTERFACE/k0+resistance+(1-INTERFACE)/k1)
                left_interface = 80.-flux_density*INTERFACE/k0
                right_interface = flux_density*(1-INTERFACE)/k1
                centres = (np.array(heat.faces[:-1])+heat.faces[1:])/2
                exact_temperature = np.where(centres < INTERFACE,
                                             80.-flux_density*centres/k0,
                                             flux_density*(1-centres)/k1)
                reference = continuum_solution((
                    Layer(INTERFACE, AREA, materials[0], 80., left_interface),
                    Layer(1-INTERFACE, AREA, materials[1], right_interface, 0.)))
                temperature = steady(heat, left=80., right=0.)
                result = solve(bar, overlap_temperatures(heat, temperature, bar), assembler=backend.assemble)
                rates = fluxes(heat, temperature, left=80., right=0.)
                h, _ = operator(heat, left=80., right=0.)
                rows.append(dict(conductivity_ratio=k_ratio, young_ratio=e_ratio,
                                 resistance_m2_K_per_W=resistance,
                                 reference_rate_W=AREA*flux_density,
                                 reference_jump_K=flux_density*resistance,
                                 reference_force_N=reference.force,
                                 force_N=float(result.element_force[0]),
                                 relative_rate_error=float(max(abs(rates-AREA*flux_density))/(AREA*flux_density)),
                                 scaled_temperature_error=float(max(abs(temperature-exact_temperature))/80.),
                                 relative_force_error=float(max(abs(result.element_force-reference.force))/abs(reference.force)),
                                 nodal_error_m=float(max(abs(result.displacement-reference.displacement(bar.nodes)))),
                                 heat_condition_number=float(np.linalg.cond(h))))
    return rows


def experiment(library=None):
    backend = NativeAssembler(library)
    temporal, spatial = refinements(backend)
    contact = contact_sweep(backend)
    root = Path(__file__).resolve().parents[1]
    sources = sorted((root/'thermomechanics').glob('*.py'))+sorted((root/'native').glob('*'))
    return dict(schema_version=1,
                scope='Independent coupled temporal and smooth-field spatial verification; heterogeneous steady contact sweep. No transient contact continuum convergence or physical validation claimed.',
                parameters=dict(area_m2=AREA, interface_m=INTERFACE, amplitude_K=AMPLITUDE,
                                end_time_s=END_TIME, thermal_conductivity_W_per_m_K=1.,
                                volumetric_capacity_J_per_m3_K=1., mechanical_nodes_m=list(grids(3, 5)[1].nodes),
                                young_Pa=[m.young for m in MATERIALS], expansion_per_K=[m.expansion for m in MATERIALS],
                                temperature_convention='increments from stress-free reference',
                                refinement_boundaries='insulated heat, fixed mechanical ends',
                                contact_boundaries_K=[80., 0.], contact_cells_per_layer=[7, 11],
                                contact_left_conductivity_W_per_m_K=13., contact_left_young_Pa=MATERIALS[0].young),
                source_sha256={p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                               for p in sources if p.is_file()},
                temporal_against_semidiscrete=temporal,
                spatial_exact_semidiscrete_against_continuum=spatial,
                contact_property_sweep=contact)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path)
    parser.add_argument('--output', type=Path, default=Path('results/coupled-refinement.json'))
    args = parser.parse_args()
    data = experiment(args.library)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({key: data[key] for key in ('temporal_against_semidiscrete', 'spatial_exact_semidiscrete_against_continuum')}, indent=2))


if __name__ == '__main__':
    main()
