"""Reproducible synthetic nonmatching transfer and transient contact experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .model import BarMesh, Material
from .heat import HeatMesh, fluxes
from .reference import Layer, continuum_solution
from .mechanics import solve
from .native import NativeAssembler
from .transfer import overlap_temperatures
from .coupled import advance


MATERIALS = (Material(70e9, 23e-6), Material(190e9, 11e-6))
AREA = 2e-4
INTERFACE = .4
TREF = 293.15


def grids(nleft, nright, resistance=0.):
    faces = np.r_[np.linspace(0., INTERFACE, 8), np.linspace(INTERFACE, 1., 12)[1:]]
    contact = np.zeros(17)
    contact[6] = resistance
    thermal = HeatMesh(faces, (13.,)*7+(4.,)*11, (2e5,)*7+(3e5,)*11, contact, AREA)
    nodes = np.r_[np.linspace(0., INTERFACE, nleft+1), np.linspace(INTERFACE, 1., nright+1)[1:]]
    mechanical = BarMesh(nodes, (MATERIALS[0],)*nleft+(MATERIALS[1],)*nright,
                         (AREA,)*(nleft+nright))
    return thermal, mechanical


def source_reference(heat_mesh, temperature):
    """Integrate original heat cells; never consume transferred element means.

    Specific to this experiment: interface=.4, constant area, zero supports.
    Exact continuum mechanics for the piecewise-constant numerical temperature.
    """
    layers = []
    for l, r, t in zip(heat_mesh.faces[:-1], heat_mesh.faces[1:], temperature, strict=True):
        material = MATERIALS[0] if (l+r)/2 < INTERFACE else MATERIALS[1]
        layers.append(Layer(r-l, AREA, material, t-TREF, t-TREF))
    return continuum_solution(layers)


def naive_means(thermal, temperature, mechanical):
    """Deliberately nonconservative baseline: linear interpolation at midpoints.

    numpy.interp bridges thermal/contact interfaces and clamps at end centres.
    It changes the reconstructed field; it is not an integral-preserving map.
    """
    centres = (np.array(thermal.faces[:-1])+thermal.faces[1:])/2
    midpoints = (np.array(mechanical.nodes[:-1])+mechanical.nodes[1:])/2
    return np.interp(midpoints, centres, temperature)-TREF


def experiment(library=None):
    backend = NativeAssembler(library)
    transfer_rows = []
    for level in (2, 4, 8, 16):
        thermal, mechanical = grids(level, level+1)
        centres = (np.array(thermal.faces[:-1])+thermal.faces[1:])/2
        # Frozen source field with a jump. Only the target grid is refined.
        temperature = TREF + np.where(centres < INTERFACE, 80.-50.*centres, 10.+20.*centres)
        reference = source_reference(thermal, temperature)
        means = overlap_temperatures(thermal, temperature, mechanical, reference_temperature=TREF)
        conservative = solve(mechanical, means, assembler=backend.assemble)
        naive = solve(mechanical, naive_means(thermal, temperature, mechanical))
        alpha = np.array([m.expansion for m in mechanical.materials])
        mapped_integral = float(np.dot(alpha*np.diff(mechanical.nodes), means))
        source_integral = sum(layer.material.expansion*layer.length*layer.temperature_left
                              for layer in reference.layers)
        scale = max(1e-12, sum(abs(layer.material.expansion*layer.length*layer.temperature_left)
                               for layer in reference.layers))
        transfer_rows.append(dict(elements=len(mechanical.materials),
                                  source_eigenstrain_integral_m=source_integral,
                                  mapped_eigenstrain_integral_m=mapped_integral,
                                  scaled_integral_error=abs(mapped_integral-source_integral)/scale,
                                  reference_force_N=reference.force,
                                  overlap_force_N=float(conservative.element_force[0]),
                                  naive_force_N=float(naive.element_force[0]),
                                  overlap_relative_force_error=float(max(abs(conservative.element_force-reference.force))/abs(reference.force)),
                                  naive_relative_force_error=float(max(abs(naive.element_force-reference.force))/abs(reference.force))))
    histories = []
    for resistance in (0., .02, 2.):
        thermal, mechanical = grids(3, 5, resistance)
        temperature = np.full(18, TREF)
        rows = []
        for iteration in range(1, 21):
            state = advance(thermal, mechanical, temperature, 500., reference_temperature=TREF,
                            thermal_left=TREF+80., thermal_right=TREF, assembler=backend.assemble)
            reference = source_reference(thermal, state.temperature)
            python = solve(mechanical, state.mean_temperature_increment)
            naive = solve(mechanical, naive_means(thermal, state.temperature, mechanical))
            q = fluxes(thermal, state.temperature, left=TREF+80., right=TREF)
            change = float(thermal.mass @ (state.temperature-temperature))
            input_energy = 500.*(q[0]-q[-1])
            force_scale = max(1., abs(reference.force))
            displacement_scale = max(1e-12, float(max(abs(reference.displacement(mechanical.nodes)))))
            rows.append(dict(time_s=iteration*500., temperature_K=state.temperature.tolist(),
                             displacement_m=state.mechanics.displacement.tolist(),
                             reference_force_N=reference.force,
                             overlap_force_N=float(state.mechanics.element_force[0]),
                             naive_force_N=float(naive.element_force[0]),
                             energy_change_J=change, boundary_energy_J=float(input_energy),
                             scaled_energy_residual=abs(change-input_energy)/max(1., abs(change), abs(input_energy)),
                             scaled_force_error=float(max(abs(state.mechanics.element_force-reference.force)))/force_scale,
                             scaled_nodal_error=float(max(abs(state.mechanics.displacement-reference.displacement(mechanical.nodes))))/displacement_scale,
                             scaled_native_python_difference=float(max(abs(state.mechanics.element_force-python.element_force)))/force_scale,
                             scaled_free_residual=float(max(abs(state.mechanics.residual[1:-1])))/force_scale))
            temperature = state.temperature
        histories.append(dict(contact_resistance_m2_K_per_W=resistance, heat_faces_m=list(thermal.faces),
                              mechanical_nodes_m=list(mechanical.nodes), steps=rows))
    root = Path(__file__).resolve().parents[1]
    sources = sorted((root/'thermomechanics').glob('*.py')) + sorted((root/'native').glob('*'))
    return dict(schema_version=1,
                scope='Supplied-cell-field transfer and mechanics verification, plus discrete energy balance; transient continuum heat accuracy and physical validity are not established here.',
                parameters=dict(interface_m=INTERFACE, area_m2=AREA, reference_temperature_K=TREF,
                                young_Pa=[m.young for m in MATERIALS], expansion_per_K=[m.expansion for m in MATERIALS],
                                conductivity_W_per_m_K=[13., 4.], capacity_J_per_m3_K=[2e5, 3e5],
                                initial_temperature_K=TREF, boundary_temperature_K=[TREF+80., TREF],
                                end_displacement_m=[0., 0.], time_step_s=500., steps=20),
                source_sha256={p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                               for p in sources if p.is_file()},
                transfer_only=transfer_rows, coupled_histories=histories,
                maximum_scaled_errors={key: max(row[key] for history in histories for row in history['steps'])
                                       for key in rows[0] if key.startswith('scaled_')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path, help='Explicit C++ library; otherwise use suite build record')
    parser.add_argument('--output', type=Path, default=Path('results/coupled-transfer.json'))
    args = parser.parse_args()
    data = experiment(args.library)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(data['maximum_scaled_errors'], indent=2))


if __name__ == '__main__':
    main()
