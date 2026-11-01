# Nonmatching thermoelastic coupling

**1 November status:** The four-page [thermoelastic synthesis](../paper/paper.pdf)
and [source](../paper/paper.tex) now cover the coupled results. The text below
retains the milestone chronology; future-tense manuscript/refinement statements
describe that earlier date. Release remains subject to acceptance review.

30 October 2026, cycle day 4. Original synthetic implementation; the released
beam paper has not yet been replaced. No component snapshot changed.

## Transfer and assumptions

The heat field is piecewise constant on thermal cells I_i. Mechanical P1
elements J_e have constant E, A and alpha. With w_ei = |I_i intersect J_e|,

```
mean_deltaT[e] = sum_i w_ei * (T[i] - Tref) / length(J_e)
integrated_eigenstrain = sum_e alpha[e] * length(J_e) * mean_deltaT[e]
```

Both meshes must have exactly identical endpoints. Mechanical material and
area jumps must be thermal faces (identical coordinate values, without
fuzzy snapping). Thermal contact jumps may lie inside a mechanical element;
the overlap integral handles the jump without smoothing it. The caller is
responsible for representing physical conductivity/capacity interfaces as
thermal faces and mechanical properties as mechanical nodes. Property values
cannot establish a physical material identity automatically.

The coupling requires the same constant area in both meshes, and a single
finite reference temperature on the same scale as the temperature field.
It performs one backward-Euler heat step, integrates eigenstrain and solves
static linear elasticity. There is no deformation feedback, latent heat,
inertia, plasticity or mechanical contact opening. Dense O(n_heat*n_mech)
overlaps and dense solves are for small verification studies.

## API and reproducible CLI

After README setup and `python suite.py native`:

```sh
python -m thermomechanics.coupled_experiment --output results/coupled-transfer.json
python -m pytest tests/test_transfer.py -q
```

`--library PATH` selects a C++ library explicitly; missing native code fails
instead of silently falling back. The CLI uses C++17 assembly with NumPy
boundary elimination/solution. It runs four target-grid levels against one
frozen thermal field and three 20-step transient contact cases. JSON records
geometry, material properties, boundary/initial conditions, temperatures,
displacements, forces, energy and normalized-LF source hashes. Reproduction
tests use scaled tolerances across compilers. CI runs the CLI on Linux and
Windows and uploads its evidence separately from committed reference data.

```python
from thermomechanics.coupled import advance
from thermomechanics.native import NativeAssembler

state = advance(heat_mesh, bar_mesh, temperatures, dt,
                reference_temperature=293.15,
                thermal_left=373.15, thermal_right=293.15,
                displacement_left=0., displacement_right=0.,
                assembler=NativeAssembler().assemble)
```

The caller supplies validated `HeatMesh` and `BarMesh` objects (see `grids`
in the experiment for a complete construction), previous cell temperatures
and dt in seconds. Thermal boundaries and an optional source in W/m^3 are
end-of-step values; omitted thermal boundaries are insulated. Set
`displacement_right=None` for a free end. Omitting `assembler` uses Python.
`CoupledState` returns the new thermal field, transferred increments and
mechanical displacements/residual/element forces. The input field is unchanged.
`overlap_temperatures` is also available independently in `transfer.py`.

## Measured evidence and limits

On the frozen discontinuous field, the independent source-cell integral is
0.0008024 m and the restrained force is -18088 N. The maximum normalized
integral error is 4.06e-16; the four grids have 5, 9, 17 and 33 elements.
The midpoint linear-interpolation baseline has relative force errors near
roundoff, near roundoff, 0.003931 and 0.005267 respectively. This outcome is
not monotone convergence: early quadrature cancellation is lost as points
approach the smoothed jump. The baseline reconstructs a different field,
bridges material/contact interfaces and clamps outside the end cell centres.
We do not claim it always performs worse.

Across 60 coupled steps, maximum scaled source-field force and nodal errors
are 1.75e-15 and 1.61e-15. Native/Python force differences are zero locally;
scaled free-node residual is 2.00e-15. Maximum discrete energy residual is
2.11e-12, using scale max(1 J, |stored energy change|, |boundary input|).
All error values are platform-dependent floating-point measurements.

At 10000 s the three computed compressive forces are -19136.55, -18692.40
and -16690.30 N for areal contact resistances 0, 0.02 and 2 m^2 K/W. These
are model outputs at one thermal resolution and time step, not converged
continuum predictions. The mechanics oracle integrates the **original heat
cells**, not transferred means; its agreement establishes transfer/assembly
consistency and mechanics accuracy for that supplied field. Discrete energy
balance does not establish thermal continuum accuracy. Independent coupled
temporal/spatial studies are scheduled for 31 October. No physical validation
or performance advantage is claimed.

## Provenance and upcoming synthesis

The method follows direct integration of the constitutive weak form and
axial compatibility, derived in [the research contract](thermoelastic-transfer.md).
Topic provenance remains MMME heat/interface resistance and TUM finite elements.
New MMME notes inspected read-only on 30 October cover segregation/phase mass
balance and implicit/theta time stepping. Only topic descriptions are recorded;
no course files, exercise values, solutions or teaching code are published.
Phase-change modelling is a possible later Transport question, outside this cycle.
The <=5-page extension synthesis is planned for 1 November with separate
thermal, transfer and mechanics error tables and complete rendered-page review.
