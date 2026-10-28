# Conservative thermal contact baseline

Original implementation in `thermomechanics/heat.py`, motivated by the previously
reviewed MMME (Modelling in Materials and Manufacturing Engineering) topic of
resistance interfaces and heat positivity. All parameters below are synthetic;
no private teaching files, solutions or code are included.

For constant cross-sectional area A and cell width h, the heat capacity is
M_i = A h_i (rho c)_i. Interior conductance is
G = A / (h_i/(2 k_i) + R + h_j/(2 k_j)), where R is areal contact
resistance in m² K/W. Boundary conductance is 2 A k/h. The oriented face
rate is Q = G(T_left - T_right), in W. Finite-volume balance is
M_i dT_i/dt = Q_left - Q_right + A h_i s_i, with source s in W/m³.
Backward Euler solves (M/dt + H) T_new = M T_old/dt + b + A h s.
Boundary values and sources are evaluated at the end of the step.

`HeatMesh(faces, conductivity, capacity, contact, area)` uses immutable tuples;
one conductivity/capacity per cell and one contact resistance per interior face
are mandatory. Interfaces must coincide with cell faces. `step`, `steady` and
`fluxes` accept `left`/`right` temperatures; `None` means insulation. A steady
solve requires at least one temperature boundary. `fluxes` returns n+1 rates
positive in increasing x. Temperatures may be signed reference increments for
linear conduction; they are not interpreted as absolute thermodynamic values.

```python
from thermomechanics.heat import HeatMesh, step
mesh = HeatMesh((0., .2, .55, 1.), (2., 8., 8.), (3., 4., 4.), (.3, 0.), .7)
temperature = step(mesh, [3., 9., 5.], .1, left=12., right=2.)
```

Run `python -m thermomechanics.heat_experiment` to reproduce
`results/thermal-contact.json`; `--output PATH` supports independent artifacts.
CI runs the experiment on Windows and Linux. Run `python suite.py test` for
the isolated component regressions and root checks.

Verification separates: (1) assembly and time-step implementation against an
independent discrete cosine eigenvalue; (2) backward Euler temporal error
against exact semidiscrete evolution; (3) thermal spatial error against the
continuum cosine **cell average**, using exact semidiscrete evolution to remove
time error; (4) heterogeneous steady rates against total series resistance.
The four spatial levels are separately checked against the assembled operator.
Contact tests use R = 0, 1e-10, .3 and 1e5. Large R approximates insulation;
infinite R is not accepted. Local/global balance includes sources and boundary
heat rates. No physical experimental validation is claimed.

Limitations: dense small systems; fixed material properties and constant area;
no radiation, convection, nonlinear conduction or two-way mechanics feedback.
Extreme property/length/dt ratios can lose accuracy in floating point; positivity
of the mathematical backward Euler operator is not a universal numerical accuracy
guarantee. The large-resistance check allows 2e-9 relative flux error because
subtraction of nearby temperatures loses digits. Current studies do not establish
second-order convergence across material jumps or on arbitrary nonuniform grids.
Nonmatching transfer and the C++ mechanics interface remain pending.
