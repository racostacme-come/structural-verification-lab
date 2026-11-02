# Coupled refinement and contact/property verification

**1 November status:** The four-page [thermoelastic synthesis](../paper/paper.pdf)
and [source](../paper/paper.tex) now cover the coupled results. The text below
retains the milestone chronology; future-tense manuscript/refinement statements
describe that earlier date. Release remains subject to acceptance review.

This original synthetic study asks whether the heat/transfer/native-mechanics
pipeline reproduces the expected temporal and continuum spatial orders, and
whether thermal contact and material contrast preserve the independent
steady reaction-force identity. It extends the same thermoelastic research
question; it is not a new component repository.

## Reproduction and architecture

After the README installation and native build steps, run:

```sh
python -m thermomechanics.refinement_experiment
python -m pytest tests/test_coupled_refinement.py -q
```

The CLI accepts `--library PATH` and `--output PATH`; the default output is
`results/coupled-refinement.json`. That JSON records all study parameters,
45 contact/property cases, four levels for each refinement, and normalized
SHA-256 hashes of the original Python/native sources. Both CI platforms
reproduce it in `validation/coupled-refinement.json`. Only mechanics assembly
is C++17: heat operators, spectral reference and linear solves use NumPy.

`spectral_temperature(mesh, initial, time)` solves the insulated semidiscrete
system through the symmetric matrix `M^(-1/2) H M^(-1/2)`. This reference shares
the heat spatial operator but does not call backward Euler. Its unequal-mass
two-cell closed form is tested separately, and the uniform-grid results are
checked against the known cosine eigenvalue. `cosine_reference(x, amplitude)`
integrates the continuum strain analytically and calls neither assembly nor
transfer. Independent Gaussian quadrature checks its antiderivative.

## Separate error questions

On the unit interval, use homogeneous synthetic conductivity and volumetric
capacity of 1 in SI units, insulated ends, and initial temperature increment
`40 cos(pi*x)` K. Initial heat values are exact **cell averages**. The continuum
temperature is `40 exp(-pi^2*t) cos(pi*x)`. Mechanics has constant area
`2e-4 m^2`, fixed ends, and the existing two materials meeting at `x=.4 m`.
The mechanical mesh stays at 3+5 elements and does not match the heat grid.
Because the expansion coefficients differ, the zero-mean cosine produces a
nonzero reaction; a homogeneous expansion coefficient would conceal that error.

For zero end displacements and no distributed axial load,

```text
N = - integral(alpha*T dx) / integral(1/(E*A) dx),
u(x) = N*integral_0^x(1/(E*A) ds) + integral_0^x(alpha*T ds).
```

* **Temporal error:** fix 20 heat cells; run 10, 20, 40 and 80 backward-Euler
  coupled steps to `t=.1 s`. Compare with exact semidiscrete heat evolution
  and independent mechanics of that piecewise-constant cell field. This
  removes continuum heat error from the time-order measurement. Temperature,
  force and nodal-displacement orders are 0.97136, 0.98538 and 0.99261.
* **Spatial error:** use 10, 20, 40 and 80 heat cells with exact semidiscrete
  time evolution. Compare heat with continuum cell averages, and coupled
  force/displacement with the analytic continuum integral. Temperature orders
  are 1.98295, 1.99575 and 1.99894. Force and maximum nodal-displacement
  orders are approximately 2.00082, 2.00021 and 2.00005. The finest-grid
  errors are 0.00189043 K, 0.154849 N and 8.74181e-9 m respectively.
* **Transfer-only error:** the separate day-4 experiment fixes the source
  field and refines only the target mesh. Its integral conservation remains
  a distinct check, not evidence of accurate heat evolution.

The spatial study measures displacement at mechanical nodes. Between-node
P1 interpolation error on the fixed mechanics mesh remains; it is not
removed by heat refinement. The pre-existing affine-temperature study
separately measures mechanical interpolation error. Reported spatial orders
apply to this smooth zero-contact case, not to arbitrary discontinuous fields.

## Contact and material contrast

Use the original 7+11 thermal cells and 3+5 mechanical elements, boundary
increments 80 and 0 K, `k_left=13 W/(m K)`, and `E_left=70 GPa`. Independently
vary `k_right/k_left` and `E_right/E_left` over 0.1, 1 and 10, and contact
resistance over 0, 1e-10, .02, 2 and 1e5 m^2 K/W. Expansion coefficients
remain 23e-6 and 11e-6 per kelvin. The reference is the exact series resistance:

```text
q'' = 80 / (.4/k_left + R + .6/k_right),
temperature jump = q''*R,   heat rate = A*q''.
```

Each material's affine temperature is independently integrated in the bar
identity. All 45 cases have relative force error below 1.2e-14 and scaled
temperature error below 1.3e-15 locally. Maximum relative heat-rate error is
2.18e-8 at high resistance, where boundary rates subtract nearly equal
temperatures. The heat matrix condition number stays below 1853 for this
bounded parameter sweep; no claim extends to infinite resistance. CI allows
2e-7 relative flux error for `R=1e5`, 1e-10 otherwise, and 1e-11 scaled
temperature/force errors. These thresholds are distinct from reproduction
tolerances for tiny platform-dependent solver residuals.

Despite near-exact reactions, maximum nodal displacement error reaches
1.452e-6 m on these fixed nonmatching meshes. Piecewise-constant reconstruction
preserves complete material integrals, but not each partial integral ending
inside a heat cell. This retained error prevents a misleading conclusion that
exact global force guarantees exact local displacement. Mechanical bonding
remains perfect even with a thermal temperature jump.

## Scope and provenance

This is numerical verification of linear one-way coupling. It does not verify
transient continuum convergence with finite contact, validate a physical
specimen, establish realistic parameter identification, measure speedup or
include nonlinear/temperature-dependent properties. The existing discrete
energy checks support conservation but do not replace continuum accuracy.
The homogeneous transient and heterogeneous steady cases cover different
parts of the contract and are reported separately.

Topic provenance remains the read-only MMME heat/interface material recorded
in the research contract, plus eligible TUM finite-element concepts. All
parameters, implementations and experiments here are original synthetic work;
no lecture files, exercises, teaching code or datasets are redistributed.
Imported module snapshots and their manifests remain unchanged. The existing
two-page beam paper is still the released synthesis; the thermoelastic
manuscript is scheduled for the next contribution and gates the new release.
