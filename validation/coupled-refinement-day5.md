# Day 5: independent coupled validation, 31 October 2026

Branch: `feat/coupled-refinement-validation`, continuing the original
27 October-2 November cycle. No merge or release in this contribution.

Added an exact semidiscrete spectral heat reference, analytic cosine-field
continuum mechanics, four-level native coupled temporal/spatial refinements,
and 45 heterogeneous steady contact/property cases. The numerical CLI and
complete measured values are `thermomechanics/refinement_experiment.py` and
`results/coupled-refinement.json`; methods/limitations are documented in
`docs/coupled-refinement.md`.

Executed locally on Windows with the existing integration environment:

* `python -m thermomechanics.refinement_experiment`: temperature/force/nodal
  temporal orders 0.97136, 0.98538, 0.99261; temperature spatial orders
  1.98295, 1.99575, 1.99894; force/nodal orders approximately 2.00082,
  2.00021, 2.00005. Spectral versus known discrete cosine residual <4.8e-12 K.
* 45 contact/property cases: maximum relative force error 1.1904e-14,
  scaled temperature error 1.2879e-15, relative heat-rate error 2.1702e-8,
  and heat condition number 1852.1. High resistance exposes cancellation
  in small boundary heat rates. Maximum local nodal error 1.4513e-6 m
  is retained and explained; exact global reactions do not imply exact
  partial strain integrals on nonmatching meshes.
* `python -m pytest tests/test_coupled_refinement.py -q`: 16 passed, including
  Gaussian-quadrature checks of the continuum oracle, unequal-mass two-cell
  spectral reference, separate orders, contact limits and full reproduction.
* `python suite.py test`: 368 passed (250 component, 118 root). Actual
  per-suite outputs are in `validation/tests.json`.
* `python suite.py native --build-root ../../.native-structural`: all five
  MSVC Release native tests passed; commands/output in `validation/native.json`.
* All three `suite.py verify` checks: 197+192+209=598 imported hashes pass.
  All 19 original repositories inspected clean; 34 indexed lecture hashes
  unchanged, and no course file modified since the prior scan. Sources stayed
  read-only. Imported snapshots and their integrity manifests were not changed.
* Baseline/native/coupled result source manifests explicitly regenerated to
  include the new study module; all previous numerical values checked exactly
  unchanged before writing. No paper artifact changed.

The workflow now reproduces the new CLI on Linux and Windows, alongside the
existing numerical/native/manuscript gates. Hosted status is recorded in the
parent portfolio after push; local success is not a hosted success claim.
The unchanged two-page paper is the released beam synthesis. Next: use these
versioned measurements for the <=5-page extension manuscript on 1 November,
compile/check/render every page, then release only if all acceptance gates pass.

These are synthetic numerical checks, not physical validation. Spatial orders
apply to a smooth zero-contact transient. Finite-contact validation here is
steady; no transient-contact continuum convergence claim is made. Mechanical
interpolation between the fixed target nodes is not eliminated by heat refinement.
