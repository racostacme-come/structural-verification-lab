# Day 4: conservative nonmatching coupling, 30 October 2026

Branch: `feat/conservative-nonmatching-coupling`, continuing the original
27 October thermoelastic cycle. No merge or release is claimed here.

Implemented overlap integration with explicit domain/property contracts,
one-way heat/transfer/mechanics state advancement, and a reproducible native
CLI with four target meshes and 60 transient steps across three contact
resistances. All synthetic parameters and source hashes are versioned in
`results/coupled-transfer.json`. The existing C++17 assembler is exercised
end to end; no new native kernel or performance improvement is claimed.

Actually executed locally (Windows, existing integration virtual environment):

* `python -m thermomechanics.coupled_experiment`: four transfer cases and
  three 20-step histories. Maximum scaled integral error 4.0536e-16;
  force/nodal errors across coupled histories 1.7463e-15/1.6065e-15;
  energy residual 2.1074e-12. Native/Python force differences zero locally.
* `python -m pytest tests/test_transfer.py -q`: 25 passed after correcting
  an overly strict convexity assertion. Its initial exact inequality failed
  by one floating-point ULP; the test now allows 8 eps times the field scale.
  No numerical implementation or experiment data were changed to hide it.
* `python suite.py test`: 352 passed (250 component, 102 root, including
  those 25 new checks). Full process outputs are in `validation/tests.json`.
* `python suite.py native --build-root ../../.native-structural`: all five
  MSVC Release CTests passed; commands/output in `validation/native.json`.
* All three suites' `suite.py verify`: 197+192+209=598 original hashes pass.
  All 19 original component working trees inspected clean. No imported
  snapshots, integrity manifest, original history or course source was edited.
* Baseline/native experiment source manifests explicitly regenerated to
  include the three new Python files; previous numerical results unchanged.
* `git diff --check` passed. Previous native milestone's hosted run
  40963142071 verified successful for Linux, Windows and manuscript.

The current branch adds the coupled CLI to Linux/Windows CI. Its hosted
status is recorded in the parent portfolio after publication; local checks
above do not imply hosted success. Existing two-page beam synthesis unchanged;
the extension paper remains scheduled for 1 November, with compilation,
<=5-page validation and full visual inspection before any release.

These checks establish source-field transfer accuracy, mechanics consistency
and discrete heat balance, not coupled continuum thermal accuracy or physical
validation. Next: independent coupled temporal/spatial refinement and contact
studies on 31 October; preserve all outcomes without tuning a convergence story.
