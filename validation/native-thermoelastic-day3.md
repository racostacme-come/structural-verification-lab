# Native thermoelastic milestone — 29 October 2026

Calendar day 3 of the cycle begun 27 October. This is a second contribution
on 29 October, after the early heat/contact contribution. Nothing is
backdated to 28 October. Branch: `feat/native-thermoelastic-assembly`.

## Changes and scope

- Original C++17 tridiagonal stiffness and thermal-load assembly, C ABI v1,
  standalone runtime invariants, and owned-buffer Python ctypes interface.
- Optional assembly callable in the existing mechanics solve; boundary
  elimination and dense linear solve remain NumPy operations.
- Seeded 36-case layered-bar experiment with explicit source hashes;
  independent continuum reaction/displacement comparisons.
- Required native/Python checks in Linux and Windows CI. Builds precede
  tests; missing native libraries fail instead of silently falling back.
- README, API/build documentation and roadmap updated. All imported sources
  and original histories preserved. No private course material added.

## Checks actually executed locally

Windows, MSVC 19.44.35219.0, Release C++17:

1. `python suite.py native --build-root <short-local-build>`: five CTests
   passed (four component regressions plus one original thermoelastic test).
2. `python -m thermomechanics.native_experiment`: 36 synthetic cases, seed
   20261029; saved `results/native-thermoelastic.json`.
3. Regenerated `results/thermoelastic-baseline.json`; all numerical data
   exactly equals the previous version, with only source hashes updated.
4. Targeted native and baseline Python tests: 38 passed.
5. `python suite.py test`: 327 passed (250 component, 77 root including
   13 new native integration checks). Raw outputs in `tests.json`.
6. Integrity verification for all three suites: 197/192/209 hashes passed.
   All 19 original component repositories were inspected and clean.
7. `git diff --check` passed. No changes to `modules/`, `suite.json` or
   `paper/`. No numerical failures occurred in this run.

Maximum scaled errors from the versioned experiment:

| Metric | Maximum |
| --- | ---: |
| C++ versus Python stiffness | 0 |
| C++ versus Python thermal load | 0 |
| Element force versus continuum | 1.0343e-14 |
| Nodal displacement versus continuum | 1.0786e-14 |
| Free-node equilibrium residual | 7.0439e-15 |
| Support force versus continuum | 7.0439e-15 |

Each gate is <=1e-11 on the documented physical scale. Roundoff-level
agreement is not a thermal-convergence, mesh-transfer, performance or
physical-validation result. The energy-gradient check and analytical native
checks are independent of the Python/C++ matrix comparison.

The earlier heat commit `25e6738` passed Linux, Windows and manuscript jobs
in hosted run 40788796353, inspected today. The new feature's hosted checks
start after push; its publication status is maintained in the portfolio.
No merge or release is attempted at this implementation milestone.

## Manuscript and next step

The existing two-page synthesis remains the released beam paper, already
compiled and reviewed in the consolidation cycle. It is unchanged today.
The extension paper is not complete: first implement nonmatching overlap
transfer and coupled CLI experiments, then independent coupled refinements,
and draft/rebuild/render the <=5-page extension synthesis. Release remains
conditional on those gates and cross-platform checks.
