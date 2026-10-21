# Validation record

Initial local run: 2026-09-20, Windows x64, Python 3.14.5, MSVC 19.44.35219,
CMake 4.4.3. Exact Python dependency versions are in `requirements-repro.txt`.
These are numerical verification results, not experimental certification.

## Observed results

| Quantity | Result |
| --- | ---: |
| Python tests | 25 passed |
| Python statement coverage | 97% (does not measure C++ coverage) |
| Native CTest | 1/1 passed in Release |
| Static campaign | 80 cases |
| Modal campaign | 10 eigensolves, 3 modes each |
| Maximum static free-DOF relative residual | 8.09e-8 |
| L/h = 100, 8 elements, full integration, computed/exact tip | 0.0195771 |
| L/h = 100, 8 elements, reduced integration, computed/exact tip | 0.9960941 |

At L/h = 200 with 80 reduced-integration elements:

| Mode | Computed Hz | Euler–Bernoulli limit Hz | Relative difference |
| --- | ---: | ---: | ---: |
| 1 | 4.112605 | 4.112609 | -0.000089% |
| 2 | 25.779595 | 25.773280 | +0.024502% |
| 3 | 72.221447 | 72.165929 | +0.076931% |

The exceptionally small first-mode difference can include cancellation between
discretization error and the physical shear/rotary-inertia correction. It does
not imply that the method is this accurate for other beams. The full-integration
first-mode frequency at the same mesh is 7.126830 Hz, about 73.29% above the EB limit.

The reduced-integration tip error decreases by a factor of four on each uniform
mesh doubling (4, 8, 16, 32 elements), checked to 0.2% relative tolerance on the
ratio. Tip-force, tip-moment and uniform-load continuum comparisons at 64 elements
are checked to 0.02% relative tolerance. Reactions are checked to 1e-7 absolute
tolerance in the test's SI units. Modal frequencies are checked within 0.2% of
the slender limit, mass orthogonality within 1e-10, and eigenpair residuals within
1e-5 relative tolerance.

## Checks executed locally

- Built and installed the C++ extension through pip.
- Ruff lint and formatting checks; clang-format dry-run checks.
- Pytest including a subprocess CLI run and complete campaign generation.
- Native CMake Release build with compiler warnings treated as errors; CTest.
- Source distribution and wheel built via `python -m build` (wheel from sdist).
- CLI campaign and standalone example executed; `pip check` passed.
- Generated figure inspected visually; numerical CSV/JSON retained in `results/`.

MSBuild emitted temporary-directory and path-spelling warnings during isolated
packaging (MSB8029/MSB8012). Compilation, linking and installation succeeded.
They are build-system warnings, not numerical test failures. No compiler source
warnings were reported. Windows is the initial local test platform; the GitHub
Actions run records separate Linux/Windows and Python-version results.

## Reproduce

Follow the README environment instructions, then run its check commands and
`shearscope campaign --output out`. Compare CSV values with relative tolerance
1e-5 and absolute tolerance 1e-7; do not require bitwise equality of plots or
linear-algebra results. Keep `results/` as the original record unless intentionally
regenerating it and documenting a change.

Untested limits include nonlinear response, damping, 3D/warping effects, thick-beam
modal continuum accuracy, and extreme conditioning beyond the supplied campaign.
