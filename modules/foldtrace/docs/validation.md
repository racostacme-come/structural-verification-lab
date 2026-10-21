# Numerical validation record

Local campaign: 2026-09-24, Windows, Python 3.14.5, MSVC 19.44.35219.0,
Release build. Python package versions are recorded in requirements-repro.txt.
Numbers below refer to the default r=0.2 campaign in results/validation.json.

| Quantity | Measured value |
| --- | ---: |
| Accepted states, including origin | 23 |
| Negative-stiffness states | 10 |
| Rejected default attempts | 0 |
| Final q (target is 2.2) | 2.247641 approximately |
| Maximum scaled equilibrium residual | 7.283453354323121e-13 |
| Maximum hyperplane residual | 1.717376241217039e-16 |
| Independent force discrepancy | 2.915983426943214e-14 |
| First analytical fold q | 0.42642777655779807 |
| Second analytical fold q | 1.573572223442202 |
| First fold load | 0.007547868508312872 |
| Maximum fold location discrepancy | 3.3306690738754696e-15 |
| Reconstruction orders | 2.02509091, 2.00272450, 2.00146466 |

The fold location comparison uses a closed-form expression in Python against
bisection of the C++ tangent, bracketed by consecutive continuation points.
The force comparison uses a separate direct strain-and-projection expression.
Neither benchmark extracts its expected value from the continuation algorithm.
Finite differences independently check tangent consistency and the energy-force
relationship at five geometries. Symmetry and stress-free states are also tested.

54 Python tests passed, with 94.71% Python statement coverage in the local run.
This coverage excludes compiled C++; it is not a whole-project coverage claim.
One native CTest executable passed in Release mode, exercising five invariants.
Tests cover full paths for ratios 0.01, 0.1, 0.2, 0.5, 1, 3, and 10; a difficult
scaled path exercises rejected attempts and step reduction. Invalid inputs,
iteration/step limits, partial-path preservation, exact-fold correction, output
schemas, CLI success/failure, and reproducibility are included.

Development checks initially reported import errors because tests were launched
before the asynchronous extension installation completed. After installation
finished, all 47 core tests passed without a solver change. The later report
added seven tests. Ruff found long lines in newly written reporting code; these
were reformatted before validation. Neither failure is presented as a passing run.

The independent interpolation benchmark halves maximum scaled path steps
0.16, 0.08, 0.04, 0.02. At each cell midpoint it compares linear interpolation of
the adjacent load values to the analytical curve. This tests useful path
reconstruction accuracy, not spatial or temporal convergence. The campaign
enforces scaled equilibrium residual <= 1e-10 and fold error <= 1e-9;
the default geometry additionally requires finest reconstruction order >= 1.8.

No wall-time speedup, general branch-tracking guarantee, dynamic response,
unconstrained stability, or material accuracy is inferred. The small model
permits strong analytical verification but is intentionally not a general FEM
solver. Results at untested geometries or extreme settings remain unverified.

The generated four-panel PNG was visually inspected for labels, overlap, and
agreement with the numerical data. The static unstable path is explicitly
distinguished from dynamic snap-through in the plot and README.

Local release checks also passed Ruff lint/format, clang-format, wheel and source
distribution builds, installation of the built wheel, all 54 tests after that
installation, the CLI campaign, the standalone example, and `pip check`.
The installed import path was verified to be in site-packages. Archive inspection
confirmed the expected C++ sources in the sdist, the compiled extension and Python
modules in the wheel, and no virtual environment, build tree, or Git directory.
The isolated MSBuild packaging step emitted MSB8029/MSB8012 infrastructure warnings
about temporary/output paths; compilation, linking, installation, and wheel
execution succeeded. C++ compiler warnings are treated as errors.

The first remote run reported deprecated Node 20 action runtimes. The workflow
was updated to commit-pinned checkout 7.0.1, setup-python 7.0.0, and upload-artifact
7.0.1, using tag commits verified through GitHub's API. The numerical code was
unchanged; the updated workflow is checked before merging the reporting branch.
