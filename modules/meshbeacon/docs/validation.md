# Numerical validation record

Generated locally on 2026-09-24 using Windows and Python 3.14.5. Machine-readable
metrics are in `results/validation.json`, with full history in `convergence.csv`.
These are measurements, not certified bounds.

Local checks passed: 46 pytest cases, 97.75% Python statement coverage, Ruff lint
and formatting, wheel and source builds, installed-wheel tests, CLI, example,
and `pip check`. Distribution contents were inspected, and the generated figure
was visually checked. The wheel import path was verified to be in site-packages.

| Check | Measured value |
| --- | ---: |
| Adaptive elements / free DOFs | 426 / 425 |
| Solves including initial mesh | 22 |
| Adaptive energy error | 0.02550963935 |
| Uniform energy error, 426 elements | 0.1642009871 |
| Adaptive / uniform error | 0.1553561876 |
| Residual estimator | 0.1250757916 |
| Effectivity | 4.903079573 |
| Sine L2 order, 128 to 256 elements | 1.999985777 |
| Sine energy order, 128 to 256 elements | 0.999989137 |
| Max nodal change, quadrature 8 to 16 | 7.03e-9 |
| Energy norm change, quadrature 16 to 24 | 2.64e-16 |

The stopping threshold is an estimator value of 0.15, not an energy-error promise.
Equal cell counts imply equal free DOFs because both endpoints are prescribed.
The exact bump solution is never consulted by the marking or refinement logic.

For `u=x*(1-x)` on a cell of length h, squared interpolation errors are `h^3/3`
in energy and `h^5/30` in L2. Tests use these closed forms on a nonuniform mesh.
For two half-length cells and `f=2`, both squared residual indicators are 0.75;
this independent hand calculation checks the estimator weights.

During development an exact-equality test for a new decimal midpoint failed by
1.11e-16. The test now tolerates floating-point rounding for the new point while
requiring exact preservation of old nodes. The solver was not changed to hide it.

See Actions for actual platform/build results. Statement coverage does not establish
mathematical correctness. Scope and quadrature limitations are in the README.
