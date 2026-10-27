# Day 1: 27 October 2026

Original prescribed-temperature baseline, branch
`feat/thermoelastic-reference-baseline`; no major release or new-paper claim.

Executed locally with Python 3.12 and the existing integration environment:

* `python -m thermomechanics --output results/thermoelastic-baseline.json`:
  four reproducible boundary/material cases and four refinement levels.
* `python -m pytest tests/test_thermomechanics.py -q`: 25 passed.
* `python suite.py test`: 284 passed (250 component, 34 root including the
  25 new thermoelastic checks). Complete outputs: `validation/tests.json`.
* `python suite.py native --build-root ../../.native-structural`: four
  existing MSVC C++17 invariant CTests passed; `validation/native.json`.
  These are regression checks, not a new thermoelastic C++ implementation.
* `python suite.py verify`: 197 unchanged imported-file hashes pass. The
  Transport 192 and Vibration 209 hashes also passed in their repositories.
* `git diff --check`: passed before commit.

Uniform restrained heating gives N=-12880 N; free heating gives tip
displacement 0.00092 m and force zero to roundoff. The original two-material
restrained case gives N=-14448.636363636362 N and interface displacement
0.00013918181818181822 m. The largest scaled force error across four cases
is 1.42e-16 (rounded up). These establish implementation agreement for
supplied fields, not thermal-discretization or physical accuracy.

| Elements | Midpoint displacement error (m) | Observed order |
| --- | --- | --- |
| 4 | 1.4375e-5 | -- |
| 8 | 3.59375e-6 | 2.0 |
| 16 | 8.984375e-7 | 2.0 |
| 32 | 2.24609375e-7 | 2.0 |

The independent quadratic reference predicts these errors analytically;
nodal error stays below 9e-19 m. This separates interpolation error from
exact nodal values. Source hashes and raw measurements live in the JSON.

Hosted CI status is recorded in the external portfolio after the push.
Linux/Windows jobs run the new tests and archive fresh baseline JSON.
The existing manuscript job continues to build/check the released two-page
beam paper; the thermoelastic synthesis is scheduled for day 6.
Next: conservative transient heat/contact with independent flux, energy,
steady-state and resistance-limit checks.
