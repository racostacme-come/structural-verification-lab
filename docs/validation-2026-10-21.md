# Consolidation validation: 21 October 2026

## Local evidence

- Windows, Python 3.12.2; 250 component Python tests and nine root integration/integrity tests passed. Complete output is in `validation/tests.json`.
- Four C++17 invariant executables built with MSVC 19.44 and passed CTest. The `native` command records configure/build/test commands and complete output in `validation/native.json`.
- 197 imported files were compared with original pinned Git blobs: 196 exact matches and one explicitly documented attributes patch. See `validation/source-audit.json` and `docs/provenance.md`.
- Fresh clones with `core.autocrlf=false` and `core.autocrlf=true` both passed all 197 snapshot hashes; `validation/checkout-audit.json` records the tested code commit. All 19 original component repositories remained clean.
- The connected cantilever experiment was reproduced. At 32 elements, relative buckling error is 7.919e-9, reduced compliance error 2.441e-4 and full-quadrature compliance error 0.7579. The scalar product is 0.822266 versus pi^2/12 = 0.822467.
- Independent tests check elastic-modulus scaling at factors 0.1, 3 and 10; geometric scaling at factors 0.25 and 4; and continuum refinement at slenderness ratios 20 and 100. These are separate from the byte-integrity checks and from experimental physical validity, which is not claimed.
- Tectonic compiled the revised synthesis to two pages. Both rendered pages were visually inspected: equations, table, references and page boundaries are legible and unclipped. No undefined references or overfull boxes. A nonfatal underfull line in the bibliography does not clip content.

## Failures found and repaired

The earlier Linux workflow failed on MeshBeacon's platform-dependent archive hashes. Canonical source hashes plus one documented LF-checkout attributes patch repair the cause while retaining provenance. No original repository was changed.

The initial local native build failed during compiler detection because the deep checkout exceeded an MSBuild path limit. The runner now accepts `--build-root`; rerunning in a short directory under the permitted project root passed all four tests. A first manuscript invocation used a relative executable path that Windows could not locate; rerunning with the absolute Tectonic executable succeeded.

## External gate

The current feature branch must pass Linux, Windows and independent pdfLaTeX CI before merging. The workflow uploads native/numerical reports and the compiled manuscript for inspection. Final workflow links and publication state are recorded in the parent portfolio and the release record after completion.
