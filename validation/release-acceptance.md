# Thermoelastic transfer: v0.2.0 acceptance review

Reviewed 2 November 2026 (Europe/Paris), calendar day 7 of the cycle begun
27 October. This is an original extension of the existing laboratory. No work
is attributed to the missed 28 October milestone; heat and native assembly
were separate contributions on 29 October.

## Decision and evidence boundary

The nine research-contract criteria below pass for the implemented scope.
The release candidate changes documentation and records fresh local checks;
numerical sources, measured JSON, manuscript and imported snapshots are
unchanged from `80da7beccf8b8c936a70877597df7eb287f51960`.
Merge/publication additionally requires successful Linux, Windows and paper
jobs on the final candidate, and equality of the tested and merged Git trees.
The existence of this document alone is not evidence of publication: the
GitHub v0.2.0 tag/release identifies the published commit.

The preceding [hosted run](https://github.com/racostacme-come/structural-verification-lab/actions/runs/41355573259)
passed all three jobs. Both numerical artifacts reproduce 381 Python tests,
five native CTests, four-level refinement/transfer studies and 45 contact
cases. All three downloaded artifact ZIP digests were verified on 1 November.
The independent four-page pdfLaTeX PDF passed all 28 manuscript hashes and
complete visual review. The final candidate must repeat these hosted gates.

## Acceptance matrix

| Contract criterion | Evidence and decision |
| --- | --- |
| 1. Free expansion and uniform restraint, scaled error <=1e-10 | Python and C++ analytical patches in `test_thermomechanics.py`, `test_native_thermomechanics.py` and `test_transfer.py`; four-case baseline maximum force error 1.42e-16 and nodal error 1.18e-16. Pass. |
| 2. Python/C++ stiffness and load <=1e-11; independent native tests | 36 seeded native cases, maximum force/nodal scaled error <1.1e-14; native stiffness/load differences zero locally. Five standalone native CTests pass, including original thermoelastic invariants. Pass. |
| 3. Weighted overlap integral <=5e-12; honest naive comparison | Four target grids, maximum scaled integral error 4.06e-16. Naive interpolation matches on two coarse grids, then reaches 0.5267% force error; this nonmonotone result is retained. Pass. |
| 4. Independent supplied-field force identity | Affine-layer oracle, randomized layered fields and independent source-cell integrals; coupled force/nodal scaled error <1.8e-15 over 60 steps. Oracle checks do not certify heat-field accuracy. Pass. |
| 5. Separate four-level studies | Backward Euler orders 0.97136, 0.98538, 0.99261 against the semidiscrete reference; smooth thermal spatial orders 1.98295, 1.99575, 1.99894; force/nodal spatial orders approach 2. Transfer-only studies hold the source field fixed. Pass. |
| 6. Energy, equilibrium, contact and bounded limits | Coupled maximum scaled energy residual 2.11e-12. Forty-five independent steady contact/property cases include R=0 and R=1e5 m2 K/W; maximum local relative heat-rate error 2.18e-8 versus 2e-7 high-R gate, force error <1.2e-14, matrix condition number <2e4. Pass within these bounds. |
| 7. Cross-platform integration and provenance | Previous exact numerical/manuscript revision passed Linux and Windows; fresh local 381 Python checks and five MSVC CTests pass. All 598 imported file hashes across three suites pass; 19 original repositories are clean. Final candidate CI and merge-tree identity remain publication gates. |
| 8. Original <=5-page synthesis | Four pages including references; source, PDF and generated tables bind to measured JSON through a 28-file manifest. Tectonic and independent CI pdfLaTeX builds passed; all pages reviewed. Fresh manifest verification and rendering on 2 November pass. |
| 9. Release scope and truthful claims | Synthetic numerical verification only; no physical validation, performance claim, plasticity, two-way feedback or transient-contact continuum convergence claim. Release only after final external gates. |

Machine-readable measurements, input SHA-256 values and fresh rendering hashes
are in [release-acceptance.json](release-acceptance.json). Numerical thresholds
are enforced in the existing numerical tests, independently of this report.

## Checks executed on 2 November

```sh
python suite.py native --build-root <short-build-directory>
python suite.py test
python manuscript.py
pdftoppm -scale-to 1400 -png paper/paper.pdf paper/renders/release
```

Fresh local results: 250 component Python checks plus 131 root checks; all five
MSVC Release native CTests passed. Full tests reran the independent references
and versioned numerical reproduction. All 28 manuscript hashes and four pages
passed. All four freshly rendered pages were visually inspected: equations,
tables, commands and references are complete, legible and without clipping.
The unchanged paper was already compiled with two independent engines on
1 November; the final hosted candidate rebuilds it again with pdfLaTeX.

Private coordination audit: 34 indexed lecture hashes unchanged, no course
files modified since the previous automation run. Course sources remained
read-only. No private lecture, exercise, solution, dataset or teaching code
is included in the release.

## Interpretation and remaining limitations

Implementation agreement, discrete temporal error, continuum error and physical
validity are distinct. Global reaction agreement alone is insufficient:
steady contact cases retain nodal displacement error up to 1.4513e-6 m.
The heterogeneous transient-contact demonstration is a discrete model output;
its continuum convergence is not established by the separate smooth-field and
steady studies. Large-resistance tests cover a finite, conditioned range.
Dense Python solves and overlap matrices target small verification problems.

The original v0.1.0 beam synthesis and all six Structural component papers
remain available with their original provenance. No imported code or integrity
manifest changed during this cycle. Further work should address separately
justified research questions rather than weaken these acceptance thresholds.
