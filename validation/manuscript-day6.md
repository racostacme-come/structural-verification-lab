# Day 6: thermoelastic synthesis and evidence gates

1 November 2026, Europe/Paris; same cycle started 27 October, target 2 November.

## Contribution

Replaced the current synthesis with an original four-page thermoelastic paper.
The released v0.1.0 beam paper remains in Git history and at its release tag.
Three numerical tables and summary macros are generated directly from the
versioned coupled-transfer/refinement JSON. The manuscript explains governing
models, overlap/native interfaces, independent references, measured convergence,
45 steady contact/property cases and the retained local displacement error.
It explicitly excludes physical validation and transient-contact continuum
convergence. Public references were checked against author/university records.

`manuscript.py` rejects incomplete/stale numerical source manifests and stale
tables. The build report binds 28 files including all 17 Python/native sources,
JSON, tables, manuscript, generator, runner, provenance manifest and PDF. The
checker accepts an external CI PDF/report without overwriting local evidence.
The paper CI job now verifies the report immediately after pdfLaTeX compilation.
Thirteen acceptance/mutation checks cover altered data, tables, source, PDF,
incomplete manifests and page-count violations. Hashes prove identity only.

## Checks actually executed

- 381 Python tests passed: 250 component and 131 root, including 13 manuscript checks.
- All five MSVC Release C++17 CTests passed after native rebuild orchestration.
- Numerical reproduction tests reran coupled transfer/refinement; versioned
  numerical source hashes and result files remain unchanged.
- All 598 imported snapshot hashes across the three suites passed; all 19 original
  component repositories remain clean. No imported files/manifests were changed.
- All 34 indexed lecture hashes unchanged, with no newly modified source files
  since the previous run. Course folders remained read-only.
- Local Tectonic compiled four pages including references, without unresolved
  references/citations, overfull boxes or LaTeX errors. All four Poppler-rendered
  pages visually reviewed; PDF/render hashes in manuscript-visual-review.json.
- Initial drafting exceeded the page limit; subsequent layout iteration exposed
  an overfull command block. Both were corrected before the successful build.
  An initial Python raw-string syntax error was also corrected. No failed draft
  is represented as passing validation.

## Publication and next step

This local record precedes the feature push. Hosted Linux/Windows numerical
jobs and independent pdfLaTeX manuscript are reviewed after push and recorded
in the persistent portfolio/ROADMAP. No merge or release is claimed here.
On 2 November perform the acceptance review, merge/release only if all gates
pass, then select the next coherent Transport question. Days 8-9 remain for
justified repairs; retain this project if a gate is unresolved.
