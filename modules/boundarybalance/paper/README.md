# Rebuild and inspect the manuscript

The paper is generated from original measurements in `examples/results`.
It is a numerical study, with no fabricated authorship or affiliation.

```bash
python -m pip install -e ".[dev]"
boundarybalance study --output examples/results
python scripts/check_samples.py
python paper/build.py --engine pdflatex
python scripts/check_manifest.py
pdftoppm -scale-to 1600 -png paper/paper.pdf paper/build/page
```

Install TeX Live packages `texlive-latex-base`, `texlive-latex-recommended`, and
`texlive-fonts-recommended`, or use Tectonic via `--engine /path/to/tectonic`.
Poppler supplies `pdftoppm`. The build script runs pdfLaTeX twice or Tectonic's
automatic passes, generates `numbers.tex` directly from JSON, rejects undefined
references and overfull boxes, enforces 1–5 pages, and records input/PDF hashes
in `validation.json`. Inspect every rendered page after changes; text extraction
and a successful compiler are not a visual layout check.

The independent CI paper job first reproduces numbers and verifies the committed
manifest, then compiles using pdfLaTeX. Its PDF, log, manifest and rendered pages
are retained in the `verified-paper` Actions artifact. Font and metadata differences
can change PDF hashes between engines. Compare the recorded input hashes to tie
each artifact to the same numerical study.

To intentionally revise results, regenerate the study and paper together. Do not
edit the generated table macros to invent or tune numerical outcomes. The PDF
includes all references in the enforced page count.
