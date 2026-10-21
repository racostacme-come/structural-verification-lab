# Academic manuscript

[paper.pdf](paper.pdf) is the reading copy; [paper.tex](paper.tex) is the original
source. It covers the model, implementation, analytical verification, restraint
experiment and limitations. Tables derive from recorded synthetic results.

Install the project's `dev` extra and Tectonic or TeX Live with pdfLaTeX and
the standard LaTeX recommended/font packages. From the repository root:

```sh
python scripts/check_samples.py
python paper/build.py --engine tectonic
# Or: python paper/build.py --engine pdflatex
python scripts/check_manifest.py
pdftoppm -scale-to 1600 -png paper/paper.pdf paper/build/page
```

The builder regenerates `numbers.tex` from `examples/results/metrics.json`,
compiles into `paper/build`, rejects overfull boxes, undefined references,
empty pages and PDFs exceeding five pages, then delivers `paper.pdf`.
`validation.json` records page count, compiler and SHA-256 hashes of inputs and
PDF. Inspect every rendered page after changes; log checks cannot judge layout.
CI independently compiles with pdfLaTeX and uploads PDF, page PNGs and log.
Different engines can yield different PDF bytes and hash records.

To change experiments, run `bucklemap study --output examples/results`, review
the data, rebuild and visually inspect the manuscript. Commit source, data,
figures, generated tables, manifest and PDF together. The reproduction checker
compares values within documented tolerances without overwriting saved inputs.
