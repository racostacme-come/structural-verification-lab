# Contributing

Use a focused branch from main and the configured Git author. Include a numerical
regression when behavior changes. Do not add course materials, private datasets,
environments, or credentials. Explain changes to regenerated sample metrics.

Install `.[dev]`. Run Ruff lint/format, pytest, the benchmark CLI, the example, and
`python -m build` before merging. Test the installed wheel too. Published branches
should pass the full CI matrix before merging. Do not rewrite shared history.

Preserve the documented residual weighting or document a deliberate change. Do
not loosen analytical tolerances to conceal failures. For timing claims distinguish
one solve from the complete adaptive history; record hardware and methodology.
