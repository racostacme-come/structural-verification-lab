# Contributing

Use a short feature branch and a focused commit. Changes to numerical methods
need an independent analytical, conservation, or convergence check. Run the checks
listed in the README before merging with `git merge --no-ff`. Record failures and
their fixes honestly. Never add course handouts, solutions, exams or private data.

Python uses Ruff; C++ uses the committed clang-format configuration. GitHub Actions
builds on Linux and Windows. Do not infer a passing CI run from local test success.
