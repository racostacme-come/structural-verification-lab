# Changelog

## v0.2.0 - Thermoelastic transfer

Original one-way thermal-to-mechanical verification on a restrained two-material
bar, integrated with the existing six-component laboratory:

- Conservative cell-centred transient conduction with finite thermal contact
  resistance and independent energy, flux and temporal/spatial checks.
- Original C++17 stiffness and thermal-load assembly through C ABI v1 and an
  explicit Python interface, with standalone native invariants.
- Material-aware overlap transfer between nonmatching meshes, a coupled CLI,
  and independent continuum force/displacement references.
- Four-level refinement and transfer experiments, 60 coupled steps and 45
  steady contact/property cases. Accurate reaction force coexists with local
  displacement reconstruction error; naive interpolation results are retained.
- Original four-page synthesis, generated measurement tables, 28-file manuscript
  integrity checks and independent Linux/Windows/pdfLaTeX CI.

See [acceptance evidence](validation/release-acceptance.md) for measured values,
verification boundaries and publication gates. This is synthetic numerical
verification, without a physical-validation or speedup claim. Native code
assembles mechanics; heat, transfer, constraints and dense solves use Python.

## v0.1.0 - Beam verification consolidation

Six versioned original components, shared beam compliance/buckling experiment,
independent analytical checks and reviewed beam synthesis. The original
[paper](https://github.com/racostacme-come/structural-verification-lab/blob/v0.1.0/paper/paper.pdf)
remains available. All component histories, notices and pinned source hashes
are preserved.
