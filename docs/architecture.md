# Architecture and compatibility

The composition boundary is a Python process with NumPy arrays and explicit model parameters. Native C++ extensions keep their distinct package names. Every module builds from its own pyproject.toml, preserving build configuration and package metadata.

The root suite.py coordinates installation, source-integrity verification, isolated Python test processes, native CMake/CTest checks, the connected experiment, and manuscript compilation. Native checks compile four separate C++17 executables without Python bindings and reject an empty CTest suite. It does not monkey-patch module internals. integration.py contains all cross-module mappings so changes in sign conventions, units or constraints are reviewable.

## Source evaluation

All included modules had clean Git worktrees when imported. Each already supplied numerical tests and an academic paper. The consolidation retains the complete tracked contents of each source commit. The original CI records are historical evidence, while validation/tests.json records the new integrated run.

## Compatibility boundaries

ShearScope uses displacement/rotation Timoshenko elements. BuckleMap uses displacement/slope Euler-Bernoulli elements. Their stiffness matrices are not interchangeable. The bridge shares geometry and EI, subtracts analytical shear compliance, and compares a scalar cantilever identity only. FoldTrace and YieldLedger solve different nonlinear/constitutive models and are not silently combined with the elastic buckling model.

## Extension policy

Design a research question first. Change module APIs only with regression tests and a migration note. Keep imported licenses, commit provenance and component papers. Do not interpret low residuals or passing tests as evidence of physical validity.
