# Structural Verification Lab

Beam discretization, stability, nonlinear equilibrium, plasticity, adaptive meshes and boundary enforcement.

This major project integrates 6 previously independent research implementations into one reproducible laboratory. Full Python and C++ source, tests, datasets, documentation, licenses and component papers are included under `modules/`; no submodule downloads are required.

## Start here

An original [thermoelastic-transfer study](docs/thermoelastic-transfer.md) is
in progress (27 October-2 November target). Day 1 provides an independent
restrained-bar continuum oracle, Python FE baseline, analytical patches and
four-level displacement interpolation study. Run `python -m thermomechanics`
after setup; [versioned baseline](results/thermoelastic-baseline.json).
Conservative heat/contact and the [C++17 assembly backend](docs/native-thermoelastic.md)
are implemented. [Conservative nonmatching transfer and the coupled CLI](docs/coupled-transfer.md)
now connect heat, transfer and C++ assembly; independent coupled refinements are next.
The paper below remains the released v0.1.0 beam synthesis until the new
study has measured results and completes manuscript validation.

Python 3.11+ (CI uses 3.12), Git, internet for initial dependency installation, and a C++17 compiler are required. On Windows use Visual Studio 2022 with Desktop development with C++ and a Windows SDK. On Linux use GCC or Clang. Use a fresh virtual environment.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python suite.py verify
python suite.py install
python suite.py native
python suite.py test
python suite.py experiment
python -m thermomechanics.native_experiment
python -m thermomechanics.coupled_experiment
```

The runner stops on install failure and reports every module's test outcome. Tests execute in separate processes to avoid duplicate test-file/import names. `validation/tests.json` contains the actual outputs; `results/integration.json` records the connected experiment. Rerunning overwrites only these generated reports. For Windows builds, if compiler discovery fails, use a Visual Studio developer shell or set `CMAKE_GENERATOR=Visual Studio 17 2022`.

Deep Windows checkouts can exceed MSBuild's path limit during compiler detection. Pass `python suite.py native --build-root <short-build-directory>` to keep generated paths short; choose an empty build directory for this suite.

## Included components

| Component | Implementation and original paper |
| --- | --- |
| shearscope | [Source and usage](modules/shearscope/README.md), [paper](modules/shearscope/paper/paper.pdf) |
| bucklemap | [Source and usage](modules/bucklemap/README.md), [paper](modules/bucklemap/paper/paper.pdf) |
| foldtrace | [Source and usage](modules/foldtrace/README.md), [paper](modules/foldtrace/paper/paper.pdf) |
| yieldledger | [Source and usage](modules/yieldledger/README.md), [paper](modules/yieldledger/paper/paper.pdf) |
| meshbeacon | [Source and usage](modules/meshbeacon/README.md), [paper](modules/meshbeacon/paper/paper.pdf) |
| boundarybalance | [Source and usage](modules/boundarybalance/README.md), [paper](modules/boundarybalance/paper/paper.pdf) |

## What is integrated

A common cantilever geometry and SI material definition connects ShearScope compliance and BuckleMap stability. A dimensionless identity tests the handoff of length, bending rigidity and boundary conditions. Full versus reduced quadrature and mesh refinement expose different discretization errors. The remaining nonlinear, plasticity, mesh and boundary modules provide complementary verification methods through the same build/test interface.

This is a modular integration, not a universal multiphysics solver. Component APIs retain their explicit assumptions and units. Not every module is physically coupled to every other module.

## Paper and evidence

Read the new [synthesis paper](paper/paper.pdf) and [LaTeX source](paper/paper.tex).
The component papers are retained as detailed supporting studies. New numerical claims are based on the connected experiment, not on simply adding together earlier performance claims.

Compile with Tectonic or pdfLaTeX:
```sh
python suite.py paper --engine tectonic
```
The command rejects unresolved references, overfull boxes and more than five pages. CI independently builds the paper. [Architecture](docs/architecture.md), [source provenance](suite.json), and [development roadmap](ROADMAP.md) describe scope and future work.

## Provenance and development

Imported snapshots are pinned to original Git commits in `suite.json`, with SHA-256 for every imported file. `verify` detects unrecorded changes. Original repositories and histories remain available; this repository starts a new integration history. Intentional module changes require a documented patch record and manifest refresh. [Provenance repair](docs/provenance.md) explains canonical source bytes and the explicit MeshBeacon attributes patch. MIT notices are retained for every component; the integration layer is also MIT licensed.

CI runs all component Python tests, integration checks, and five native C++17 invariant tests on Linux and Windows. Component Python extensions are built during install; the original thermoelastic library is built by `native` before tests. The `native` command writes commands, return codes and compiler/CTest output to `validation/native.json`. Independent integration tests check modulus scaling, geometrically similar beams and refinement toward analytical continuum limits. Sanitizers and full historical reproduction campaigns are outside this gate. Recorded results are synthetic numerical verification, not experimental physical validation.
# Thermal contact extension

The original [conservative heat/contact solver](docs/thermal-contact.md) adds
backward Euler conduction, independent energy/flux checks and separate temporal
and spatial refinement studies. Reproduce with
`python -m thermomechanics.heat_experiment`. This is an in-progress extension;
the existing synthesis paper describes the released beam study.
