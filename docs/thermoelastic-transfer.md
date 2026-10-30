# Thermoelastic transfer: research contract

Status 30 October, calendar day 4: baseline, conservative heat/contact,
[native C++17 assembly](native-thermoelastic.md) and
[nonmatching transfer with coupled CLI](coupled-transfer.md) are implemented.
Independent coupled refinements are next. Two contributions occurred on
29 October; no contribution is backdated to 28 October.
The day-1 descriptions below retain their original scope.

Cycle started 27 October 2026 (Europe/Paris). Target release: 2 November;
earliest 31 October; 3-4 November only for justified repairs. This is an
original extension of Structural Verification Lab, developed on
`feat/thermoelastic-reference-baseline`. The v0.1.0 paper remains the reviewed
beam-consolidation paper; it does not yet describe or validate this extension.

## Research question

How do thermal contact resistance and transfer between nonmatching thermal
and mechanical meshes affect transient reaction force and displacement in
a restrained two-material bar?

The model is one-dimensional, small-strain, quasi-static linear elasticity,
perfect mechanical bonding, piecewise constant properties, and prescribed
or subsequently computed temperatures. Thermal contact can be imperfect
without allowing a mechanical displacement jump. All examples are original
synthetic data, not identified properties or experimental validation. No
plasticity, inertia, thermally induced property changes, two-way coupling,
lateral constraints, or contact opening are included.

## Equations, units, and signs

Coordinates and displacement use m, area m^2, Young's modulus Pa, force N,
temperature increments K, thermal expansion 1/K. Expansion `alpha` is not
thermal diffusivity. Positive force means tension; positive displacement
points right. With no distributed axial loads,

\[
N'=0,\qquad N=EA(u'-\alpha\Delta T),\qquad
N=\frac{u_R-u_L-\int_0^L\alpha\Delta T\,dx}
        {\int_0^L(EA)^{-1}\,dx}.
\]

The support forces on the bar are `(-N, N)`. A traction-free right end has
`N=0`, and the left displacement fixes rigid translation. Displacement at
any coordinate follows by integrating `N/(EA)+alpha*deltaT` from the left.
The independent oracle integrates affine temperatures analytically in each
layer; displacement is generally quadratic even when nodal FE values are exact.

For an element of length h, constant EA and exact mean temperature increment,

\[
K_e=\frac{EA}{h}\begin{bmatrix}1&-1\\-1&1\end{bmatrix},\qquad
f_e^{th}=EA\alpha\overline{\Delta T}\begin{bmatrix}-1\\1\end{bmatrix}.
\]

Eliminate prescribed displacements including the `K_fc u_c` term. Compute
support reactions as `K u-f`, never by negating an internal-force convention
without checking the end orientation. Both constrained nodes and the
traction-free right end are supported, including a single-element mesh.

## Architecture and API

* `thermomechanics/model.py`: validated immutable Material and BarMesh SI
  contracts. Every material and area jump must coincide with a mechanical
  node. Temperature input contains **element means**, not nodal values.
* `reference.py`: Layer and analytical compatibility oracle, independent of
  FE assembly and mesh construction; affine temperature within each layer.
* `mechanics.py`: original dense NumPy P1 assembly and boundary elimination;
  independent oracle tests establish a future C++ comparison reference.
* `__main__.py`: deterministic original experiments and normalized-LF source
  SHA-256 hashes. No sibling-repository imports or course-source dependencies.
* Planned: conservative cell-centred heat, finite thermal contact resistance,
  overlap eigenstrain integration, and a C++17 assembly kernel with Python
  interface and standalone native tests. These are not implemented on day 1.

From the repository root after `python suite.py install`:

```sh
python -m thermomechanics --output results/thermoelastic-baseline.json
python -m pytest tests/test_thermomechanics.py -q
```

Without `--output`, the CLI prints JSON. An explicit output path is
overwritten. `solve(mesh, mean_delta_temperature, left=0, right=0)` returns
displacements, the full residual, and element forces. Set `right=None` for
a free right end. `continuum_solution(layers, left=0, right=0)` returns the
constant force and evaluates displacement at coordinates in `[0,L]`.
The dense implementation targets small verification problems; no conditioning
or scalability claim is made for extreme material contrasts or fine meshes.

## Acceptance and error separation

1. Free expansion and fully restrained patches: scaled error <=1e-10.
2. Independent Python/C++ stiffness and thermal load: scaled error <=1e-11,
   plus native invariants without Python.
3. Nonmatching overlap transfer: material-weighted integrated eigenstrain
   error <=5e-12 on a documented relative/absolute scale. Compare naive point
   interpolation, recording outcomes even when it performs similarly.
4. Reactions match the continuum identity for **supplied** temperature fields.
   This cannot establish the accuracy of the future heat solver.
5. Four-level temporal, thermal spatial and transfer-only studies where
   applicable: target first-order backward Euler and second-order smooth
   thermal spatial convergence. Investigate departures, do not tune data.
6. Independent energy balance, equilibrium, finite-contact steady flux and
   zero/large resistance limits within documented conditioning bounds.
7. Linux/Windows regression tests, C++ builds and all snapshot hashes pass;
   synthetic end-to-end experiments reproduce.
8. Replace the synthesis at `paper/paper.tex` and `paper/paper.pdf` with the
   measured extension study, <=5 pages including references; compile in CI,
   reject undefined references/layout errors, inspect all rendered pages.
9. Release only after all applicable gates pass. Day 9 blockers require a
   recovery plan and continued work on this project.

Day-1 results distinguish FE implementation agreement from continuous
displacement interpolation error. Four boundary/material cases and four
refinement levels are recorded in `results/thermoelastic-baseline.json`.
For the affine field `deltaT=80*x` K with L=1 m and alpha=23e-6/K, the
exact midpoint P1 error is `alpha*80*h^2/8`. Exact element means give exact
nodal displacement (up to roundoff); between nodes the error is second order.
This is a mechanics interpolation study, **not thermal spatial convergence**.
Roundoff metrics are platform-dependent; CI uses scales/tolerances, not
byte-for-byte floating-point JSON equality.

## Seven-day work plan and manuscript

| Day/date | Contribution and gate |
| --- | --- |
| 1 / 27 Oct | Contracts, independent continuum oracle, Python baseline and patches |
| 2 / 28 Oct | Conservative transient heat and contact; flux, energy, steady checks |
| 3 / 29 Oct | C++17 assembly, Python interface and native invariants |
| 4 / 30 Oct | Nonmatching overlap transfer and coupled CLI; failure cases |
| 5 / 31 Oct | Separate refinement studies, contact sweeps, Linux/Windows evidence |
| 6 / 1 Nov | Reproduction, original <=5-page synthesis and full visual review |
| 7 / 2 Nov | Acceptance review, tested merge/release, next coherent question |
| 8-9 / 3-4 Nov | Justified repairs only; report blockers and retain the project |

The planned synthesis will connect the heat resistance network, eigenstrain
integrals and restrained force identity; measured tables will separate heat,
transfer and mechanics errors. There is no new-paper completion claim today.

## Topic provenance, licensing and derivation

DTU course 41747, **Modelling in Materials and Manufacturing Engineering**
(MMME), Lecture 4 topics reviewed read-only on 27 October: explicit heat
positivity, resistance-based interface reconstruction and the initial-time
limit of matched interface weights. TUM finite-element topics inform the P1
weak form. No course files, exercise parameters, teaching code, solutions,
exams or datasets are included. The public repository records topics only.
The equations above are derived from equilibrium, the constitutive law and
integrated compatibility; no novelty or external physical validation is claimed.
The original integration layer's MIT license applies. All six component
snapshots and their integrity/provenance manifest remain unchanged.
