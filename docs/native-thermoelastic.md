# Original C++17 thermoelastic assembly

This is day 3 of the thermoelastic extension (29 October 2026). The kernel
and interface are original MIT-licensed code, outside the six preserved
component snapshots. The 19 original component repositories are unchanged.
The research question and equations are in [the contract](thermoelastic-transfer.md).
MMME resistance-network topics motivate the subsequent thermal coupling;
TUM finite-element topics motivate the mechanics weak form. Course files,
teaching code and exercise data are not dependencies or published assets.

## Build and reproduce

After the README environment setup, run from this repository:

```sh
python suite.py native
python suite.py test
python -m thermomechanics.native_experiment
```

The native command builds four preserved component invariant executables and
the new shared library plus standalone thermoelastic CTest. On Windows, a
short `--build-root` avoids MSBuild path limits. CMake records the exact Release
library path; `suite.py` copies it to ignored `build/native-library.txt`.
Do not copy that machine-specific record between machines. Rebuild instead.
Linux and Windows CI build before testing; missing libraries fail the native
Python tests, rather than silently selecting the Python backend.

Standalone builds need only CMake >=3.20 and a C++17 compiler:

```sh
cmake -S native -B build/thermo -DCMAKE_BUILD_TYPE=Release
cmake --build build/thermo --config Release
ctest --test-dir build/thermo -C Release --output-on-failure
```

For standalone builds, pass the library path in `library-Release.txt` to
`NativeAssembler(path)` or set `THERMOMECHANICS_LIBRARY` to that path.
Explicit constructor paths take priority, then the environment, then the
suite build record. No directory scanning or automatic compilation occurs.

```python
from thermomechanics import BarMesh, Material, solve
from thermomechanics.native import NativeAssembler

native = NativeAssembler()
mesh = BarMesh((0., .4, 1.), (Material(70e9, 23e-6), Material(190e9, 11e-6)),
               (2e-4, 1.5e-4))
k, thermal_load = native.assemble(mesh, [60., 25.])
result = solve(mesh, [60., 25.], left=0., right=0., assembler=native.assemble)
```

`solve` retains the independent Python assembler by default. Supplying the
native callable replaces assembly only: constrained-node elimination, dense
linear solution, and force recovery still use Python/NumPy. This is not an
independent C++ linear solver or a claimed performance improvement.

## Kernel and interface contract

For each element, the C++ loop computes `k=EA/h` and
`f=EA*alpha*mean(T-Tref)`. It scatters `(k,k)` to the global diagonal, `-k`
to the symmetric off-diagonal, and `(-f,+f)` to the thermal load. SI units
and tension-positive reactions match the Python contract. Every material
or area discontinuity must be a mechanical node. Expansion coefficients
may be negative; modulus, area and element length must be positive.

`native/assembly.h` defines C ABI version 1. All input/output buffers have
documented lengths, must be valid and nonoverlapping, and belong to the
caller. The C API cannot determine the allocation extent behind a raw
pointer. The Python wrapper validates temperature shape and provides owned,
contiguous float64 buffers, explicit argument/result types, and a retained
library handle. Status 1/2 rejects counts, pointers or nonfinite/invalid
inputs; status 3 rejects nonrepresentable arithmetic. Outputs must be
discarded on error. The wrapper translates status codes to Python exceptions.
No C++ allocation or exception crosses the C boundary. Finite but badly
conditioned systems and tiny thermal-load underflow are not diagnosed.

Native assembly has O(n) storage and work; the wrapper reconstructs a dense
matrix with O(n^2) storage. This implementation targets small verification
problems; no scalability, robustness at arbitrary property contrasts or
performance claim is made. The C++ library has no mutable global state.

The FFI follows the official [ctypes documentation](https://docs.python.org/3/library/ctypes.html),
including explicit function signatures and the C calling convention.
The build uses CMake's [shared library target](https://cmake.org/cmake/help/latest/command/add_library.html).

## Verification and evidence

The standalone Release CTest uses runtime checks, not disabled `assert`s.
It checks hand-calculated heterogeneous assembly, signed support reactions,
rigid translation, elastic energy, free expansion, restraint, negative/zero
expansion, invalid inputs and arithmetic failures. Python checks add an
independent finite-difference energy gradient and noncontiguous input handling.

The [seeded experiment](../thermomechanics/native_experiment.py) runs 36
original synthetic cases: 1/2/7/24 elements, three heterogeneous realizations,
and free/fixed/imposed right displacement. Layer temperatures are affine;
exact element means feed assembly. Native stiffness/load are compared with
Python; forces, reactions and nodal displacements are compared with the
independent analytical layer-integral oracle. Both prescribed displacements
are nonzero in the imposed-end case, checking elimination of constrained terms.

[Versioned results](../results/native-thermoelastic.json) contain normalized-LF
SHA-256 hashes for the Python and C++ source inputs. CI rejects stale hashes
and scaled errors above 1e-11; it permits bounded platform variation in
roundoff metrics. The earlier baseline JSON was regenerated for the added
source hashes; its numerical results are unchanged.

These checks establish assembly consistency and analytical agreement for
supplied element temperatures. They do not measure heat discretization,
nonmatching-transfer error, between-node displacement error or physical
validity. Those are separate studies. The released two-page beam paper is
unchanged; the extension synthesis remains scheduled after coupled results.
