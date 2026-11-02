"""Explicit C++17 backend via a versioned C ABI; no automatic Python fallback."""
import ctypes
import os
from pathlib import Path
import numpy as np
from .model import element_temperatures


class NativeAssembler:
    """Load an explicit library, environment path, or suite.py build record.

    Instances own the library handle. Calls allocate independent contiguous
    float64 buffers; the C++ kernel has no mutable global state. Assembly is
    O(n) natively; reconstruction of the Python dense matrix is O(n^2).
    """

    def __init__(self, library=None):
        path = library if library is not None else os.environ.get('THERMOMECHANICS_LIBRARY')
        if not path:
            record = Path(__file__).resolve().parents[1]/'build/native-library.txt'
            if record.exists():
                path = record.read_text(encoding='utf-8').strip()
        if not path:
            raise RuntimeError('Build native/CMakeLists.txt and set THERMOMECHANICS_LIBRARY or pass library=path')
        self.path = Path(path).resolve(strict=True)
        self._library = ctypes.CDLL(str(self.path))
        version = self._library.tm_abi_version
        version.argtypes = []
        version.restype = ctypes.c_int
        if version() != 1:
            raise RuntimeError('Unsupported thermoelastic native ABI (expected 1)')
        array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED'))
        self._assemble = self._library.tm_assemble
        self._assemble.argtypes = [ctypes.c_size_t] + [array]*8
        self._assemble.restype = ctypes.c_int

    def assemble(self, mesh, delta_temperature):
        dt = element_temperatures(mesh, delta_temperature)
        # Native calls retain only independently owned input buffers.
        arrays = [np.array(values, dtype=np.float64, order='C', copy=True) for values in
                  (mesh.nodes, [m.young for m in mesh.materials], mesh.areas,
                   [m.expansion for m in mesh.materials], dt)]
        n = len(mesh.materials)
        diagonal, off_diagonal, load = np.empty(n+1), np.empty(n), np.empty(n+1)
        status = self._assemble(n, *arrays, diagonal, off_diagonal, load)
        if status in (1, 2):
            raise ValueError(f'Native thermoelastic assembly rejected input (status {status})')
        if status == 3:
            raise FloatingPointError('Native thermoelastic arithmetic overflow/underflow')
        if status:
            raise RuntimeError(f'Unknown native assembly status {status}')
        stiffness = np.diag(diagonal) + np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
        return stiffness, load
