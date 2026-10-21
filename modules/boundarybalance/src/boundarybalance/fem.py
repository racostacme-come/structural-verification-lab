"""P1 Poisson assembly, strong/penalty/Nitsche solves and continuum norms.

All quantities are dimensionless. Callback inputs have shape (points, 2).
Scalar callbacks return (points,); gradient callbacks return (points, 2).
The supported geometry is a uniformly triangulated unit square.
"""

from collections.abc import Callable
from dataclasses import dataclass
from numbers import Integral

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu, spsolve

Field = Callable[[np.ndarray], np.ndarray]


def _integer(value, name, minimum, maximum=None):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be <= {maximum}")


def _values(function, points, shape):
    values = np.asarray(function(points), dtype=float)
    if values.shape != shape or not np.isfinite(values).all():
        raise ValueError(f"callback must return finite values of shape {shape}")
    return values


@dataclass(frozen=True)
class Mesh:
    """CCW triangles and their outward-oriented boundary edges.

    Construct through :func:`unit_square`. Edges store (triangle, local a, local b).
    Arrays are owned by the mesh; callers must not mutate them.
    """

    points: np.ndarray
    triangles: np.ndarray
    edges: np.ndarray
    boundary: np.ndarray

    @property
    def interior(self):
        return np.setdiff1d(np.arange(len(self.points)), self.boundary)


def unit_square(n: int) -> Mesh:
    """Split each of n x n squares along its southwest-to-northeast diagonal."""
    _integer(n, "n", 2)
    axis = np.linspace(0, 1, n + 1)
    xx, yy = np.meshgrid(axis, axis)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    triangles = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, d, c = a + 1, a + n + 1, a + n + 2
            triangles.extend(((a, b, c), (a, c, d)))
    triangles = np.array(triangles, dtype=int)
    owners = {}
    for t, nodes in enumerate(triangles):
        for a, b in ((0, 1), (1, 2), (2, 0)):
            key = tuple(sorted((int(nodes[a]), int(nodes[b]))))
            if key in owners:
                del owners[key]
            else:
                owners[key] = (t, a, b)
    edges = np.array(list(owners.values()), dtype=int)
    boundary = np.unique([triangles[t, [a, b]] for t, a, b in edges])
    return Mesh(points, triangles, edges, boundary)


def triangle_rule(order=6):
    """Tensor Gauss-Legendre under a Duffy map; weights sum to 1/2."""
    _integer(order, "quadrature order", 2, 20)
    x, w = leggauss(order)
    x, w = (x + 1) / 2, w / 2
    r, s = np.meshgrid(x, x, indexing="ij")
    weights = (w[:, None] * w[None, :] * (1 - r)).ravel()
    r, s = r.ravel(), ((1 - r) * s).ravel()
    return np.column_stack((1 - r - s, r, s)), weights


def _geometry(mesh):
    vertices = mesh.points[mesh.triangles]
    jac = np.stack((vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), -1)
    det = np.linalg.det(jac)
    if np.any(det <= 0):
        raise ValueError("mesh triangles must have positive orientation and area")
    grad = np.einsum("ij,tjk->tik", np.array([[-1, -1], [1, 0], [0, 1]]), np.linalg.inv(jac))
    return vertices, det, grad


@dataclass(frozen=True)
class System:
    """Assembled operators: C_ij = integral(phi_i * normal_gradient(phi_j))."""

    mesh: Mesh
    stiffness: csr_matrix
    flux: csr_matrix
    boundary_mass: csr_matrix
    load: np.ndarray
    boundary_load: np.ndarray
    adjoint_load: np.ndarray
    prescribed: np.ndarray

    def matrix(self, method="nitsche", gamma=10.0):
        """Return a full weak matrix; strong constraints use interior elimination."""
        if method not in ("nitsche", "penalty", "strong"):
            raise ValueError("method must be nitsche, penalty or strong")
        if not np.isscalar(gamma) or not np.isfinite(gamma) or gamma <= 0:
            raise ValueError("gamma must be a finite positive scalar")
        if method == "strong":
            return self.stiffness.copy()
        matrix = self.stiffness + gamma * self.boundary_mass
        if method == "nitsche":
            matrix = matrix - self.flux - self.flux.T
        return matrix.tocsr()


def assemble(mesh: Mesh, source: Field, boundary: Field, order=6) -> System:
    """Assemble -Laplacian(u)=source with prescribed trace on the entire boundary.

    Boundary h is the edge length, NOT the triangle diameter. The matrices and
    boundary operators are exact for P1 basis functions. Callback integrals use
    the specified quadrature; callers should assess quadrature convergence.
    """
    basis, weights = triangle_rule(order)
    vertices, det, grad = _geometry(mesh)
    size = len(mesh.points)
    points = np.einsum("qi,tij->tqj", basis, vertices)
    f = _values(source, points.reshape(-1, 2), (points.shape[0] * points.shape[1],))
    f = f.reshape(points.shape[:2])
    local_load = np.einsum("tq,qi,q,t->ti", f, basis, weights, det)
    load = np.zeros(size)
    np.add.at(load, mesh.triangles.ravel(), local_load.ravel())
    local_k = det[:, None, None] / 2 * np.einsum("tik,tjk->tij", grad, grad)
    rows = np.repeat(mesh.triangles, 3, axis=1).ravel()
    cols = np.tile(mesh.triangles, (1, 3)).ravel()
    stiffness = coo_matrix((local_k.ravel(), (rows, cols)), shape=(size, size)).tocsr()
    erows, ecols, cdata, bdata = [], [], [], []
    boundary_load, adjoint_load = np.zeros(size), np.zeros(size)
    z, w = leggauss(order)
    z, w = (z + 1) / 2, w / 2
    for t, a, b in mesh.edges:
        nodes = mesh.triangles[t]
        tangent = vertices[t, b] - vertices[t, a]
        length = np.linalg.norm(tangent)
        normal = np.array([tangent[1], -tangent[0]]) / length
        normal_grad = grad[t] @ normal
        edge_basis = np.zeros((len(z), 3))
        edge_basis[:, a], edge_basis[:, b] = 1 - z, z
        g = _values(boundary, edge_basis @ vertices[t], (len(z),))
        trace_integral = length * (w @ edge_basis)
        local_c = np.outer(trace_integral, normal_grad)
        local_b = edge_basis.T @ (w[:, None] * edge_basis)
        erows.extend(np.repeat(nodes, 3))
        ecols.extend(np.tile(nodes, 3))
        cdata.extend(local_c.ravel())
        bdata.extend(local_b.ravel())
        boundary_load[nodes] += edge_basis.T @ (w * g)
        adjoint_load[nodes] -= normal_grad * length * (w @ g)
    flux = coo_matrix((cdata, (erows, ecols)), shape=(size, size)).tocsr()
    mass = coo_matrix((bdata, (erows, ecols)), shape=(size, size)).tocsr()
    prescribed = _values(boundary, mesh.points[mesh.boundary], (len(mesh.boundary),))
    return System(mesh, stiffness, flux, mass, load, boundary_load, adjoint_load, prescribed)


def critical_penalty(system: System) -> float:
    """Mesh-specific SPD threshold from an interior Schur complement.

    For D=K-C-C.T, S=D_bb-D_bi K_ii^-1 D_ib. The largest generalized
    eigenvalue of (-S,B_bb) is gamma_star. A(gamma) is SPD iff gamma>gamma_star
    in exact arithmetic. This is a floating-point diagnostic, not a certificate.
    Boundary-dense work is limited to 512 boundary nodes (n <= 128).
    """
    b, i = system.mesh.boundary, system.mesh.interior
    if len(b) > 512:
        raise ValueError("Schur diagnostic supports at most 512 boundary nodes")
    d = (system.stiffness - system.flux - system.flux.T).tocsr()
    dib = d[i][:, b].toarray()
    schur = d[b][:, b].toarray() - dib.T @ splu(d[i][:, i].tocsc()).solve(dib)
    schur = (schur + schur.T) / 2
    mass = system.boundary_mass[b][:, b].toarray()
    return float(eigh(-schur, mass, eigvals_only=True, subset_by_index=[len(b) - 1] * 2)[0])


@dataclass(frozen=True)
class Solution:
    values: np.ndarray
    backward_error: float
    gamma_star: float | None


def solve(system: System, method="nitsche", gamma=10.0) -> Solution:
    """Solve with sparse LU; Nitsche requires a margin above the measured SPD threshold.

    No result is returned for a nonfinite solution. The residual is evaluated on
    free equations only for strong constraints and on all equations otherwise.
    """
    matrix = system.matrix(method, gamma)
    rhs = system.load.copy()
    threshold = None
    if method == "strong":
        b, i = system.mesh.boundary, system.mesh.interior
        values = np.zeros(len(rhs))
        values[b] = system.prescribed
        rhs = rhs[i] - matrix[i][:, b] @ values[b]
        matrix = matrix[i][:, i]
        unknown = spsolve(matrix, rhs)
        values[i] = unknown
    else:
        rhs += gamma * system.boundary_load
        if method == "nitsche":
            threshold = critical_penalty(system)
            if gamma <= threshold + 1e-8 * max(1, abs(threshold)):
                raise ValueError(f"gamma={gamma:g} must exceed SPD threshold {threshold:.10g}")
            rhs += system.adjoint_load
        unknown = spsolve(matrix, rhs)
        values = unknown
    if not np.isfinite(values).all():
        raise ArithmeticError("nonfinite linear solve")
    scale = np.max(np.asarray(abs(matrix).sum(axis=1))) * np.linalg.norm(unknown, np.inf)
    scale += np.linalg.norm(rhs, np.inf)
    residual = np.linalg.norm(matrix @ unknown - rhs, np.inf)
    return Solution(values, float(residual / scale if scale else residual), threshold)


def errors(mesh: Mesh, values, exact: Field, gradient: Field, order=8) -> dict[str, float]:
    """Quadrature-measured continuum L2, H1 seminorm and boundary L2 errors."""
    values = np.asarray(values, dtype=float)
    if values.shape != (len(mesh.points),) or not np.isfinite(values).all():
        raise ValueError("values must contain one finite scalar per node")
    basis, weights = triangle_rule(order)
    vertices, det, grad = _geometry(mesh)
    points = np.einsum("qi,tij->tqj", basis, vertices)
    flat = points.reshape(-1, 2)
    truth = _values(exact, flat, (len(flat),)).reshape(points.shape[:2])
    truth_grad = _values(gradient, flat, flat.shape).reshape(points.shape)
    nodal = values[mesh.triangles]
    uh = np.einsum("qi,ti->tq", basis, nodal)
    duh = np.einsum("ti,tij->tj", nodal, grad)
    l2 = np.einsum("tq,q,t->", (uh - truth) ** 2, weights, det)
    h1 = np.einsum("tqj,q,t->", (duh[:, None, :] - truth_grad) ** 2, weights, det)
    z, w = leggauss(order)
    z, w = (z + 1) / 2, w / 2
    boundary_error = 0.0
    for t, a, b in mesh.edges:
        p = (1 - z[:, None]) * vertices[t, a] + z[:, None] * vertices[t, b]
        difference = (1 - z) * nodal[t, a] + z * nodal[t, b] - _values(exact, p, (len(z),))
        boundary_error += np.linalg.norm(vertices[t, b] - vertices[t, a]) * (w @ difference**2)
    return {
        "l2": float(np.sqrt(l2)),
        "h1": float(np.sqrt(h1)),
        "boundary_l2": float(np.sqrt(boundary_error)),
    }
