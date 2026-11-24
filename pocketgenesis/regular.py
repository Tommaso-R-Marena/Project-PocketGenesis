"""UNTRUSTED solver: regular (weighted Delaunay) triangulation, alpha complex, cavities.

Floating point is used freely here; every output that matters is re-verified by the exact
checker (``pocketgenesis.checker``). Witness points for Cech triangles are recomputed exactly
in rational arithmetic from the combinatorial structure.

Theory: Edelsbrunner, "The union of balls and its dual shape" (1995); Theorem B5 in
mathematics/proofs/P3_duality_certificates.md.
"""
from __future__ import annotations

from fractions import Fraction
import itertools

import numpy as np
from scipy.spatial import ConvexHull

from . import checker as ck


def regular_tetrahedra(P: np.ndarray, w: np.ndarray):
    """Tetrahedra of the regular triangulation of weighted points (P_i, w_i) (power distance
    ||z-P_i||^2 - w_i), via the lower convex hull of the lifted points (P_i, |P_i|^2 - w_i).

    An apex far above the lifted set makes the hull full-dimensional even for cospherical,
    equal-weight inputs; facets touching the apex are discarded. Falls back to joggling (QJ),
    which yields a triangulation of a perturbed input -- acceptable because every certificate
    is re-verified exactly on the original data."""
    n = len(P)
    lifted = np.column_stack([P, (P * P).sum(1) - w])
    c = lifted.mean(0)
    L = lifted - c
    span = np.ptp(L, axis=0).max() + 1.0
    apex = np.zeros(4)
    apex[3] = 10.0 * span * span + 10.0 * span
    pts = np.vstack([L, apex])
    tets = None
    for opts in ("Qt Qc", "QJ"):
        try:
            hull = ConvexHull(pts, qhull_options=opts)
        except Exception:
            continue
        simp = hull.simplices
        keep = (simp < n).all(1) & (hull.equations[:, 3] < 0)
        tets = np.sort(simp[keep], axis=1)
        break
    if tets is None:
        raise RuntimeError("regular triangulation failed")
    v = np.abs(np.einsum("ij,ij->i", np.cross(P[tets[:, 1]] - P[tets[:, 0]], P[tets[:, 2]] - P[tets[:, 0]]),
                         P[tets[:, 3]] - P[tets[:, 0]]))
    return tets[v > 1e-10]


def orthocenter_float(P: np.ndarray, w: np.ndarray, idx) -> tuple[np.ndarray, float]:
    """Power centre of the weighted points idx inside their affine hull, and its power."""
    idx = list(idx)
    x0 = P[idx[0]]
    D = np.array([P[k] - x0 for k in idx[1:]])          # m x 3
    rhs = np.array([0.5 * ((P[k] @ P[k]) - (x0 @ x0) + w[idx[0]] - w[k]) for k in idx[1:]])
    # z = x0 + D^T lam ; condition D z = rhs  ->  (D D^T) lam = rhs - D x0
    G = D @ D.T
    lam = np.linalg.solve(G, rhs - D @ x0)
    z = x0 + D.T @ lam
    return z, float((z - x0) @ (z - x0) - w[idx[0]])


def orthocenter_exact(Pq, wq, idx):
    """Exact rational power centre in the affine hull (same equations as above)."""
    idx = list(idx)
    x0 = Pq[idx[0]]
    D = [ck.sub(Pq[k], x0) for k in idx[1:]]
    m = len(D)
    rhs = [(ck.norm2(Pq[k]) - ck.norm2(x0) + wq[idx[0]] - wq[k]) / 2 - ck.dot(D[a], x0)
           for a, k in enumerate(idx[1:])]
    G = [[ck.dot(D[a], D[c]) for c in range(m)] for a in range(m)]
    # Gaussian elimination over Q
    M = [G[a][:] + [rhs[a]] for a in range(m)]
    for col in range(m):
        piv = next(r for r in range(col, m) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for r in range(m):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [M[r][c] - f * M[col][c] for c in range(m + 1)]
    lam = [M[a][m] / M[a][a] for a in range(m)]
    z = x0
    for a in range(m):
        z = ck.add(z, ck.scale(lam[a], D[a]))
    return z


class AlphaComplex:
    """Alpha complex at alpha = 0 of balls B(x_i, R_i) (weights R_i^2) and its cavities."""

    def __init__(self, P: np.ndarray, R: np.ndarray):
        self.P = np.asarray(P, float)
        self.R = np.asarray(R, float)
        self.w = self.R ** 2
        self.tets = regular_tetrahedra(self.P, self.w)
        self._build()

    def _build(self):
        P, w = self.P, self.w
        tri_cof: dict[tuple, list] = {}
        for t_id, t in enumerate(self.tets):
            for f in itertools.combinations(t, 3):
                tri_cof.setdefault(tuple(f), []).append(t_id)
        self.tri_cof = tri_cof
        self.tet_size = np.empty(len(self.tets))
        self.tet_center = np.empty((len(self.tets), 3))
        for t_id, t in enumerate(self.tets):
            z, s = orthocenter_float(P, w, t)
            self.tet_center[t_id] = z
            self.tet_size[t_id] = s
        # triangle filtration values f(sigma) at alpha = 0 (Edelsbrunner-Muecke attachment rule)
        self.tri_val = {}
        for tri, cof in tri_cof.items():
            z, s = orthocenter_float(P, w, tri)
            attached = False
            for t_id in cof:
                (opp,) = set(self.tets[t_id]) - set(tri)
                if (z - P[opp]) @ (z - P[opp]) - w[opp] < s - 1e-12:
                    attached = True
            if attached:
                self.tri_val[tri] = min(self.tet_size[t_id] for t_id in cof)
            else:
                self.tri_val[tri] = s
            # a triangle cannot enter before its cofaces' values imply... (monotone closure)
        self.tet_in = self.tet_size <= 0.0
        self.tri_in = {tri: (v <= 0.0) for tri, v in self.tri_val.items()}
        # union-find over tetrahedra not in K; node -1 = outside
        parent = {}

        def find(a):
            while parent.setdefault(a, a) != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        OUT = -1
        find(OUT)
        for t_id in range(len(self.tets)):
            if not self.tet_in[t_id]:
                find(t_id)
        for tri, cof in tri_cof.items():
            if self.tri_in[tri]:
                continue
            cells = [c for c in cof if not self.tet_in[c]]
            if len(cof) == 1:
                cells.append(OUT)        # hull triangle: other side is the outside
            for c in cells[1:]:
                union(cells[0], c)
        self._find = find
        comps: dict[int, list[int]] = {}
        for t_id in range(len(self.tets)):
            if not self.tet_in[t_id]:
                comps.setdefault(find(t_id), []).append(t_id)
        out_root = find(OUT)
        self.cavities = [v for k, v in comps.items() if k != out_root]

    def cavity_boundary(self, cav):
        """Triangles with exactly one side in the cavity (its boundary 2-cycle)."""
        cs = set(cav)
        bd = []
        for tri, cof in self.tri_cof.items():
            inside = sum(1 for c in cof if c in cs)
            if inside == 1:
                bd.append(tri)
        return bd

    def locate(self, y: np.ndarray):
        """Index of a tetrahedron containing y (barycentric test), or None."""
        P = self.P
        T = self.tets
        A = P[T[:, 0]]
        M = np.stack([P[T[:, 1]] - A, P[T[:, 2]] - A, P[T[:, 3]] - A], axis=2)
        try:
            lam = np.linalg.solve(M, (y - A)[:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            return None
        ok = (lam >= -1e-12).all(1) & (lam.sum(1) <= 1 + 1e-12)
        hits = np.nonzero(ok)[0]
        return int(hits[0]) if len(hits) else None


def exact_witness(Pq, wq, tri, alpha: AlphaComplex):
    """A rational point in the three balls of `tri` (Cech witness), or None."""
    cands = [tuple(tri)]
    for t_id in alpha.tri_cof.get(tuple(tri), []):
        cands.append(tuple(int(v) for v in alpha.tets[t_id]))
    for idx in cands:
        z = orthocenter_exact(Pq, wq, idx)
        if all(ck.norm2(ck.sub(z, Pq[i])) <= wq[i] for i in tri):
            return z
    return None
