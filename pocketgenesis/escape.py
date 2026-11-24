"""UNTRUSTED grid engine: clearance, escape clearance, basins, maximin paths.

Union-find flooding over a regular grid processed in decreasing clearance (merge-tree
construction). Gives E_grid, D_grid, and all basins (Theorem F6) in O(N log N). Output is
advisory: certificates built from it are verified exactly by pocketgenesis.checker.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def clearance(X: np.ndarray, r: np.ndarray, Y: np.ndarray, k: int = 24, tree=None) -> np.ndarray:
    """c_X at points Y (float). Uses the k nearest centres, exact if radii spread is small
    relative to inter-atomic spacing (and always an upper bound of the true minimum)."""
    tree = tree or cKDTree(X)
    k = min(k, len(X))
    d, idx = tree.query(Y, k=k)
    if k == 1:
        d, idx = d[:, None], idx[:, None]
    return (d - r[idx]).min(1)


class Grid:
    def __init__(self, X: np.ndarray, r: np.ndarray, h: float, pad: float, lo=None, hi=None):
        self.X, self.r, self.h = np.asarray(X, float), np.asarray(r, float), float(h)
        lo = np.min(self.X - self.r[:, None], 0) - pad if lo is None else np.asarray(lo, float)
        hi = np.max(self.X + self.r[:, None], 0) + pad if hi is None else np.asarray(hi, float)
        self.shape = tuple(int(np.ceil((hi[d] - lo[d]) / h)) + 1 for d in range(3))
        self.lo = lo
        axes = [lo[d] + h * np.arange(self.shape[d]) for d in range(3)]
        G = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
        self.points = G
        self.tree = cKDTree(self.X)
        self.c = clearance(self.X, self.r, G, tree=self.tree)

    def index(self, p) -> int:
        ijk = np.clip(np.rint((np.asarray(p) - self.lo) / self.h).astype(int), 0, np.array(self.shape) - 1)
        return int(np.ravel_multi_index(tuple(ijk), self.shape))

    def _prep(self):
        if hasattr(self, "_ijk"):
            return
        N = len(self.c)
        self._ijk = np.stack(np.unravel_index(np.arange(N), self.shape), 1).astype(np.int64)
        sh = self.shape
        self._strides = (sh[1] * sh[2], sh[2], 1)
        self._bnd = ((self._ijk == 0) | (self._ijk == np.array(sh) - 1)).any(1)

    def neighbors(self, v: int):
        self._prep()
        ijk = self._ijk[v]
        for d in range(3):
            st = self._strides[d]
            if ijk[d] > 0:
                yield v - st
            if ijk[d] < self.shape[d] - 1:
                yield v + st

    def is_boundary(self, v: int) -> bool:
        self._prep()
        return bool(self._bnd[v])

    def flood(self):
        """Compute E_grid for every node and the list of basins (max c, mouth e, argmax, size)."""
        N = len(self.c)
        order = np.argsort(-self.c, kind="stable")
        parent = np.full(N, -1, dtype=np.int64)     # -1 = not yet inserted
        bulk = np.zeros(N, dtype=bool)              # per root: connected to the grid boundary
        members: dict[int, list[int]] = {}
        topnode: dict[int, int] = {}
        E = np.full(N, np.nan)
        basins = []

        def find(a):
            root = a
            while parent[root] != root:
                root = parent[root]
            while parent[a] != root:
                parent[a], a = root, parent[a]
            return root

        sh = self.shape
        for v in order:
            v = int(v)
            cv = self.c[v]
            parent[v] = v
            members[v] = [v]
            topnode[v] = v
            bulk[v] = self.is_boundary(v)
            roots = {v}
            for u in self.neighbors(v):
                if parent[u] >= 0:
                    roots.add(find(u))
            bulk_any = any(bulk[rt] for rt in roots)
            # merge all into the largest
            roots = list(roots)
            if bulk_any:
                for rt in roots:
                    if not bulk[rt]:
                        mem = members.pop(rt, [])
                        if rt != v and mem:
                            basins.append((float(self.c[topnode[rt]]), float(cv), int(topnode[rt]), len(mem)))
                        E[mem] = cv
                    else:
                        members.pop(rt, None)
                main = roots[0]
                for rt in roots[1:]:
                    parent[rt] = main
                bulk[main] = True
                if np.isnan(E[v]):
                    E[v] = cv
            else:
                roots.sort(key=lambda t: -len(members.get(t, [])))
                main = roots[0]
                for rt in roots[1:]:
                    parent[rt] = main
                    members[main].extend(members.pop(rt))
                    if self.c[topnode[rt]] > self.c[topnode[main]]:
                        topnode[main] = topnode[rt]
        self.E = E
        self.D = self.c - E
        self.basins = sorted(basins, key=lambda b: -(b[0] - b[1]))
        return self.basins

    def maximin_path(self, y, level: float):
        """BFS path (grid nodes) from y's node to the boundary inside {c >= level}."""
        from collections import deque
        s = self.index(y)
        if self.c[s] < level:
            return None
        prev = {s: -1}
        dq = deque([s])
        while dq:
            v = dq.popleft()
            if self.is_boundary(v):
                path = []
                while v != -1:
                    path.append(v)
                    v = prev[v]
                return [self.points[p] for p in reversed(path)]
            for u in self.neighbors(v):
                if u not in prev and self.c[u] >= level:
                    prev[u] = v
                    dq.append(u)
        return None


def compress_polyline(pts):
    """Merge collinear consecutive grid steps (the polygon is unchanged as a set)."""
    if len(pts) <= 2:
        return pts
    out = [pts[0]]
    for k in range(1, len(pts) - 1):
        d1 = pts[k] - out[-1]
        d2 = pts[k + 1] - pts[k]
        if np.linalg.norm(np.cross(d1, d2)) > 1e-9 or d1 @ d2 < 0:
            out.append(pts[k])
    out.append(pts[-1])
    return out
