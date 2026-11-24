"""Synthetic geometries with analytically known answers (regression tests and examples)."""
from __future__ import annotations

from fractions import Fraction
import math

import numpy as np


def octahedral_cage(rho=10.0, r=1.5):
    X = np.array([[rho, 0, 0], [-rho, 0, 0], [0, rho, 0], [0, -rho, 0], [0, 0, rho], [0, 0, -rho]], float)
    return X, np.full(6, r)


def octahedral_gate(rho=10.0, r=1.5):
    """Exact escape clearance of the cage centre: rho*sqrt(2/3) - r (Example S1)."""
    return rho * math.sqrt(2.0 / 3.0) - r


def icosahedron_dirs():
    phi = (1 + 5 ** 0.5) / 2
    V = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [(0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)]
    V = np.array(V, float)
    return V / np.linalg.norm(V, axis=1)[:, None]


ICOSA_COVER_DEG = 37.37736814064969


def icosahedral_cage(rho=10.0, r=1.5):
    return icosahedron_dirs() * rho, np.full(12, r)


def icosahedral_gate(rho=10.0, r=1.5):
    return rho * math.sin(math.radians(ICOSA_COVER_DEG)) - r


def example_S1(eps=0.0, rho=10.0, r=1.5, d=9.0, rs=1.5):
    """Octahedral cage + chamber-control atom; eps>0 gives the perturbed configuration Y."""
    X, rad = octahedral_cage(rho, r)
    u = np.ones(3) / math.sqrt(3)
    X = np.vstack([X, -(d + eps) * u])
    rad = np.concatenate([rad + eps, [rs]])
    return X, rad


def covering_construction(dirs: np.ndarray, a: float, r: float):
    """Theorem C4 (<=): balls of radius r at distance a + r along covering directions."""
    return dirs * (a + r), np.full(len(dirs), r)


def fejes_toth_min(theta_rad: float) -> int:
    """Least N allowed by Fejes Toth's covering bound for caps of angular radius theta."""
    c = math.cos(theta_rad)
    for N in range(3, 10 ** 6):
        w = N * math.pi / (6 * (N - 2))
        if c <= (1 / math.tan(w)) / math.sqrt(3) + 1e-15:
            return N
    raise ValueError


def theta_star(a: float, b: float, r: float) -> float:
    return math.asin((r + b) / (a + r))


def np_gadget(clauses, nvars, M=24):
    """Theorem C5 reduction: returns (fixed shell atoms, groups, windows, a, b) with integer data.

    clauses: list of tuples of nonzero ints (DIMACS literals). Shell: unit balls at even lattice
    points on the boundary of [0,2M]^3 except windows on z = 0. Each group j (variable j+1) has
    two states; each state is a list of ball centres (radius 1)."""
    assert M % 2 == 0
    m = len(clauses)
    slots = []
    for x in range(8, 2 * M - 7, 8):
        for y in range(8, 2 * M - 7, 8):
            slots.append((x, y, 0))
    assert len(slots) >= m, "increase M"
    windows = slots[:m]
    wset = set(windows)
    shell = []
    for x in range(0, 2 * M + 1, 2):
        for y in range(0, 2 * M + 1, 2):
            for z in range(0, 2 * M + 1, 2):
                if x in (0, 2 * M) or y in (0, 2 * M) or z in (0, 2 * M):
                    if (x, y, z) not in wset:
                        shell.append((x, y, z))
    park = [p for p in shell if p[2] == 2 * M]        # top face: harmless duplicate positions
    groups = []
    pk = 0
    for j in range(1, nvars + 1):
        T, F = [], []
        for ci, cl in enumerate(clauses):
            for lit in cl:
                if abs(lit) != j:
                    continue
                w = windows[ci]
                p = park[pk % len(park)]
                pk += 1
                if lit > 0:
                    T.append(w)
                    F.append(p)
                else:
                    T.append(p)
                    F.append(w)
        groups.append((T, F))
    a = Fraction(M - 2)
    b = Fraction(1, 2)
    return shell, groups, windows, a, b


def np_gadget_config(shell, groups, assignment):
    pts = list(shell)
    for (T, F), val in zip(groups, assignment):
        pts += T if val else F
    X = np.array(pts, float)
    return X, np.ones(len(pts))
