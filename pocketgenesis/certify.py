"""UNTRUSTED certificate construction (solver side). Outputs are verified by checker.verify."""
from __future__ import annotations

from fractions import Fraction
import itertools
import random

import numpy as np

from . import checker as ck
from .regular import AlphaComplex, exact_witness
from .escape import Grid, clearance, compress_polyline


def fq(x: float, den: int = 10**6) -> Fraction:
    """Rational approximation of a float on a fixed grid (exact thereafter)."""
    return Fraction(int(round(float(x) * den)), den)


def s(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else str(fr.numerator)


def sv(v) -> list:
    return [s(Fraction(t)) for t in v]


def atoms_rational(X, r, den: int = 1000):
    return [(fq(p[0], den), fq(p[1], den), fq(p[2], den), fq(ri, den)) for p, ri in zip(X, r)]


def atoms_json(Aq):
    return [sv(a) for a in Aq]


def float_atoms(Aq):
    X = np.array([[float(a[0]), float(a[1]), float(a[2])] for a in Aq])
    r = np.array([float(a[3]) for a in Aq])
    return X, r


def best_chamber_point(X, r, cands, wall: float = -np.inf, iters=300, step=0.25, seed=0):
    """Pattern search maximising c from candidate seeds (float, untrusted). A move is accepted
    only if sampled clearance along the move stays above `wall`, so the search does not hop
    across the barrier {c <= wall} (soundness never depends on this)."""
    rng = np.random.default_rng(seed)
    cands = list(cands)
    if not cands:
        return None, -np.inf
    cvals = clearance(X, r, np.array(cands))
    order = np.argsort(-cvals)[:8]
    best_p, best_c = None, -np.inf
    ts = np.linspace(0, 1, 9)[1:]
    for j0 in order:
        p = np.array(cands[j0], float)
        cp = cvals[j0]
        if cp <= wall:
            continue
        st = step
        for _ in range(iters):
            trial = p + rng.normal(size=(12, 3)) * st
            ct = clearance(X, r, trial)
            j = int(np.argmax(ct))
            if ct[j] > cp:
                seg = p[None] + ts[:, None] * (trial[j] - p)[None]
                if clearance(X, r, seg).min() > wall:
                    p, cp = trial[j], ct[j]
                    continue
            st *= 0.85
            if st < 1e-4:
                break
        if cp > best_c:
            best_p, best_c = p, cp
    return best_p, best_c


def generic_ray(Pq, T, y, seed=0, tries=200):
    rnd = random.Random(seed)
    for _ in range(tries):
        u = (Fraction(rnd.randint(-997, 997)), Fraction(rnd.randint(-997, 997)), Fraction(rnd.randint(-997, 997)))
        if ck.norm2(u) == 0:
            continue
        par = ck.winding_parity(Pq, T, y, u)
        if par is not None:
            return u, par
    return None, None


def cavity_pocket_certificates(Aq, b: Fraction, a_min: float | None = None, eps=None, eta=None,
                               shrink=None, den=10**4, max_cavities=50):
    """Certificates of P_{a,b} for every cavity of U_b (Theorem B4/B5, R2 when eps/eta given).

    The alpha complex is computed for radii r_i + b - tol_i (shrunken balls, R2), so that every
    Cech witness also satisfies the robust inequality.  Returns list of (cert, verdict).
    """
    X, r = float_atoms(Aq)
    n = len(Aq)
    tol = np.zeros(n)
    if eps is not None:
        tol += np.array([float(e) for e in eps])
    if eta is not None:
        tol += np.array([float(e) for e in eta])
    Rb = r + float(b) - tol
    if (Rb <= 0).any():
        raise ValueError("tolerances exceed inflated radii")
    alpha = AlphaComplex(X, Rb)
    Pq = [a[:3] for a in Aq]
    rq = [a[3] for a in Aq]
    tolq = [Fraction(0)] * n
    if eps is not None:
        tolq = [tolq[i] + Fraction(eps[i]) for i in range(n)]
    if eta is not None:
        tolq = [tolq[i] + Fraction(eta[i]) for i in range(n)]
    wq = [(rq[i] + b - tolq[i]) ** 2 for i in range(n)]
    out = []
    for cav in alpha.cavities[:max_cavities]:
        T = alpha.cavity_boundary(cav)
        if not T:
            continue
        seeds = [alpha.tet_center[t] for t in cav] + [X[list(alpha.tets[t])].mean(0) for t in cav]
        # keep seeds inside the cavity's tetrahedra hull region roughly: use centroids of cavity tets
        yf, cy = best_chamber_point(X, r, seeds, wall=float(b))
        if yf is None:
            continue
        y = tuple(fq(v, den) for v in yf)
        # chamber value: largest rational a (on the den grid) with c(y) >= a + tol robustly
        cyv = float(clearance(X, r, np.array([[float(v) for v in y]]))[0])
        a = Fraction(int(np.floor((cyv - float(max(tol, default=0)) - 1e-6) * den)), den)
        if a_min is not None and float(a) < a_min:
            continue
        if not a > b:
            continue
        W = []
        okw = True
        for tri in T:
            z = exact_witness(Pq, wq, tri, alpha)
            if z is None:
                okw = False
                break
            W.append(z)
        if not okw:
            out.append((None, {"verdict": "UNKNOWN", "reason": "no exact witness for some triangle"}))
            continue
        u, par = generic_ray(Pq, T, y)
        if u is None:
            out.append((None, {"verdict": "UNKNOWN", "reason": "no generic ray"}))
            continue
        cert = {"atoms": atoms_json(Aq),
                "claim": {"type": "pocket", "a": s(a), "b": s(b), "y": sv(y),
                          "barrier": {"type": "cech", "T": [list(map(int, t)) for t in T],
                                      "W": [sv(z) for z in W], "ray": sv(u)}}}
        if eps is not None:
            cert["eps"] = [s(Fraction(e)) for e in eps]
        if eta is not None:
            cert["eta"] = [s(Fraction(e)) for e in eta]
        out.append((cert, ck.verify(cert)))
    return out


def escape_lower_certificate(Aq, y, level_guess: float, grid: Grid | None = None, h=0.4, pad=6.0,
                             eps=None, bulk=None, den=10**4, decrements=8):
    """Path certificate E(y) >= b' (B2/R1b) from a grid maximin path; b' lowered until valid."""
    X, r = float_atoms(Aq)
    grid = grid or Grid(X, r, h, pad)
    o = X.mean(0) if bulk is None else np.asarray(bulk[0], float)
    tol = 0.0 if eps is None else max(float(e) for e in eps)
    rho = float(np.max(np.linalg.norm(X - o, axis=1))) + tol
    oq = tuple(fq(v, den) for v in o)
    Rq = fq(rho + 1.0, den)
    for k in range(decrements):
        lev = level_guess - k * grid.h
        pts = grid.maximin_path(np.array([float(v) for v in y]), lev)
        if pts is None:
            continue
        pts = compress_polyline(pts)
        # continue radially from the last grid point to the bulk sphere
        last = pts[-1]
        dirv = last - o
        dirv = dirv / (np.linalg.norm(dirv) + 1e-12)
        # F1(a): beyond the centres' radius, moving outward never decreases clearance
        extra = max(0.0, (rho + 2.0) - float(np.linalg.norm(last - o)))
        far = last + dirv * extra
        poly = [tuple(Fraction(v) for v in y)] + [tuple(fq(v, den) for v in p) for p in pts[1:]]
        if extra > 0:
            poly.append(tuple(fq(v, den) for v in far))
        # exact level: try lev - h (segments may dip by h/2 below node values)
        b_try = fq(lev - grid.h, den)
        cert = {"atoms": atoms_json(Aq),
                "claim": {"type": "escape_lower", "y": sv(y), "b": s(b_try), "path": [sv(p) for p in poly],
                          "bulk": {"o": sv(oq), "R": s(Rq)}}}
        if eps is not None:
            cert["eps"] = [s(Fraction(e)) for e in eps]
        v = ck.verify(cert)
        if v["verdict"] == "PROVED":
            return cert, v
    return None, {"verdict": "UNKNOWN", "reason": "no valid path found"}


def seg_clearance(X, r, p, q):
    """Float clearance min_i (dist(x_i,[p,q]) - r_i) of a segment (closed form)."""
    d = q - p
    dd = d @ d
    if dd == 0:
        return float((np.linalg.norm(X - p, axis=1) - r).min())
    t = np.clip(((X - p) @ d) / dd, 0.0, 1.0)
    proj = p[None] + t[:, None] * d[None]
    return float((np.linalg.norm(X - proj, axis=1) - r).min())


def optimize_path(X, r, poly, iters=400, step=0.3, seed=0):
    """Increase the bottleneck min_k segclear(seg_k) of a polygon by local moves of interior
    vertices and greedy shortcutting (untrusted; the result is re-verified exactly)."""
    rng = np.random.default_rng(seed)
    P = [np.array(p, float) for p in poly]

    def bott(P):
        return min(seg_clearance(X, r, P[k], P[k + 1]) for k in range(len(P) - 1))

    # shortcutting
    k = 0
    while k < len(P) - 2:
        cur = min(seg_clearance(X, r, P[k], P[k + 1]), seg_clearance(X, r, P[k + 1], P[k + 2]))
        if seg_clearance(X, r, P[k], P[k + 2]) >= cur - 1e-12:
            del P[k + 1]
        else:
            k += 1
    for it in range(iters):
        segs = [seg_clearance(X, r, P[k], P[k + 1]) for k in range(len(P) - 1)]
        kmin = int(np.argmin(segs))
        improved = False
        for v in (kmin, kmin + 1):
            if v == 0 or v == len(P) - 1:
                continue
            for _ in range(8):
                cand = P[v] + rng.normal(size=3) * step
                old = min(segs[v - 1], segs[v])
                new = min(seg_clearance(X, r, P[v - 1], cand), seg_clearance(X, r, cand, P[v + 1]))
                if new > old:
                    P[v] = cand
                    improved = True
                    break
        if not improved:
            # split the bottleneck segment to add freedom
            mid = 0.5 * (P[kmin] + P[kmin + 1])
            P.insert(kmin + 1, mid)
            step *= 0.9
        if step < 1e-4:
            break
    return P, bott(P)


def escape_lower_certificate_opt(Aq, y, level_guess: float, h=0.4, pad=6.0, eps=None, den=10**6,
                                 iters=600):
    """Grid path -> optimized polygon -> exact certificate at the optimized bottleneck."""
    cert0, v0 = escape_lower_certificate(Aq, y, level_guess, h=h, pad=pad, eps=eps, decrements=30)
    if cert0 is None:
        return None, v0
    X, r = float_atoms(Aq)
    tol = np.zeros(len(r)) if eps is None else np.array([float(e) for e in eps])
    poly = [np.array([float(c) for c in ck.qv(p)]) for p in cert0["claim"]["path"]]
    P, bt = optimize_path(X, r + tol, poly, iters=iters)
    P[0] = np.array([float(v) for v in y])
    pq = [tuple(Fraction(v) for v in y)] + [tuple(fq(v, den) for v in p) for p in P[1:]]
    best = (cert0, v0)
    for slack in (1e-6, 1e-4, 1e-3, 1e-2):
        b_try = fq(bt - slack, den)
        if b_try <= Fraction(cert0["claim"]["b"]):
            break
        cert = dict(cert0)
        cert["claim"] = dict(cert0["claim"], b=s(b_try), path=[sv(p) for p in pq])
        v = ck.verify(cert)
        if v["verdict"] == "PROVED":
            return cert, v
    return best
