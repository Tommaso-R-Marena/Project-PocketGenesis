"""Trusted certificate checker for PocketGenesis (exact rational arithmetic only).

TRUSTED BASE: CPython integers and ``fractions.Fraction``, and this file. No floating point,
no square roots, no external libraries. Every check compares squared rational quantities.

Each ``check_*`` function returns ``(ok, reason)``. ``verify(cert)`` dispatches a JSON-like
certificate and returns a verdict in {"PROVED", "REFUTED", "UNKNOWN"} together with the list of
facts it established. ``UNKNOWN`` is returned whenever any check fails; it is never converted.

Mathematical references (mathematics/proofs/):
  chamber           -> D4.1, R1(a)            : c_X(y) >= a for all X in U
  path              -> B2, R1(b)              : E_X(y) >= b' for all X in U
  box barrier       -> B3, R1(c)              : E_X(y) <= b for all X in U
  Cech barrier      -> B4, R2, R3             : E_X(y) <= b for all X in U (moving barrier)
  pocket            -> chamber + barrier, a>b : P_{a,b}(X) for all X in U
  absence (octree)  -> R1, F1(c), Algorithm CL tests T1/T2 : no chamber witness for any X in U
"""
from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

Q = Fraction
Vec = tuple  # 3-tuple of Fractions


# ----------------------------------------------------------------------------- parsing

def q(v) -> Fraction:
    """Parse an exact rational. Floats are rejected on purpose."""
    if isinstance(v, bool):
        raise TypeError("bool is not a rational")
    if isinstance(v, Fraction):
        return v
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str):
        return Fraction(v.strip())
    raise TypeError(f"refusing non-exact value {v!r} of type {type(v).__name__}")


def qv(v: Sequence) -> Vec:
    if len(v) != 3:
        raise ValueError("expected a 3-vector")
    return (q(v[0]), q(v[1]), q(v[2]))


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(s: Fraction, a: Vec) -> Vec:
    return (s * a[0], s * a[1], s * a[2])


def dot(a: Vec, b: Vec) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm2(a: Vec) -> Fraction:
    return dot(a, a)


def cross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def det3(a: Vec, b: Vec, c: Vec) -> Fraction:
    return dot(a, cross(b, c))


class Model:
    """Nominal atoms X0 with per-atom position/radius tolerances (eps_i, eta_i).

    The certified statements hold for every X in U(X0; eps, eta) (DEFINITIONS D1.2).
    eps = eta = 0 gives statements about X0 itself.
    """

    def __init__(self, atoms: Iterable, eps: Iterable | None = None, eta: Iterable | None = None):
        self.x: list[Vec] = []
        self.r: list[Fraction] = []
        for a in atoms:
            self.x.append(qv(a[:3]))
            ri = q(a[3])
            if ri < 0:
                raise ValueError("negative radius")
            self.r.append(ri)
        n = len(self.x)
        self.eps = [q(e) for e in eps] if eps is not None else [Q(0)] * n
        self.eta = [q(e) for e in eta] if eta is not None else [Q(0)] * n
        if len(self.eps) != n or len(self.eta) != n:
            raise ValueError("tolerance vector length mismatch")
        if any(e < 0 for e in self.eps) or any(e < 0 for e in self.eta):
            raise ValueError("negative tolerance")
        for i in range(n):
            if self.r[i] - self.eta[i] < 0:
                raise ValueError("radius tolerance exceeds radius (radii must stay >= 0)")

    @property
    def n(self) -> int:
        return len(self.x)

    def tol(self, i: int) -> Fraction:
        return self.eps[i] + self.eta[i]


# ----------------------------------------------------------------------------- primitives

def seg_dist2(x: Vec, p: Vec, qq: Vec) -> Fraction:
    """Exact squared distance from point x to the closed segment [p, qq]."""
    d = sub(qq, p)
    dd = norm2(d)
    if dd == 0:
        return norm2(sub(x, p))
    t = dot(sub(x, p), d) / dd
    if t < 0:
        t = Q(0)
    elif t > 1:
        t = Q(1)
    return norm2(sub(x, add(p, scale(t, d))))


def box_dist2(x: Vec, lo: Vec, hi: Vec) -> Fraction:
    """Exact squared distance from x to the axis-aligned box [lo, hi]."""
    s = Q(0)
    for k in range(3):
        if x[k] < lo[k]:
            s += (lo[k] - x[k]) ** 2
        elif x[k] > hi[k]:
            s += (x[k] - hi[k]) ** 2
    return s


def ge_dist(d2: Fraction, s: Fraction) -> bool:
    """dist >= s  given d2 = dist^2 (dist >= 0)."""
    return s <= 0 or d2 >= s * s


def le_dist(d2: Fraction, s: Fraction) -> bool:
    """dist <= s  given d2 = dist^2."""
    return s >= 0 and d2 <= s * s


def resolve_point(m: Model, spec) -> tuple[Vec, Fraction]:
    """A point given explicitly, or intrinsically as an affine combination of atom centres.

    Returns (nominal point, worst-case displacement over U).  For an explicit point the
    displacement is 0; for sum_k alpha_k x_k (sum alpha = 1) it is sum |alpha_k| eps_k (R3).
    """
    if isinstance(spec, dict) and "affine" in spec:
        idx = spec["affine"]["indices"]
        co = [q(c) for c in spec["affine"]["coeffs"]]
        if len(idx) != len(co) or sum(co) != 1:
            raise ValueError("affine combination must have coefficients summing to 1")
        p = (Q(0), Q(0), Q(0))
        disp = Q(0)
        for i, c in zip(idx, co):
            p = add(p, scale(c, m.x[i]))
            disp += abs(c) * m.eps[i]
        return p, disp
    return qv(spec), Q(0)


# ----------------------------------------------------------------------------- basic checks

def check_bulk(m: Model, o: Vec, R: Fraction) -> tuple[bool, str]:
    """R >= max_i (||x_i - o|| + eps_i): the bulk sphere encloses every possible centre (F1)."""
    if R < 0:
        return False, "negative bulk radius"
    for i in range(m.n):
        s = R - m.eps[i]
        if s < 0 or norm2(sub(m.x[i], o)) > s * s:
            return False, f"bulk sphere does not enclose atom {i}"
    return True, "bulk ok"


def check_chamber(m: Model, y: Vec, a: Fraction, extra: Fraction = Q(0)) -> tuple[bool, str]:
    """c_X(y') >= a for all X in U and all y' with ||y'-y|| <= extra (R1a, R3)."""
    for i in range(m.n):
        s = a + m.r[i] + m.tol(i) + extra
        if not ge_dist(norm2(sub(y, m.x[i])), s):
            return False, f"atom {i} intrudes the chamber"
    return True, "chamber ok"


def check_path(m: Model, path: Sequence[Vec], b: Fraction, o: Vec, R: Fraction) -> tuple[bool, str]:
    """E_X(path[0]) >= b for all X in U (B2, R1b)."""
    if len(path) < 1:
        return False, "empty path"
    ok, why = check_bulk(m, o, R)
    if not ok:
        return False, why
    if norm2(sub(path[-1], o)) < R * R:
        return False, "path does not reach the bulk"
    segs = list(zip(path[:-1], path[1:])) or [(path[0], path[0])]
    for (p, pq) in segs:
        lo = tuple(min(p[k], pq[k]) for k in range(3))
        hi = tuple(max(p[k], pq[k]) for k in range(3))
        for i in range(m.n):
            s = b + m.r[i] + m.tol(i)
            if s <= 0:
                continue
            if box_dist2(m.x[i], lo, hi) >= s * s:   # exact, sound prefilter
                continue
            if seg_dist2(m.x[i], p, pq) < s * s:
                return False, f"segment {p}->{pq} too close to atom {i}"
    return True, "path ok"


def cell_corners(g: Vec, h: Fraction, k: Sequence[int]):
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                yield (g[0] + h * (k[0] + dx), g[1] + h * (k[1] + dy), g[2] + h * (k[2] + dz))


def check_box_barrier(m: Model, y: Vec, b: Fraction, g: Vec, h: Fraction,
                      A: Iterable, S: Iterable, o: Vec, R: Fraction) -> tuple[bool, str]:
    """E_X(y) <= b for all X in U (B3, R1c).  S entries are (i, j, k, atom_index)."""
    if h <= 0:
        return False, "nonpositive grid spacing"
    ok, why = check_bulk(m, o, R)
    if not ok:
        return False, why
    Aset = {tuple(int(t) for t in k) for k in A}
    Smap = {}
    for s in S:
        key = (int(s[0]), int(s[1]), int(s[2]))
        Smap[key] = int(s[3])
    if Aset & set(Smap):
        return False, "A and S overlap"
    # (i) y in some cell of A
    if not any(all(g[d] + h * k[d] <= y[d] <= g[d] + h * (k[d] + 1) for d in range(3)) for k in Aset):
        return False, "y is not in the interior region A"
    # (ii) face-closure
    for k in Aset:
        for d in range(3):
            for sgn in (-1, 1):
                nb = list(k)
                nb[d] += sgn
                nb = tuple(nb)
                if nb not in Aset and nb not in Smap:
                    return False, f"cell {k} has an unclosed face toward {nb}"
    # (iii) shell cells inside a single shrunken ball
    for key, i in Smap.items():
        if not (0 <= i < m.n):
            return False, "bad atom index in shell"
        s = m.r[i] + b - m.tol(i)
        for c in cell_corners(g, h, key):
            if not le_dist(norm2(sub(c, m.x[i])), s):
                return False, f"shell cell {key} not inside ball {i}"
    # (iv) interior region strictly inside the bulk sphere
    R2 = R * R
    for k in Aset:
        for c in cell_corners(g, h, k):
            if norm2(sub(c, o)) >= R2:
                return False, f"interior cell {k} reaches the bulk"
    return True, "box barrier ok"


def _ray_triangle(y: Vec, u: Vec, A: Vec, B: Vec, C: Vec):
    """Exact ray/triangle test. Returns 1 (proper crossing), 0 (miss) or None (degenerate)."""
    e1, e2 = sub(B, A), sub(C, A)
    nrm = cross(e1, e2)
    if norm2(nrm) == 0:
        return None   # degenerate realized triangle: rejected (conservative)
    den = dot(nrm, u)
    if den == 0:
        # ray parallel to the plane; degenerate only if it lies in the plane
        return None if dot(nrm, sub(y, A)) == 0 else 0
    t = dot(nrm, sub(A, y)) / den
    if t < 0:
        return 0
    if t == 0:
        return None
    P = add(y, scale(t, u))
    # barycentric coordinates of P (exact)
    w = sub(P, A)
    d00, d01, d11 = dot(e1, e1), dot(e1, e2), dot(e2, e2)
    d20, d21 = dot(w, e1), dot(w, e2)
    den2 = d00 * d11 - d01 * d01
    beta = (d11 * d20 - d01 * d21) / den2
    gamma = (d00 * d21 - d01 * d20) / den2
    if beta < 0 or gamma < 0 or beta + gamma > 1:
        return 0
    if beta == 0 or gamma == 0 or beta + gamma == 1:
        return None   # hits an edge or vertex: not in general position
    return 1


def winding_parity(points: Sequence[Vec], T: Sequence[Sequence[int]], y: Vec, u: Vec):
    """Mod-2 crossing number of the ray y + t u (t > 0) with the realized triangles (B4b).
    Returns 0/1, or None if the ray is not in general position."""
    if norm2(u) == 0:
        return None
    par = 0
    for t in T:
        res = _ray_triangle(y, u, points[t[0]], points[t[1]], points[t[2]])
        if res is None:
            return None
        par ^= res
    return par


def is_mod2_cycle(T: Sequence[Sequence[int]]) -> bool:
    count: dict[tuple[int, int], int] = {}
    for t in T:
        a, b2, c = sorted(int(v) for v in t)
        if a == b2 or b2 == c:
            return False
        for e in ((a, b2), (a, c), (b2, c)):
            count[e] = count.get(e, 0) + 1
    return all(v % 2 == 0 for v in count.values())


def check_cech_barrier(m: Model, y_spec, b: Fraction, a: Fraction, T, W, u: Vec) -> tuple[bool, str]:
    """E_X(y) <= b and c_X(y) >= a for all X in U, via a moving Cech 2-cycle (B4, R2, R3).

    W[t] is the witness of triangle T[t] (explicit point or affine spec).  Requires a > b.
    """
    if not a > b:
        return False, "need a > b"
    if not T:
        return False, "empty cycle"
    if not is_mod2_cycle(T):
        return False, "T is not a mod-2 cycle"
    y, ey = resolve_point(m, y_spec)
    ok, why = check_chamber(m, y, a, extra=ey)
    if not ok:
        return False, why
    if len(W) != len(T):
        return False, "one witness per triangle required"
    for t, wspec in zip(T, W):
        p, ep = resolve_point(m, wspec)
        for i in t:
            i = int(i)
            s = m.r[i] + b - m.tol(i) - ep
            if not le_dist(norm2(sub(p, m.x[i])), s):
                return False, f"witness of triangle {tuple(t)} not in shrunken ball {i}"
    par = winding_parity(m.x, T, y, u)
    if par is None:
        return False, "ray direction not in general position"
    if par != 1:
        return False, "cycle does not enclose the chamber point (parity 0)"
    return True, "Cech barrier ok"


# ----------------------------------------------------------------------------- absence certificates

def _box_of(node):
    lo = qv(node["lo"])
    hi = qv(node["hi"])
    if not all(lo[k] < hi[k] for k in range(3)):
        raise ValueError("degenerate box")
    return lo, hi


def _children_ok(lo, hi, kids) -> bool:
    mid = tuple((lo[k] + hi[k]) / 2 for k in range(3))
    want = set()
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                clo = (lo[0] if dx == 0 else mid[0], lo[1] if dy == 0 else mid[1], lo[2] if dz == 0 else mid[2])
                chi = (mid[0] if dx == 0 else hi[0], mid[1] if dy == 0 else hi[1], mid[2] if dz == 0 else hi[2])
                want.add((clo, chi))
    got = set()
    for kd in kids:
        got.add(_box_of(kd))
    return got == want and len(kids) == 8


def check_absence(m: Model, a: Fraction, b: Fraction, tree: dict, o: Vec, bmargin: Fraction) -> tuple[bool, str]:
    """No X in U has a chamber witness y (c_X(y) >= a and E_X(y) <= b), a > b.

    ``tree`` is an octree over a box containing B(o, max_i ||x_i - o|| + eps_i); outside it
    every point has E = c (F1c) and cannot be a witness.  Each leaf carries a test:
      {"test": "T1", "atom": i}               chamber infeasible on the whole leaf (some atom
                                              always within a + r_i of every point)
      {"test": "T2", "id": k, "parent": k'}  leaf robustly clear at level b + bmargin
                                              (bmargin > 0) and joined, through intersecting
                                              T2 leaves with decreasing depth, to a root leaf
                                              containing a point outside the bulk sphere.
      root T2 leaves have "parent": null and must contain a point with ||z-o|| >= R_cover.
    """
    if not a > b:
        return False, "need a > b"
    if bmargin <= 0:
        return False, "bmargin must be positive"
    lo, hi = _box_of(tree)
    # the cover box must contain the ball B(o, rho) with rho >= ||x_i - o|| + eps_i for all i
    hw = min(min(o[k] - lo[k], hi[k] - o[k]) for k in range(3))
    if hw < 0:
        return False, "o not inside cover box"
    for i in range(m.n):
        s = hw - m.eps[i]
        if s < 0 or norm2(sub(m.x[i], o)) > s * s:
            return False, f"cover box does not contain the reach of atom {i}"
    Rcov = hw   # bulk radius: every centre of every X in U lies within hw of o
    leaves = []
    stack = [tree]
    while stack:
        nd = stack.pop()
        lo_, hi_ = _box_of(nd)
        if "children" in nd and nd["children"]:
            if not _children_ok(lo_, hi_, nd["children"]):
                return False, "malformed octree node"
            stack.extend(nd["children"])
        else:
            leaves.append((lo_, hi_, nd))
    t2 = {}
    for lo_, hi_, nd in leaves:
        test = nd.get("test")
        if test == "T1":
            i = int(nd["atom"])
            # every y in the leaf: ||y - x_i|| <= ||y_c - x0_i|| + rhoY + eps_i < a + r0_i - eta_i
            # exact sufficient test: all 8 corners strictly inside B(x0_i, a + r0_i - eps_i - eta_i)
            s = a + m.r[i] - m.tol(i)
            if s <= 0:
                return False, "T1 radius nonpositive"
            for c in _corners(lo_, hi_):
                if norm2(sub(c, m.x[i])) >= s * s:
                    return False, f"T1 fails on leaf {lo_}-{hi_}"
        elif test == "T2":
            # robust clearance >= b + bmargin on the whole leaf: dist(x0_i, leaf) >= r_i + b + bmargin + tol
            for i in range(m.n):
                s = m.r[i] + b + bmargin + m.tol(i)
                if box_dist2(m.x[i], lo_, hi_) < s * s:
                    return False, f"T2 leaf {lo_}-{hi_} not clear of atom {i}"
            t2[nd["id"]] = (lo_, hi_, nd.get("parent"), int(nd.get("depth", -1)))
        else:
            return False, f"leaf without a valid test: {nd.get('test')}"
    # connectivity of T2 leaves to the outside (strictly decreasing depth along parents)
    for k, (lo_, hi_, par, dep) in t2.items():
        if par is None:
            if dep != 0:
                return False, "root T2 leaf must have depth 0"
            far = max(norm2(sub(c, o)) for c in _corners(lo_, hi_))
            if far < Rcov * Rcov:
                return False, "root T2 leaf does not reach the bulk sphere"
        else:
            if par not in t2:
                return False, "dangling parent"
            plo, phi, _, pdep = t2[par]
            if pdep != dep - 1:
                return False, "depth must decrease by one along parents"
            if not all(max(lo_[d], plo[d]) <= min(hi_[d], phi[d]) for d in range(3)):
                return False, "parent leaf does not intersect child leaf"
    return True, "absence ok"


def _corners(lo: Vec, hi: Vec):
    for x in (lo[0], hi[0]):
        for yy in (lo[1], hi[1]):
            for z in (lo[2], hi[2]):
                yield (x, yy, z)


# ----------------------------------------------------------------------------- dispatcher

def verify(cert: dict) -> dict:
    """Verify a certificate. Returns {"verdict", "facts", "reason"}."""
    try:
        m = Model(cert["atoms"], cert.get("eps"), cert.get("eta"))
        claim = cert["claim"]
        kind = claim["type"]
        facts: list[str] = []
        scope = "all X in U(X0;eps,eta)" if (any(m.eps) or any(m.eta)) else "X0"
        if kind == "chamber":
            y, ey = resolve_point(m, claim["y"])
            ok, why = check_chamber(m, y, q(claim["a"]), extra=ey)
            if ok:
                facts.append(f"c_X(y) >= {claim['a']} for {scope}")
        elif kind == "escape_lower":
            y, ey = resolve_point(m, claim["y"])
            if ey != 0:
                return {"verdict": "UNKNOWN", "facts": [], "reason": "path start must be explicit"}
            path = [qv(p) for p in claim["path"]]
            if path[0] != y:
                return {"verdict": "UNKNOWN", "facts": [], "reason": "path must start at y"}
            ok, why = check_path(m, path, q(claim["b"]), qv(claim["bulk"]["o"]), q(claim["bulk"]["R"]))
            if ok:
                facts.append(f"E_X(y) >= {claim['b']} for {scope}")
        elif kind == "escape_upper_box":
            y = qv(claim["y"])
            gr = claim["grid"]
            ok, why = check_box_barrier(m, y, q(claim["b"]), qv(gr["g"]), q(gr["h"]), claim["A"], claim["S"],
                                        qv(claim["bulk"]["o"]), q(claim["bulk"]["R"]))
            if ok:
                facts.append(f"E_X(y) <= {claim['b']} for {scope}")
        elif kind == "pocket":
            a, b = q(claim["a"]), q(claim["b"])
            if not a > b:
                return {"verdict": "UNKNOWN", "facts": [], "reason": "need a > b"}
            bar = claim["barrier"]
            if bar["type"] == "cech":
                ok, why = check_cech_barrier(m, claim["y"], b, a, bar["T"], bar["W"], qv(bar["ray"]))
            elif bar["type"] == "box":
                y, ey = resolve_point(m, claim["y"])
                if ey != 0:
                    return {"verdict": "UNKNOWN", "facts": [], "reason": "box barrier needs explicit y"}
                ok, why = check_chamber(m, y, a)
                if ok:
                    gr = bar["grid"]
                    ok, why = check_box_barrier(m, y, b, qv(gr["g"]), q(gr["h"]), bar["A"], bar["S"],
                                                qv(bar["bulk"]["o"]), q(bar["bulk"]["R"]))
            else:
                return {"verdict": "UNKNOWN", "facts": [], "reason": "unknown barrier type"}
            if ok:
                facts.append(f"P_{{a={claim['a']}, b={claim['b']}}} holds (witness y) for {scope}")
        elif kind == "absence":
            a, b = q(claim["a"]), q(claim["b"])
            ok, why = check_absence(m, a, b, claim["tree"], qv(claim["o"]), q(claim["bmargin"]))
            if ok:
                facts.append(f"no chamber witness for P_{{a={claim['a']}, b={claim['b']}}} for {scope}")
                return {"verdict": "REFUTED", "facts": facts, "reason": why,
                        "note": "the claim P_{a,b} is refuted on the whole region"}
        else:
            return {"verdict": "UNKNOWN", "facts": [], "reason": f"unknown claim type {kind}"}
        return {"verdict": "PROVED" if ok else "UNKNOWN", "facts": facts, "reason": why}
    except Exception as exc:  # malformed certificates never yield a positive verdict
        return {"verdict": "UNKNOWN", "facts": [], "reason": f"malformed certificate: {exc}"}
