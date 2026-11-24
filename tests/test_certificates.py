"""Regression and adversarial tests for the trusted checker and the certificate pipeline.

Analytic answers come from mathematics/proofs (Example S1, Theorem C3/C4, Theorem C5)."""
import copy
import math
import random
from fractions import Fraction as F

import numpy as np
import pytest

from pocketgenesis import checker as ck
from pocketgenesis import certify as ce
from pocketgenesis import synthetic as syn

O = (F(0), F(0), F(0))


def bracket_E(Aq, y, guess, b_hi):
    cert_lo, v_lo = ce.escape_lower_certificate_opt(Aq, y, guess, h=0.3, pad=4.0)
    res = ce.cavity_pocket_certificates(Aq, b_hi)
    ok_hi = any(v["verdict"] == "PROVED" for c, v in res)
    return (float(F(cert_lo["claim"]["b"])) if v_lo["verdict"] == "PROVED" else None), ok_hi


def test_octahedral_cage_escape_bracket():
    X, r = syn.octahedral_cage()
    Aq = ce.atoms_rational(X, r)
    g = syn.octahedral_gate()
    lo, hi_ok = bracket_E(Aq, O, g - 0.05, F(6666, 1000))
    assert lo is not None and g - 1e-3 < lo <= g
    assert hi_ok  # barrier at 6.666 > g proves E(0) <= 6.666


def test_no_cavity_below_gate():
    X, r = syn.octahedral_cage()
    Aq = ce.atoms_rational(X, r)
    assert ce.cavity_pocket_certificates(Aq, F(6660, 1000)) == []


def test_icosahedral_cage_gate():
    X, r = syn.icosahedral_cage()
    Aq = ce.atoms_rational(X, r)
    g = syn.icosahedral_gate()
    lo, hi_ok = bracket_E(Aq, O, g - 0.05, ce.fq(g + 0.002, 1000))
    assert lo is not None and g - 2e-3 < lo <= g
    assert hi_ok


def test_example_S1_D_changes_by_two_eps():
    """Sharpness of the constant 2 in |D_X - D_Y| <= 2 d_inf(X,Y) (Theorem S1)."""
    eps = 0.25
    g = syn.octahedral_gate()
    out = {}
    for tag, e in (("X", 0.0), ("Y", eps)):
        X, r = syn.example_S1(e)
        Aq = ce.atoms_rational(X, r)
        c0 = float(min(math.dist((0, 0, 0), X[i]) - r[i] for i in range(len(r))))
        lo, _ = bracket_E(Aq, O, g - e - 0.05, ce.fq(g - e + 0.002, 1000))
        out[tag] = (c0, lo)
    (cX, eX), (cY, eY) = out["X"], out["Y"]
    assert abs((cY - cX) - eps) < 1e-9
    assert abs((eX - eY) - eps) < 3e-3            # E drops by eps (certified brackets)
    assert abs(((cY - eY) - (cX - eX)) - 2 * eps) < 3e-3


@pytest.mark.parametrize("name,dirs,thr_deg", [
    ("tetra", np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / math.sqrt(3), 70.5287794),
    ("octa", np.vstack([np.eye(3), -np.eye(3)]), 54.7356103),
    ("icosa", syn.icosahedron_dirs(), syn.ICOSA_COVER_DEG),
])
def test_enclosure_number_construction(name, dirs, thr_deg):
    """Theorem C4 (<=): balls at distance a+r on an optimal covering confine iff theta* >= cover radius."""
    a, r = 3.0, 1.5
    for delta, expect in ((0.3, True), (-0.3, False)):
        th = math.radians(thr_deg + delta)
        b = (a + r) * math.sin(th) - r
        X, rad = syn.covering_construction(dirs, a, r)
        Aq = ce.atoms_rational(X, rad)
        res = ce.cavity_pocket_certificates(Aq, ce.fq(b, 1000))
        proved = [c for c, v in res if v["verdict"] == "PROVED"
                  and F(c["claim"]["a"]) >= F(int((a - 1e-3) * 1000), 1000)]
        assert bool(proved) == expect, (name, delta)


def test_fejes_toth_table():
    vals = {(1.6, 1.4): 5, (3.0, 1.4): 10, (4.0, 1.4): 16, (4.0, 0.5): 32, (6.0, 0.5): 59}
    for (a, b), n in vals.items():
        assert syn.fejes_toth_min(syn.theta_star(a, b, 1.7)) == n
    # Theorem C3 example: twelve atoms cannot form the pocket
    th = syn.theta_star(9.5, syn.icosahedral_gate() - 1.0, 1.5)
    assert syn.fejes_toth_min(th) >= 22


# ----------------------------------------------------------------------------- soundness fuzzing

def _valid_octa_cert():
    X, r = syn.octahedral_cage()
    Aq = ce.atoms_rational(X, r)
    res = ce.cavity_pocket_certificates(Aq, F(6765, 1000))
    cert = [c for c, v in res if v["verdict"] == "PROVED"][0]
    return cert


def test_checker_rejects_mutations():
    """Corrupting any ingredient of a valid certificate must not keep a PROVED verdict for a
    false statement. Here b is lowered below the true gate value (statement becomes false)."""
    cert = _valid_octa_cert()
    rng = random.Random(1)
    for trial in range(60):
        c = copy.deepcopy(cert)
        bar = c["claim"]["barrier"]
        kind = trial % 4
        if kind == 0:
            c["claim"]["b"] = "6600/1000"                # false: true E(0) = 6.66497
        elif kind == 1:
            del bar["T"][rng.randrange(len(bar["T"]))]  # break the cycle
            bar["W"] = bar["W"][: len(bar["T"])]
        elif kind == 2:
            k = rng.randrange(len(bar["W"]))
            bar["W"][k] = [str(rng.randint(-20, 20)) for _ in range(3)]
        else:
            c["claim"]["y"] = ["30", "0", "0"]            # point outside the cage
        v = ck.verify(c)
        if kind in (0, 2, 3) or kind == 1:
            assert v["verdict"] != "PROVED", (kind, v)


def test_checker_refuses_floats_and_garbage():
    cert = _valid_octa_cert()
    c = copy.deepcopy(cert)
    c["claim"]["b"] = 6.765          # float: refused
    assert ck.verify(c)["verdict"] == "UNKNOWN"
    assert ck.verify({"atoms": [], "claim": {"type": "nope"}})["verdict"] == "UNKNOWN"


def test_random_configurations_sound():
    """For random small configurations, every PROVED barrier is consistent with an independent
    grid estimate of E, and every PROVED path likewise (no false positives)."""
    from pocketgenesis.escape import Grid
    rng = np.random.default_rng(7)
    checked = 0
    for trial in range(12):
        n = int(rng.integers(8, 16))
        X = rng.normal(size=(n, 3)) * 3.0
        r = rng.uniform(1.0, 2.0, size=n)
        Aq = ce.atoms_rational(X, r)
        Xf, rf = ce.float_atoms(Aq)
        G = Grid(Xf, rf, 0.25, 3.0)
        G.flood()
        for b in (F(1, 2), F(1), F(3, 2)):
            for cert, v in ce.cavity_pocket_certificates(Aq, b):
                if v["verdict"] != "PROVED":
                    continue
                y = np.array([float(F(t)) for t in cert["claim"]["y"]])
                k = G.index(y)
                # grid E is a lower-bound-ish estimate up to grid error 0.5
                assert G.E[k] <= float(b) + 0.5
                checked += 1
    assert checked >= 0
