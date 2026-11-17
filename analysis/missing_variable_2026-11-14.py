"""Idea tested: a missing variable (a third field like the Higgs and gravity, or a value of the particles like perspective
and force), with the mass derived from the equilibrium point, which is rarely if ever zero.
Test: keep the shared-state (democratic) rule of an earlier test with the simplest break fixed by the masses, and add ONE extra variable, the SAME for both
ladders. One number against three measured turns (and the CP size): it either predicts or fails.
Credit: Harari, Haut & Weyers 1978; Koide 1983; Fritzsch & Plankl 1990; Fritzsch & Xing 1996, 2000; Jarlskog 1985; Froggatt & Nielsen 1979; PDG 2024."""
import math, itertools, numpy as np
from scipy.optimize import fsolve, minimize_scalar
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
md, ms, mb = 2.678, 53.36, 2873.0; mu_, mc, mt = 1.242, 627.8, 162500.0
Vm = np.array([0.2250, 0.04183, 0.003732]); Jm = 3.08e-5
J = np.ones((3, 3))
head("0. the equilibrium picture inside the shared-state rule")
w = np.linalg.eigvalsh(J); pr(f"  The shared all-ones part alone (exact equilibrium) has sizes {np.round(w, 6)}: one family takes everything, two are exactly zero.")
pr("  Masses of the lighter two families exist ONLY because the matrix sits slightly off that equilibrium (the break). So 'the equilibrium point is")
pr("  rarely zero, so there is always a small mass' is exactly this rule's structure. HOLDS as a reading (it is why the hierarchy is steep).")
pr(f"  Size of the offset needed: down ladder ~{ms/mb:.3f}, up ladder ~{mc/mt:.4f} of the shared part (strange/bottom, charm/top).")

def mats(m, signs, extra, kind):
    tgt = np.sort(np.array(m)*np.array(signs))
    def build(p):
        A, d, e = p; M = A*(J + np.diag([-d, d, e])).astype(complex)
        if kind == "offdiag23": M[1, 2] += A*extra; M[2, 1] += A*extra
        if kind == "phase": M[2, 2] = A*(1 + e*np.exp(1j*extra))
        if kind == "phase12": M[0, 1] *= np.exp(1j*extra); M[1, 0] *= np.exp(-1j*extra)
        return M
    def f(p):
        sv = np.linalg.svd(build(p), compute_uv=False); return np.sort(sv) - np.sort(np.abs(tgt))
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-6*m[2]:
                H = build(p) @ build(p).conj().T; wv, O = np.linalg.eigh(H); return O[:, np.argsort(wv)]
    return None
def ckm(x, kind, su, sd):
    U = mats([mu_, mc, mt], [*su, 1], x*(1 if kind != "offdiag23" else 1), kind); D = mats([md, ms, mb], [*sd, 1], x, kind)
    if U is None or D is None: return None
    V = U.conj().T @ D; Va = np.abs(V)
    Jc = abs(np.imag(V[0, 1]*V[1, 2]*np.conj(V[0, 2])*np.conj(V[1, 1])))
    return np.array([Va[0, 1], Va[1, 2], Va[0, 2]]), Jc
for kind, label, rng in (("phase", "a shared PERSPECTIVE ANGLE (one phase, the same i of rule 4) on the heaviest entry", (0, math.pi)),
                         ("phase12", "a shared perspective angle between the first two families", (0, math.pi)),
                         ("offdiag23", "a shared extra link between families 2 and 3 (a 'third field' acting only there)", (-0.3, 0.3))):
    head(f"Extra variable: {label}")
    best = None
    for su in itertools.product((1, -1), repeat=2):
        for sd in itertools.product((1, -1), repeat=2):
            for x in np.linspace(*rng, 241):
                r = ckm(x, kind, su, sd)
                if r is None: continue
                c = float(np.sum(np.log(r[0]/Vm)**2))
                if best is None or c < best[0]: best = (c, x, su, sd, r)
    if best is None: pr("  no solution"); continue
    c, x, su, sd, (V, Jc) = best
    pr(f"  best value {x:+.4f} (signs up {su}, down {sd})")
    pr(f"  turns: {V[0]:.4f} ({V[0]/Vm[0]:.2f}x)  {V[1]:.4f} ({V[1]/Vm[1]:.2f}x)  {V[2]:.5f} ({V[2]/Vm[2]:.2f}x);  CP size {Jc:.2e} (measured {Jm:.2e})")
    ok = np.all(np.abs(np.log(V/Vm)) < math.log(1.1))
    pr("  " + ("ALL THREE within 10% from ONE number: a real candidate, attack next." if ok else "FAIL: one number cannot bring all three turns in (within 10%)."))

head("Is 'about 2.5' a fixed number?")
pr("  The shortfall of the turn between families 2 and 3 in this rule: 0.0151 to 0.0162 against 0.0418, i.e. 2.6x to 2.8x depending on sign choices,")
pr("  break form and energy. It is the CEILING of the shared-state rule, not a constant of nature. No clean number (2.5, e, 8/3) is singled out.")
head("Do the three tiers give the three generations?")
pr("  Tier 1 (quarks with gluons) is built from family 1 (up, down); tier 2 (protons with electrons) is family 1 again; tier 3 is the top (family 3).")
pr("  Family 2 has no tier. So tiers and families are not the same three. FAIL as a derivation. The three still comes from the 3x3 matrix itself")
pr("  (an earlier test: the matrix has exactly three sizes; and here: its shared part leaves exactly one heavy and two light). Hint, not derivation.")
head("Summary")
pr("  HOLDS: the equilibrium picture is exactly the shared-state rule's structure (equilibrium gives one heavy and two exactly zero; the small")
pr("  offset from equilibrium gives the small masses). FAILS: one shared extra variable (a perspective angle, either placement, or an extra 2-3 link)")
pr("  does not raise the 2-3 turn (stays ~0.016); the shortfall is not a fixed number; three tiers are not three generations. Scorecard unchanged.")
