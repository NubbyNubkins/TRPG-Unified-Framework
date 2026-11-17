"""The open problem. (1) A rule for the fold's timing (phase): fixed candidate phases vs the best phase, with chance control.
(2) One heavy part for both jobs? E6 gives one heavy down-type colour triplet D per family (Guersey, Ramond & Sikivie 1976). Does adding
n copies of D (vector-like, (3,1,-1/3)) make the three force strengths meet (one loop, from the top's energy), and at what mass?
Compare with the mass the fold implies if its link equals the Higgs value (g = v)."""
import math, itertools, numpy as np
from scipy.optimize import fsolve, brentq
def pr(*a): print(*a, flush=True)
J1 = np.ones((3, 3)); mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857]); SG = list(itertools.product((1, -1), repeat=2))
def solve(m, signs, k, u):
    tgt = np.sort(np.array(m)*np.array(signs)); Pu = np.outer(u, u.conj())
    build = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]])) - k*p[0]*Pu
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p)); return O[:, np.argsort(np.abs(w))]
    return None
UP = [solve(mu, [*s, 1], 0.0, np.zeros(3)) for s in SG]
def at(k, phi):
    u = np.array([-1, -1, 2*np.exp(1j*phi)])/math.sqrt(6); best = None
    for sd in SG:
        D = solve(md, [*sd, 1], k, u)
        if D is None: continue
        for U in UP:
            if U is None: continue
            V = U.conj().T @ D; A = np.abs(V); t = np.array([A[0,1], A[1,2], A[0,2], A[2,0]])
            w = math.exp(max(abs(np.log(t/Vm))))
            if best is None or w < best[0]: best = (w, t)
    return best
ks = np.arange(0.10, 0.401, 0.01)
def best_k(phi):
    r = [(at(k, phi), k) for k in ks]; r = [(b[0], k, b[1]) for b, k in r if b]
    return min(r) if r else (np.inf, None, None)
pr("1. The fold's timing: best squeeze at each phase (worst of the four tilts)")
phis = np.arange(0, 181, 5); curve = []
for ph in phis:
    w, k, t = best_k(math.radians(ph)); curve.append(w)
    if ph % 15 == 0: pr(f"  phase {ph:3d} deg: best k {k}, worst {w:.2f}x")
curve = np.array(curve); i = int(np.argmin(curve)); pr(f"  Best phase on the 5-degree grid: {phis[i]} deg, worst {curve[i]:.2f}x")
for name, ph in (("60 deg (sixth of a turn)", 60), ("90 deg (quarter turn)", 90), ("120 deg (third of a turn, three placements)", 120),
                 ("2/9 rad = 12.7 deg (the lepton twist)", math.degrees(2/9))):
    w, k, t = best_k(math.radians(ph)); share = np.mean(curve <= w*1.005)
    pr(f"  {name:44s}: k {k:.2f}, worst {w:.2f}x, share of phases as good {100*share:.0f}%")
pr("\n2. One heavy part for both jobs: D quarks and the three strengths (one loop from the top's energy)")
g0 = np.array([0.4626, 0.6478, 1.1666]); a0 = 4*math.pi/g0**2; mt = 172.57
b = np.array([41/10, -19/6, -7.0]); dbD = np.array([4/15, 0.0, 2/3])
def inv(mu_, n, MD):
    L = math.log(mu_/mt); x = a0 - b/(2*math.pi)*L
    if n and mu_ > MD: x = x - n*dbD/(2*math.pi)*math.log(mu_/MD)
    return x
def mismatch(n, MD):
    f12 = lambda lm: (lambda x: x[0] - x[1])(inv(math.exp(lm), n, MD))
    lm = brentq(f12, math.log(1e3), math.log(1e19)); x = inv(math.exp(lm), n, MD)
    return x[2] - x[0], math.exp(lm)
d0, M0 = mismatch(0, 1e30); pr(f"  Standard Model alone: 1/alpha3 - 1/alpha12 = {d0:+.2f} at {M0:.1e} GeV (zero means the three meet)")
for n in (1, 2, 3):
    f = lambda lM: mismatch(n, math.exp(lM))[0]
    try:
        lM = brentq(f, math.log(200), math.log(1e15)); MD = math.exp(lM); _, MG = mismatch(n, MD)
        pr(f"  {n} D quark(s): the strengths meet exactly if D weighs {MD:.2e} GeV; they meet at {MG:.2e} GeV")
    except ValueError: pr(f"  {n} D quark(s): no mass makes them meet")
v = 246.22; p0 = 2.873/3
for name, k in (("real fold k = 0.12", 0.12), ("phased fold k = 0.26", 0.26)):
    pr(f"  Fold with link = Higgs value v: M_D = v^2/(k p0) = {v**2/(k*p0):.2e} GeV ({name})")
pr("  Proton decay needs the meeting point above about 5e15 GeV (Super-Kamiokande); LHC needs D above about 1.5 TeV.")
