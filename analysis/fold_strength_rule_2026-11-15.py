"""Target: a rule that fixes the fold strength (how strongly the far-side coloured part links to the quarks).
1. The strength in matrix terms: k p0 = g^2 / M_D, so k = (g / sqrt(p0 M_D))^2: the link measured against the balance point of the two scales.
2. Which ladder: the 27 of E6 holds one heavy down-type partner D per family and no heavy up-type partner (Guersey, Ramond & Sikivie 1976;
   Achiman & Stech 1978), so the fold acts on the down ladder only.
3. Candidate rules for the amount, each evaluated with NO fitting, against a chance control (how much of the whole range does as well).
4. What the far-side scale implies for the first-row leak (Cabibbo angle anomaly)."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
Vm = np.array([0.2250, 0.04183, 0.003732]); J = np.ones((3, 3))
mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
u0 = np.array([-1, -1, 2])/math.sqrt(6); P = np.outer(u0, u0)
def solve_ladder(m, signs, k=0.0):
    tgt = np.sort(np.array(m)*np.array(signs))
    build = lambda p: p[0]*(J + np.diag([-p[1], p[1], p[2]])) - k*p[0]*P
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p)); return O[:, np.argsort(np.abs(w))], p
    return None, None
def turns(U, D): V = np.abs(U.T @ D); return np.array([V[0, 1], V[1, 2], V[0, 2]])
def worst(t): return math.exp(max(abs(np.log(t/Vm))))
SG = list(itertools.product((1, -1), repeat=2))
UP = [solve_ladder(mu, [*su, 1])[0] for su in SG]
def at(k):
    best = None
    for sd in SG:
        D, p = solve_ladder(md, [*sd, 1], k)
        if D is None: continue
        for U in UP:
            if U is None: continue
            t = turns(U, D); w = worst(t)
            if best is None or w < best[0]: best = (w, t, p)
    return best
head("1. The fold strength in matrix terms")
p0 = 2873.0/3
pr(f"  Level repulsion: k p0 = g^2 / M_D  ->  k = (g / sqrt(p0 M_D))^2. The ladder scale p0 ~ m_b/3 = {p0:.0f} MeV at the top's energy.")
pr("  sqrt(p0 M_D) is TRPG's balance point between the ladder and the heavy part, so k is the squared link measured against that balance.")
pr(f"  The best folds (earlier tests, down only) have |k| about 0.11-0.12, so g = {math.sqrt(0.116):.2f} x balance (sign = the heavy part's mass sign).")
head("2. Which ladder: down only (E6 27 has a heavy D and no heavy U)")
pr("  27 = 16 + 10 + 1 under SO(10); the 10 holds one down-type colour triplet pair (D, Dbar) and a Higgs-like doublet pair; no up-type partner.")
pr("  So level repulsion from the far-side coloured part can only reach the down ladder: V1 (down only) is the variant the matrix allows.")
head("3. The worst tilt across the whole range of fold strengths (down only)")
ks = np.concatenate([np.linspace(-0.5, -0.0025, 200), [0.0], np.linspace(0.0025, 0.5, 200)])
curve = []
for k in ks:
    b = at(k); curve.append(b[0] if b else np.inf)
curve = np.array(curve)
for k, w in zip(ks[::20], curve[::20]): pr(f"  k {k:+.3f}: worst {w:.2f}x")
i = int(np.argmin(curve)); pr(f"  Best on the grid: k {ks[i]:+.3f}, worst {curve[i]:.2f}x. No fold: worst {curve[200]:.2f}x.")
def chance(w): return np.mean(curve <= w*1.005)
for w in (1.15, 1.20, 1.25): pr(f"  Share of the range [-0.5, 0.5] with worst <= {w:.2f}x: {100*chance(w):.0f}%")
head("4. Candidate rules for the amount (no fitting), with chance control")
g3 = 1.1666; a_s = g3**2/(4*math.pi)
rules = [("1/9 (half the 2/9 twist)", 1/9), ("alpha_s at the top's energy", a_s), ("1/8 (one gluon of eight)", 1/8),
         ("1/6 (u's weight on one light family)", 1/6), ("twist 2/9", 2/9)]
for name, k in rules:
    for s in (-1, 1):
        b = at(s*k)
        if b: pr(f"  {name:38s} k {s*k:+.4f}: tilts {b[1][0]:.4f} {b[1][1]:.5f} {b[1][2]:.5f}  worst {b[0]:.2f}x  chance {100*chance(b[0]):.0f}%")
pr("  Rule A: the fold replaces the ladder's own third-placement break (e = 0): one break, not two.")
bestA = None
for sd in SG:
    tgt = np.sort(np.array(md)*np.array([*sd, 1]))
    f = lambda q: np.sort(np.linalg.eigvalsh(q[0]*(J + np.diag([-q[1], q[1], 0])) - q[2]*q[0]*P)) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for k0 in (-0.4, -0.1, 0.1, 0.4):
            q, info, ok, _ = fsolve(f, [md[2]/3, d0, k0], full_output=True)
            if ok == 1 and np.max(np.abs(f(q))) < 1e-7*md[2]:
                w, O = np.linalg.eigh(q[0]*(J + np.diag([-q[1], q[1], 0])) - q[2]*q[0]*P); D = O[:, np.argsort(np.abs(w))]
                for U in UP:
                    if U is None: continue
                    t = turns(U, D); ww = worst(t)
                    if bestA is None or ww < bestA[0]: bestA = (ww, t, q[2])
if bestA: pr(f"  Rule A: k {bestA[2]:+.4f}: tilts {bestA[1][0]:.4f} {bestA[1][1]:.5f} {bestA[1][2]:.5f}  worst {bestA[0]:.2f}x  chance {100*chance(bestA[0]):.0f}%")
else: pr("  Rule A: no solution (the measured down masses cannot be reached with e = 0 along this direction).")
tries = 2*len(rules) + 1
for w in (1.15, 1.20):
    f = chance(w); pr(f"  Look-elsewhere: {tries} tries; chance at least one lands at <= {w:.2f}x by luck = {100*(1-(1-f)**tries):.0f}%")
head("5. The far-side scale and the first-row leak")
for name, MD in (("far-side pairs", 2.6e7), ("an earlier test exact far-side mass", 2.127e8)):
    th2 = 0.116*p0/(MD*1e3); pr(f"  {name}: M_D = {MD:.2e} GeV -> mixing^2 = k p0 / M_D = {th2:.1e}; first-row leak = mixing^2/6 = {th2/6:.1e}")
pr("  Measured shortfall 1.6e-3. If the fold's heavy part sits at the far-side scale, it leaks only about 1e-10 to 1e-9: the first row must close.")
pr("  Prediction: the 2.2 sigma 'Cabibbo angle anomaly' is not a leak into TRPG's heavy D; it must shrink as the data improve.")
