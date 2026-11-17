"""Can the bag rules (one shape in three placements, rule 4's one i, the far-side fold on the down side only) produce the
four-zero pattern (Fritzsch & Xing 2003): in the hierarchical basis M = [[0, C, 0], [C*, B~, B], [0, B*, A]] for both up and down?
1. Write the democratic base (J + diag(-d, d, e)) and the fold in the hierarchical basis f1=(1,-1,0)/sqrt2, f2=(1,1,-2)/sqrt6, f3=(1,1,1)/sqrt3.
2. Count free numbers. 3. Test the closest rules-only version: keep only the first-second part of the split (so (1,3) = 0)."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
F = np.array([[1, 1, 1], [-1, 1, 1], [0, -2, 1]], float)/np.array([math.sqrt(2), math.sqrt(6), math.sqrt(3)])
J = np.ones((3, 3)); u = np.array([-1, -1, 2])/math.sqrt(6)
np.set_printoptions(precision=4, suppress=True)
pr("1. The pieces in the hierarchical basis")
pr("  J (one shape, all placements equal):\n", F.T @ J @ F)
pr("  split diag(-1, 1, 0) of the first two placements:\n", F.T @ np.diag([-1, 1, 0.]) @ F)
pr("  third-placement break diag(0, 0, 1):\n", F.T @ np.diag([0, 0, 1.]) @ F)
pr("  fold u u^T along (-1,-1,2):\n", F.T @ np.outer(u, u) @ F)
pr("  => (1,1) is always 0 (Gatto's relation built in); the fold touches only (2,2) (it is four-zero's B~ freedom);")
pr("     but the split gives (1,3) = sqrt2 x (1,2), never 0 unless the light masses vanish. The democratic base is NOT four-zero.")
mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857]); SG = list(itertools.product((1, -1), repeat=2))
E = lambda i, j: np.outer(np.eye(3)[i], np.eye(3)[j])
def build(p, k):   # hierarchical basis: four-zero form with TRPG ratios for the third-placement break
    p0, c, e = p
    return p0*(3*E(2, 2) + c*(E(0, 1) + E(1, 0)) + e*((2/3)*E(1, 1) - (math.sqrt(2)/3)*(E(1, 2) + E(2, 1)) + (1/3)*E(2, 2)) - k*E(1, 1))
def solve(m, signs, k):
    tgt = np.sort(np.array(m)*np.array(signs)); f = lambda p: np.sort(np.linalg.eigvalsh(build(p, k))) - tgt
    for c0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, c0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p, k)); return O[:, np.argsort(np.abs(w))]
    return None
UP = [U for U in (solve(mu, [*s, 1], 0.0) for s in SG) if U is not None]
pr(f"\n2. Free numbers. Quark flavour numbers: 6 masses + 3 tilts + 1 twist = 10 (the four mixing numbers are what is to be explained).")
pr("  Four-zero alone: per ladder 4 sizes (A, B, B~, C), masses fix 3 -> 1 free each; 2 timings -> 4 free for 4 mixing numbers: 0 predictions.")
pr("  With TRPG rules: up has no fold -> its B~ fixed by the break's ratios (0 free); down fold k (1 free); one i fixes timing a = 90 deg;")
pr("  timing b free (1). Total 2 free for 4 mixing numbers -> 2 predictions. The cost: one added reading (the split acts only between the")
pr("  first two placements, not on the third).")
pr(f"  Up ladder solutions found: {len(UP)}")
ks = np.round(np.arange(-0.6, 0.601, 0.02), 3); angs = np.radians(np.arange(0, 360, 15)); res = {}
for k in ks:
    Ds = [D for D in (solve(md, [*s, 1], k) for s in SG) if D is not None]
    for a in angs:
        for b in angs:
            P = np.diag([np.exp(1j*a), np.exp(1j*b), 1.0]); best = (np.inf, None)
            for D in Ds:
                for U in UP:
                    A_ = np.abs(U.T @ P @ D); t = np.array([A_[0,1], A_[1,2], A_[0,2], A_[2,0]]); w = math.exp(max(abs(np.log(t/Vm))))
                    if w < best[0]: best = (w, t)
            res[(k, round(math.degrees(a)), round(math.degrees(b)))] = best
vals = np.array([v[0] for v in res.values()])
pr("\n3. Rules-only four-zero version against the four tilts")
kb = min(res, key=lambda q: res[q][0]); w, t = res[kb]
pr(f"  best overall (k, a, b all free): k {kb[0]:+.2f}, a {kb[1]}, b {kb[2]}: tilts {t[0]:.4f} {t[1]:.5f} {t[2]:.5f} {t[3]:.5f}, worst {w:.3f}x")
for name, sel in (("a = 90 (one i), k and b free", lambda q: q[1] in (90, 270)), ("a = 90, b = 0", lambda q: q[1] in (90, 270) and q[2] == 0),
                  ("a = 90, no fold, b free", lambda q: q[1] in (90, 270) and q[0] == 0)):
    sub = {q: v for q, v in res.items() if sel(q)}; q = min(sub, key=lambda z: sub[z][0]); w, t = sub[q]
    pr(f"  {name:30s}: k {q[0]:+.2f}, a {q[1]}, b {q[2]}: tilts {t[0]:.4f} {t[1]:.5f} {t[2]:.5f} {t[3]:.5f}, worst {w:.3f}x; share of grid as good {100*np.mean(vals <= w*1.003):.2f}%")
pr(f"  share of grid within 1.10x: {100*np.mean(vals <= 1.10):.2f}%; within 1.05x: {100*np.mean(vals <= 1.05):.2f}%")
