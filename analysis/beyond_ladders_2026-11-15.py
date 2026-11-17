"""Question: do the masses have to be read as ladders, or is something else being perceived?
Fact: nine charged masses (and their ratios; inside one charge type the ratios do not change with energy). Picture: the ladder (one charge
type, three families) and the cone (square roots, 45 degrees). Tests of other ways of perceiving the same masses, each with chance controls:
A. any three masses grouped together (by family, by generation, mixed), either sign on the lightest root: which groups sit on a 45-degree cone?
B. a different shape: the cone built on mass^p instead of sqrt(mass); which power makes each charge type 'self = relationship'?
Credit: Koide 1981; Rodejohann & Zhang 2011 (mixed groups); Kartavtsev 2011 (other powers); Xing, Zhang & Zhou 2008/2012; PDG 2024."""
import math, itertools, numpy as np
from scipy.optimize import brentq
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
footA = {"e": 0.51099895, "mu": 105.6583755, "tau": 1776.86, "u": 2.16, "c": 1273.0, "t": 172570.0, "d": 4.70, "s": 93.5, "b": 4183.0}
footB = {"e": 0.48657, "mu": 102.718, "tau": 1746.24, "u": 1.242, "c": 627.8, "t": 162500.0, "d": 2.678, "s": 53.36, "b": 2873.0}   # quarks at the top's energy; leptons at the Z (they barely run between the Z and the top)
def Q(ms, signs): r = np.array(signs)*np.sqrt(ms); return np.sum(ms)/np.sum(r)**2
def scan(m, tol):
    keys = list(m); hits = []
    for tri in itertools.combinations(keys, 3):
        ms = np.array([m[k] for k in tri]); o = np.argsort(ms); ms = ms[o]; tri = tuple(np.array(tri)[o])
        for sg in ((1, 1, 1), (-1, 1, 1)):
            q = Q(ms, sg)
            if abs(q - 2/3) < tol: hits.append((tri, sg[0], q))
    return hits
rng = np.random.default_rng(21)
for name, m in (("A (usual quoted masses)", footA), ("B (all at the top's energy)", footB)):
    head(f"A. Every group of three, footing {name}: which sit on a 45-degree cone (Q = 2/3)?")
    for tol in (1e-3, 1e-4):
        h = scan(m, tol)
        pr(f"  within {tol:g} of 2/3: {len(h)} groups: " + "; ".join(f"{'-' if s < 0 else ''}{'/'.join(t)} ({q:.6f})" for t, s, q in h))
        lo, hi = math.log(min(m.values())), math.log(max(m.values())); c = []
        for _ in range(2000):
            rm = dict(zip(m, np.exp(rng.uniform(lo, hi, 9)))); c.append(len(scan(rm, tol)))
        c = np.array(c); pr(f"      random nine-mass sets: average {c.mean():.2f} such groups; {100*np.mean(c >= len(h)):.0f}% have at least {len(h)}")
    tri_gen = [("e", "u", "d"), ("mu", "c", "s"), ("tau", "t", "b")]
    pr("  Grouped by generation instead of by charge type (no sign flip / lightest flipped):")
    for g in tri_gen:
        ms = np.sort(np.array([m[k] for k in g])); pr(f"    {'/'.join(g):10s} Q = {Q(ms, (1,1,1)):.4f} / {Q(ms, (-1,1,1)):.4f}   (2/3 = 0.6667)")

head("B. A different shape: the 'self = relationship' cone built on mass^p (p = 1/2 is the usual cone)")
def Qp(ms, p): x = np.array(ms)**p; return np.sum(x**2)/np.sum(x)**2
for name, m in (("A", footA), ("B", footB)):
    out = []
    for fam in (("e", "mu", "tau"), ("d", "s", "b"), ("u", "c", "t")):
        ms = [m[k] for k in fam]; p = brentq(lambda p: Qp(ms, p) - 2/3, 0.01, 2); out.append(p)
    pr(f"  footing {name}: power that makes each charge type exact: leptons {out[0]:.4f}, down {out[1]:.4f}, up {out[2]:.4f}")
pr("  Every set of three has such a power (the cone measure runs smoothly from 1/3 to 1 as the power grows), so a power existing says nothing.")
pr("  A pattern would be the powers following the charges or one simple rule. Leptons 1/2 (charge 1), down about 0.40 (charge 1/3), up about 0.31")
pr("  (charge 2/3): not ordered by charge (the down sits between). The quark powers do not depend on energy (ratios do not run), so they are real")
pr("  properties of the masses, but no simple rule links them.")
