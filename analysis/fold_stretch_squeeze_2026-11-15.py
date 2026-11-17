"""Idea tested: expansion and
shrinkage acting on the quark masses. The fold M -> M - k M0 u u^T IS a squeeze of the shape along one direction u (k > 0) or a stretch (k < 0).
Tests, all along the same direction u = (-1,-1,2)/sqrt6 (an earlier test: the shadow of the ladders' own third-placement break, the TM1 shape):
 V1 down squeezed only; V2 both ladders squeezed by the same amount; V3 down squeezed, up stretched by the same amount;
 V4 each ladder its own amount (2 numbers: fit the 2-3 and 1-2 tilts, predict the 1-3 tilt). Chance control: random directions.
Cause: level repulsion (second-order mixing) with the heavy far-side coloured part: a light state coupled to a heavy one is pushed down by
coupling^2 / heavy mass (Schur 1917; Feshbach 1958; the seesaw idea of Minkowski 1977, Gell-Mann, Ramond & Slansky 1979; for quarks
Barger et al. 1986, Branco & Lavoura 1986). Credit also: Fritzsch & Xing (democratic matrices); Harrison, Perkins & Scott 2002."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
Vm = np.array([0.2250, 0.04183, 0.003732]); J = np.ones((3, 3))
mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
def solve_ladder(m, signs, fold=None):
    tgt = np.sort(np.array(m)*np.array(signs))
    def build(p):
        M = p[0]*(J + np.diag([-p[1], p[1], p[2]]))
        if fold is not None and fold[1] != 0: u, k = fold; M = M - k*p[0]*np.outer(u, u)
        return M
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p)); return O[:, np.argsort(np.abs(w))]
    return None
def turns(U, D): V = np.abs(U.T @ D); return np.array([V[0, 1], V[1, 2], V[0, 2]])
def worst(t): return math.exp(max(abs(np.log(t/Vm))))
SG = list(itertools.product((1, -1), repeat=2))
def one_number(u, mode, ks=np.concatenate([np.linspace(-0.4, -0.005, 80), np.linspace(0.005, 0.4, 80)])):
    """one strength k; mode 'down': fold only down; 'same': both k; 'opposite': down k, up -k. Set the 2-3 tilt exactly; return best."""
    best = None
    for sd in SG:
        for su in SG:
            prev = None
            for k in ks:
                D = solve_ladder(md, [*sd, 1], (u, k)); ku = {"down": 0, "same": k, "opposite": -k}[mode]
                U = solve_ladder(mu, [*su, 1], (u, ku))
                if D is None or U is None: prev = None; continue
                t = turns(U, D)
                if prev is not None and (prev[1][1] - Vm[1])*(t[1] - Vm[1]) < 0 and abs(k - prev[0]) < 0.02:
                    w = (Vm[1] - prev[1][1])/(t[1] - prev[1][1]); tt = prev[1] + w*(t - prev[1]); kk = prev[0] + w*(k - prev[0])
                    if best is None or worst(tt) < worst(best[1]): best = (kk, tt, sd, su)
                prev = (k, t)
    return best
u0 = np.array([-1, -1, 2])/math.sqrt(6)
head("1. What the fold is, geometrically")
pr("  M -> M - k M0 u u^T shrinks the shape along one direction u by the fraction k (k < 0 stretches it), leaving the two directions")
pr("  across u untouched. So the 'expansion and shrinkage' and the fold are one operation: a squeeze or stretch of the hidden matrix")
pr("  along one direction. A turn cannot change the mass ratios; a squeeze can, and it changes the turn between the ladders too.")
head("2. Squeeze and stretch along u = (-1,-1,2), one number")
for mode, label in (("down", "V1 down squeezed only"), ("same", "V2 both squeezed equally"), ("opposite", "V3 down squeezed, up stretched equally")):
    b = one_number(u0, mode)
    if b is None: pr(f"  {label:40s}: no amount gives the measured 2-3 tilt"); continue
    pr(f"  {label:40s}: amount {b[0]:+.4f}  1-2 tilt {b[1][0]:.4f} ({b[1][0]/Vm[0]:.2f}x)  1-3 tilt {b[1][2]:.5f} ({b[1][2]/Vm[2]:.2f}x)  worst {worst(b[1]):.2f}x")
head("3. V4: each ladder its own amount (fit the 2-3 and 1-2 tilts; the 1-3 tilt is the prediction)")
best4 = None
for sd in SG:
    for su in SG:
        def eqs(x):
            D = solve_ladder(md, [*sd, 1], (u0, x[0])); U = solve_ladder(mu, [*su, 1], (u0, x[1]))
            if D is None or U is None: return [9, 9]
            t = turns(U, D); return [math.log(t[1]/Vm[1]), math.log(t[0]/Vm[0])]
        for kd0 in (-0.2, -0.1, 0.1, 0.2):
            for ku0 in (-0.2, -0.05, 0.05, 0.2):
                x, info, ok, _ = fsolve(eqs, [kd0, ku0], full_output=True)
                if ok == 1 and max(abs(np.array(eqs(x)))) < 1e-6:
                    D = solve_ladder(md, [*sd, 1], (u0, x[0])); U = solve_ladder(mu, [*su, 1], (u0, x[1])); t = turns(U, D)
                    if best4 is None or abs(math.log(t[2]/Vm[2])) < abs(math.log(best4[1][2]/Vm[2])): best4 = (x, t, sd, su)
if best4:
    x, t, sd, su = best4
    pr(f"  down amount {x[0]:+.4f}, up amount {x[1]:+.4f}: 2-3 and 1-2 tilts exact; predicted 1-3 tilt {t[2]:.5f} vs {Vm[2]} ({t[2]/Vm[2]:.2f}x)")
    rng = np.random.default_rng(9); N = 60; hits = 0; done = 0
    for i in range(N):
        u = rng.normal(size=3); u /= np.linalg.norm(u); bb = None
        for sd in SG:
            for su in SG:
                def eqs(x):
                    D = solve_ladder(md, [*sd, 1], (u, x[0])); U = solve_ladder(mu, [*su, 1], (u, x[1]))
                    if D is None or U is None: return [9, 9]
                    tt = turns(U, D); return [math.log(tt[1]/Vm[1]), math.log(tt[0]/Vm[0])]
                for kd0, ku0 in ((0.1, 0.05), (-0.1, -0.05), (0.2, -0.1)):
                    xx, info, ok, _ = fsolve(eqs, [kd0, ku0], full_output=True)
                    if ok == 1 and max(abs(np.array(eqs(xx)))) < 1e-6:
                        D = solve_ladder(md, [*sd, 1], (u, xx[0])); U = solve_ladder(mu, [*su, 1], (u, xx[1])); tt = turns(U, D)
                        e = abs(math.log(tt[2]/Vm[2])); bb = e if bb is None else min(bb, e)
        done += bb is not None; hits += (bb is not None and bb <= abs(math.log(t[2]/Vm[2])))
    pr(f"  Control: {N} random directions, same two-number fit: {hits} ({100*hits/N:.0f}%) predict the 1-3 tilt at least as well ({done} had solutions).")
else:
    pr("  No pair of amounts fits both tilts along this direction.")
