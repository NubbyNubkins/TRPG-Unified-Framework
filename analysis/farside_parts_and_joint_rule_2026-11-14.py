"""Tested: the six coloured far-side parts as real particles (strengths), then one joint rule for the
up and down quark ladders (masses + tilt).
Credit: Georgi, Quinn & Weinberg 1974; Georgi & Glashow 1974 (the coloured partner of the Higgs pair); Dimopoulos, Raby & Wilczek 1981; Amaldi, de Boer
& Furstenau 1991; Harari, Haut & Weyers 1978 (democratic matrices); Koide 1983, 1989; Fritzsch & Plankl 1990; Fritzsch & Xing 1996, 2000;
Super-Kamiokande 2020 (proton lifetime); PDG 2024; FLAG 2024."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
mZ = 91.1876; ainv = 127.951; a3inv = 1/0.1180; s2m = 0.23122
BSM = np.array([41/10, -19/6, -7.0])
def predict(extra=()):
    """one loop; extra = list of (mass, db vector). Strong = weak at M_U, hypercharge must meet there too. Returns (sin^2 at Z, M_U, mismatch)."""
    def run(a0, L, i):
        tot = BSM[i]*L
        for M, db in extra:
            lM = math.log(M/mZ)
            if lM < L: tot += db[i]*(L - lM)
        return a0 - tot/(2*math.pi)
    def resid(a2):
        lo, hi = 0.0, 80.0
        for _ in range(200):
            mid = (lo + hi)/2; f = run(a2, mid, 1) - run(a3inv, mid, 2)
            fl = run(a2, lo, 1) - run(a3inv, lo, 2); (lo, hi) = (mid, hi) if fl*f > 0 else (lo, mid)
        L = (lo + hi)/2; return run((3/5)*(ainv - a2), L, 0) - run(a2, L, 1), L
    lo, hi = 20.0, 40.0
    for _ in range(100):
        mid = (lo + hi)/2; (lo, hi) = (mid, hi) if resid(lo)[0]*resid(mid)[0] > 0 else (lo, mid)
    a2 = (lo + hi)/2; return a2/ainv, mZ*math.exp(resid(a2)[1])

head("PART A. The six coloured far-side parts as real particles")
pr("  What they are (44.2): with the colourless weak pair (the Higgs) they fill the far side's 10 parts. Six = a coloured triplet of charge 1/3,")
pr("  spin 0 (three colours x particle/antiparticle). This is exactly the Higgs pair's coloured partner in the 1974 unified models.")
db_T = np.array([(3/5)*(1/3)*3*(1/3)**2, 0.0, (1/3)*(1/2)])
pr(f"  Its effect on the running: hypercharge +{db_T[0]:.4f}, weak 0, strong +{db_T[2]:.4f}")
s0, M0 = predict(); pr(f"  without it: sin^2 = {s0:.4f} ({100*(s0/s2m - 1):+.1f}%), meeting at {M0:.1e} GeV")
for M in (1e3, 1e6, 1e10, 1e13):
    s, MU = predict([(M, db_T)]); pr(f"  with it from {M:7.0e} GeV: sin^2 = {s:.4f} ({100*(s/s2m - 1):+.1f}%), meeting at {MU:.1e} GeV")
pr("  It makes the gap WIDER at every mass. Reason: it slows the strong force's weakening (strong +1/6) and touches")
pr("  the weak force not at all, so the weak and strong lines meet later, where hypercharge has already passed. Same for ANY number of copies.")
pr("  ATTACK from the other side (could they be spin 1/2 instead?): a coloured triplet of spin 1/2 adds hypercharge 2/15 and strong 2/3:")
s, MU = predict([(1e3, 2*db_T*np.array([1, 0, 2]) )]); pr(f"    from 1 TeV: sin^2 = {s:.4f} ({100*(s/s2m - 1):+.1f}%). Also wider.")
pr("  And they cannot be light in any case: a coloured charge-1/3 piece that touches quarks and leptons makes the proton decay. Super-Kamiokande")
pr("  (proton lasts > 2.4e34 years) puts it above ~1e11 GeV even with the smallest (electron-sized) couplings. So as real particles: FAIL twice.")
pr("  What closes the gap instead, for the record: things carrying the WEAK charge between the Z and the meeting. Required change: b_weak - b_strong")
need = []
for M in (1e3,):
    for k in np.linspace(0, 6, 601):
        s, MU = predict([(M, np.array([0.0, k, 0.0]))])
        if abs(s - s2m) < 5e-4: need.append(k)
pr(f"  must grow by about {np.mean(need):.2f} from 1 TeV (pure weak-charged additions). The far side has no such part beyond the Higgs pair. FAIL.")

head("PART B. One joint rule for the up and down ladders, tilt = angle between them")
md, ms, mb = 2.678, 53.36, 2873.0; mu_, mc, mt = 1.242, 627.8, 162500.0       # one energy (the top's), MeV
Vm = {'us': 0.2250, 'cb': 0.04183, 'ub': 0.003732}
pr("  B1. The shared-state rule. TRPG reads each ladder as a cone around one axis, and the two ladders as hung from one shared state (47.4).")
pr("      Written as matrices: each ladder = (all-ones matrix: the shared state, every family touching every family equally) + a small break.")
pr("      This is the 'democratic' form (Harari, Haut & Weyers 1978; Koide; Fritzsch & Plankl 1990). Both ladders share the SAME all-ones part;")
pr("      they differ only in the break. The tilt then comes out; nothing is fitted to it. Break used: (-d, +d, e) on the diagonal (the simplest that")
pr("      gives three different masses; Fritzsch & Xing 1996). Three masses fix three numbers per ladder.")
J = np.ones((3, 3))
def ladder(m, signs):
    tgt = np.sort(np.array(m)*np.array(signs))
    def f(p):
        A, d, e = p; ev = np.linalg.eigvalsh(A*(J + np.diag([-d, d, e]))); return ev - tgt
    best = None
    for A0 in (m[2]/3,):
        for d0 in (0.001, 0.01, 0.05):
            for e0 in (-0.1, 0.01, 0.2):
                p, info, ok, _ = fsolve(f, [A0, d0, e0], full_output=True)
                if ok == 1 and np.max(np.abs(f(p))) < 1e-6*m[2]:
                    best = p; break
            if best is not None: break
        if best is not None: break
    if best is None: return None
    w, O = np.linalg.eigh(best[0]*(J + np.diag([-best[1], best[1], best[2]])))
    order = np.argsort(np.abs(w)); return O[:, order], best
res = []
for su in itertools.product((1, -1), repeat=2):
    for sd in itertools.product((1, -1), repeat=2):
        U = ladder([mu_, mc, mt], [*su, 1]); D = ladder([md, ms, mb], [*sd, 1])
        if U is None or D is None: continue
        V = np.abs(U[0].T @ D[0]); res.append((su, sd, V[0, 1], V[1, 2], V[0, 2]))
pr(f"      {len(res)} sign choices solve. Results (|V_us|, |V_cb|, |V_ub|), measured ({Vm['us']}, {Vm['cb']}, {Vm['ub']}):")
for su, sd, us, cb, ub in res: pr(f"        up signs {su}, down signs {sd}:  {us:.4f}  {cb:.4f}  {ub:.5f}")
bestr = min(res, key=lambda r: sum(math.log(r[i]/Vm[k])**2 for i, k in zip((2, 3, 4), ('us', 'cb', 'ub'))))
pr(f"      best: {bestr[2]:.4f} ({bestr[2]/Vm['us']:.2f}x), {bestr[3]:.4f} ({bestr[3]/Vm['cb']:.2f}x), {bestr[4]:.5f} ({bestr[4]/Vm['ub']:.2f}x)")
pr("  B2. The tilt as the plain angle between the two ladders (each ladder as the arrow (sqrt m1, sqrt m2, sqrt m3)):")
a = np.sqrt([mu_, mc, mt]); b = np.sqrt([md, ms, mb]); ang = math.degrees(math.acos(a @ b/np.linalg.norm(a)/np.linalg.norm(b)))
pr(f"      angle = {ang:.2f} degrees; Cabibbo angle 13.0, V_cb angle 2.40, V_ub angle 0.21. Matches none (and it uses the masses' own sizes).")
pr("  B3. The cone form of the same idea: the two cones' axes tilted by the measured tilt. an earlier test C5 showed the same-axis start gives zero tilt;")
pr("      putting the measured tilt in by hand gives back the measured tilt and no masses: it relates, it does not produce. Not a rule.")

head("Summary")
pr("  A: the six coloured far-side parts as real particles widen the strength gap and would break the proton. FAIL. The strengths need weak-charged")
pr("     states between the Z and the meeting energy, which the main line does not have. Force groups stays PARTIAL.")
pr("  B: the shared-state (democratic) joint rule is the right SHAPE (both ladders from one shared part, the tilt comes out unfitted) but its")
pr("     numbers fail; the plain angle fails. Quark masses stay PARTIAL.")

head("ATTACK on B1: give the shared-state rule its most general diagonal break (one free number per ladder) and let it FIT the tilt")
def ladder_x(m, signs, x):
    tgt = np.sort(np.array(m)*np.array(signs))
    f = lambda p: np.linalg.eigvalsh(p[0]*(J + np.diag([x - p[1], x + p[1], p[2]]))) - tgt
    for d0 in (0.001, 0.01, 0.05, 0.2):
        for e0 in (-0.5, -0.1, 0.01, 0.2, 0.5):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-6*m[2]:
                w, O = np.linalg.eigh(p[0]*(J + np.diag([x - p[1], x + p[1], p[2]]))); return O[:, np.argsort(np.abs(w))]
    return None
xs = np.concatenate([np.linspace(-0.3, -0.03, 10), np.linspace(-0.02, 0.02, 41), np.linspace(0.03, 0.3, 10)]); best = (1e9,)
for su in itertools.product((1, -1), repeat=2):
    for sd in itertools.product((1, -1), repeat=2):
        Us = {x: ladder_x([mu_, mc, mt], [*su, 1], x) for x in xs}; Ds = {x: ladder_x([md, ms, mb], [*sd, 1], x) for x in xs}
        for xu, U in Us.items():
            if U is None: continue
            for xd, D in Ds.items():
                if D is None: continue
                V = np.abs(U.T @ D); c = sum(math.log(V[i, j]/Vm[k])**2 for (i, j), k in (((0, 1), 'us'), ((1, 2), 'cb'), ((0, 2), 'ub')))
                if c < best[0]: best = (c, xu, xd, V[0, 1], V[1, 2], V[0, 2])
pr(f"  best over 2 free numbers (one per ladder) and all signs: |V_us| {best[3]:.4f} ({best[3]/Vm['us']:.2f}x), |V_cb| {best[4]:.4f} ({best[4]/Vm['cb']:.2f}x),")
pr(f"  |V_ub| {best[5]:.5f} ({best[5]/Vm['ub']:.2f}x)   (break shifts {best[1]:+.3f}, {best[2]:+.3f})")
pr("  Verdict: even with a free number per ladder (break shifts from -0.3 to +0.3) the turn to the third family stays 0.39x too small. The first and")
pr("  outer turns can be fitted, the middle one cannot. The shared-state rule with a diagonal break is RULED OUT, not only unfitted.")
