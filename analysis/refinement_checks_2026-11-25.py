"""Refinement checks of 25 November 2026, each fixed before running.
1. Can ANY flavour-blind Z_N (N = 2..18) give the far-side weight rule (D mass from Phi at dimension 4, weak-pair
   mass first at Phi^2) while forbidding Delta B = 1 operators with a lepton and cancelling the anomalies?
2. Does the Higgs self-coupling equal 1/8 (m_H = v/2) at mu = v, or at the wall v/sqrt2? Two-loop running of
   higgs_at_the_grain_2026-11-17.py with PDG 2026 inputs.
3. Do the quark and lepton masses sit on integer loop rungs ln(m_wall/m) / ln(16 pi^2)? Compared with chance.
4. The paper's own pairing-gap band (m_s = 100 MeV, a4 from alpha_s run to 3 mu) with an independent solver
   (cfl_gap_and_tides_2026-11-25.py): heaviest mass and Lambda_1.4 at B = 57.1, 60.1, 63.7 MeV/fm^3.
"""
import itertools, math, os, random, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

print("1. FAR-SIDE WEIGHT RULE FROM A DISCRETE CHARGE, ANY Z_N")
found = []
for N in range(2, 19):
    for qQ, ql, h, dD, dL in itertools.product(range(N), repeat=5):
        qu, qd, qe = (qQ + h) % N, (qQ - h) % N, (ql - h) % N
        if any(o % N == 0 for o in (3 * qQ + ql, 2 * qu + qd + qe, 2 * qQ - qu - qe)): continue
        if dD % N == 0 or not any((dD + s) % N == 0 for s in (1, -1)): continue
        if dL % N == 0 or any((dL + s) % N == 0 for s in (1, -1)) or not any((dL + s) % N == 0 for s in (2, -2)): continue
        if (3 * (2 * qQ - qu - qd) + 3 * dD) % N or (3 * (3 * qQ + ql) + 3 * dL) % N: continue
        found.append((N, qQ, ql, h, dD, dL))
print(f"   solutions found: {len(found)}")
print("   Reason: Yukawa invariance makes the quarks' colour sum vanish, so the colour anomaly is 3 x (D step) = +-3,")
print("   which vanishes mod N only for N = 3; for N = 3, Phi^2 has the charge of Phi*, so the weak pairs get a mass at")
print("   dimension 4 and the rule is lost. The gravitational condition was not even needed.")

print("\n2. HIGGS SELF-COUPLING 1/8 (m_H = v/2)")
src = open(os.path.join(HERE, "higgs_at_the_grain_2026-11-17.py"), encoding="utf8").read()
g = {}; exec(src[:src.index('print("check against')], g)
kw = dict(mt=172.60, a3=0.1180, mw=80.3625)
lam = lambda mu, mh=125.13: g["lam_at"](mu, mh, **kw)[4]
print(f"   lambda(m_t) = {lam(172.57):.5f}, lambda(v/sqrt2 = 174.10) = {lam(174.10):.5f}, lambda(v = 246.22) = {lam(246.22):.5f}")
print(f"   lambda = 1/8 at mu = {brentq(lambda m: lam(m) - 0.125, 173, 2000):.1f} GeV (+-2 GeV from m_H, m_t)")
print(f"   m_H needed for lambda(v) = 1/8: {brentq(lambda m: g['lam_at'](246.22, m, **kw)[4] - 0.125, 115, 135):.2f} GeV vs 125.13 +- 0.11")

print("\n3. LOOP RUNGS ln(174.10 GeV / m) / ln(16 pi^2), MSbar masses at M_Z (MeV)")
m = {"b": 2855.0, "tau": 1746.17, "c": 619.0, "mu": 102.718, "s": 53.16, "d": 2.67, "u": 1.27, "e": 0.48657}
L = math.log(16 * math.pi ** 2); r = {k: math.log(174104.0 / v) / L for k, v in m.items()}
dist = lambda x: abs(x - round(x)); obs = sum(map(dist, r.values())) / len(r)
random.seed(1); lo, hi = min(r.values()), max(r.values()); n = 100000
p = sum(sum(dist(random.uniform(lo, hi)) for _ in r) / len(r) <= obs for _ in range(n)) / n
print("   rungs:", ", ".join(f"{k} {v:.2f}" for k, v in r.items()))
print(f"   mean distance to an integer rung {obs:.3f}; random masses over the same span do as well {p:.0%} of the time")

print("\n4. PAIRING-GAP BAND OF SECTION 9, INDEPENDENT SOLVER")
import importlib
C = importlib.import_module("cfl_gap_and_tides_2026-11-25")
code = os.path.join(HERE, "..", "code", "neutron_stars_no_knobs_2026-11-14.py")
s = open(code, encoding="utf8").read(); ns = {}
exec(s[s.index("def alpha_s"):s.index("def eos")], {"math": math}, ns); alpha_s = ns["alpha_s"]
class RunEOS:
    def __init__(o, B, gap, ms=100.0, X=3.0):
        mu = np.linspace(150, 1300, 6000)
        a4 = np.array([1 - 2 * alpha_s(X * x / 1000) / math.pi for x in mu])
        P = 3 / (4 * math.pi ** 2) * (a4 * mu ** 4 - (ms ** 2 - 4 * gap ** 2) * mu ** 2) / C.HC ** 3 - B
        E = mu * np.gradient(P, mu) - P
        i = int(np.argmax(P > 0)); o.mu0 = np.interp(0, P[i-1:i+1], mu[i-1:i+1]); o.P, o.E = P[i-1:], E[i-1:]
    def eps(o, P): return float(np.interp(max(P, 0), o.P, o.E))
    def deps(o, P, h=0.05): return (o.eps(P + h) - o.eps(max(P - h, 0))) / (P + h - max(P - h, 0))
for B in (57.1, 60.1, 63.7):
    for gap in (82, 97):
        e = RunEOS(B, gap); rows = []
        for Pc in np.geomspace(5, 2500, 60):
            try: rows.append(C.star(e, Pc))
            except Exception: pass
        M = np.array([x[0] for x in rows]); Lm = np.array([x[2] for x in rows]); i = int(M.argmax())
        L14 = float(np.exp(np.interp(1.4, M[:i+1], np.log(Lm[:i+1]))))
        print(f"   B={B} gap={gap} MeV: M_max={M[i]:.3f} Msun, Lambda_1.4={L14:.0f}, E/A={3*e.mu0:.1f} MeV", flush=True)
