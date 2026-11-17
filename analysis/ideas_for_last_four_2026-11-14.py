"""Ideas that could close the last four open items, each tested.
Idea 2 first: the fold's direction is the neutrino direction (TM1 pattern), so nothing about the direction is chosen.
Credit: Schur 1917; Feshbach 1958; Harrison, Perkins & Scott 2002 (tribimaximal); Xing & Zhou 2007, Albright & Rodejohann 2009 (TM1);
Fritzsch & Xing (democratic mass matrix); NuFIT 6.0; PDG 2024; Rothstein, Babu & Seckel 1993; Akhmedov, Berezhiani, Mohapatra & Senjanovic 1993
(gravity breaks global symmetries; majoron mass); Kamionkowski & March-Russell 1992, Holman et al. 1992 (the quality problem);
Lazarides & Shafi 1991; Murayama, Suzuki, Yanagida & Yokoyama 1993 (the settling field drives the start)."""
import math, itertools, numpy as np
from scipy.optimize import fsolve, brentq
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
Vm = np.array([0.2250, 0.04183, 0.003732]); J = np.ones((3, 3))
def solve_ladder(m, signs, fold=None):
    tgt = np.sort(np.array(m)*np.array(signs))
    def build(p):
        M = p[0]*(J + np.diag([-p[1], p[1], p[2]]))
        if fold is not None: u, k = fold; M = M - k*p[0]*np.outer(u, u)
        return M
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p)); return O[:, np.argsort(np.abs(w))]
    return None
def turns(U, D): V = np.abs(U.T @ D); return np.array([V[0, 1], V[1, 2], V[0, 2]])
mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)       # MeV at the top's energy (as an earlier test)
Us = {su: solve_ladder(mu, [*su, 1]) for su in itertools.product((1, -1), repeat=2)}
def exact_cb(u):
    """for every sign choice, set the 2-3 turn exactly by the fold strength; return all solutions (k, sd, su, turns)."""
    out = []
    ks = np.concatenate([np.linspace(-0.4, -0.005, 80), np.linspace(0.005, 0.4, 80)])
    for sd in itertools.product((1, -1), repeat=2):
        Ds = [(k, solve_ladder(md, [*sd, 1], (u, k))) for k in ks]
        for su, U in Us.items():
            if U is None: continue
            prev = None
            for k, D in Ds:
                if D is None: prev = None; continue
                t = turns(U, D)
                if prev is not None and (prev[1][1] - Vm[1])*(t[1] - Vm[1]) < 0 and abs(k - prev[0]) < 0.02:
                    w = (Vm[1] - prev[1][1])/(t[1] - prev[1][1]); out.append((prev[0] + w*(k - prev[0]), sd, su, prev[1] + w*(t - prev[1])))
                prev = (k, t)
    return out
def worst(t): return math.exp(max(abs(np.log(t/Vm))))

head("IDEA 2. The fold follows the neutrino direction (no direction chosen; strength set by the 2-3 tilt; the 1-2 and 1-3 tilts predicted)")
s13 = math.sqrt(0.02225); th = math.asin(math.sqrt(3)*s13)
e111 = np.array([1, 1, 1])/math.sqrt(3); e011 = np.array([0, 1, -1])/math.sqrt(2)
cands = {"TM1 first column (2,-1,-1)": np.array([2, -1, -1])/math.sqrt(6), "(-1,2,-1)": np.array([-1, 2, -1])/math.sqrt(6),
         "(-1,-1,2)": np.array([-1, -1, 2])/math.sqrt(6),
         "TM1 second column (measured theta13)": math.cos(th)*e111 + math.sin(th)*e011,
         "TM1 third column (measured theta13)": -math.sin(th)*e111 + math.cos(th)*e011}
pr(f"  TM1 rotation angle from the measured theta13: {math.degrees(th):.2f} degrees. Basis: the same cube basis as an earlier test (democratic axis = (1,1,1)).")
best_all = []
for name, u in cands.items():
    sols = exact_cb(u)
    if not sols: pr(f"  {name:38s}: no strength gives the measured 2-3 tilt"); continue
    b = min(sols, key=lambda s: worst(s[3])); best_all.append((name, b))
    pr(f"  {name:38s}: strength {b[0]:+.4f}  1-2 tilt {b[3][0]:.4f} ({b[3][0]/Vm[0]:.2f}x)  1-3 tilt {b[3][2]:.5f} ({b[3][2]/Vm[2]:.2f}x)  worst {worst(b[3]):.2f}x  [{len(sols)} sign solutions]")
rng = np.random.default_rng(7); N = 200; ws = []
for i in range(N):
    u = rng.normal(size=3); u /= np.linalg.norm(u); sols = exact_cb(u)
    if sols: ws.append(min(worst(s[3]) for s in sols))
    else: ws.append(99)
ws = np.array(ws)
bw = min(worst(b[3]) for _, b in best_all) if best_all else 99
pr(f"  Control: {N} random directions, same procedure: {100*np.mean(ws <= bw):.0f}% do at least as well as the best neutrino direction ({bw:.2f}x);")
pr(f"  {100*np.mean(ws < 1.10):.0f}% get all three tilts within 10%; {100*np.mean(ws < 1.05):.0f}% within 5%.")
