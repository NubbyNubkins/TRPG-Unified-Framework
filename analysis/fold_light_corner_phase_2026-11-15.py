"""Fix the base so the main tilt is right. In the democratic shape the hidden matrices' light corner (in the hierarchical
basis) is automatically empty, so the main tilt is |sqrt(d/s) - e^{i sigma} sqrt(u/c)| (Fritzsch 1977; Fritzsch & Xing 1999):
0.180 if the two light corners are in step, 0.229 if a quarter turn apart, 0.269 if opposite. Test: real fold on the down ladder plus a
relative timing between the two ladders, P = diag(e^{ia}, e^{ib}, 1). Rule candidate: a = 90 deg (rule 4's one i is a quarter turn)."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
J1 = np.ones((3, 3)); mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857]); SG = list(itertools.product((1, -1), repeat=2))
u0 = np.array([-1, -1, 2])/math.sqrt(6); P0 = np.outer(u0, u0)
F = np.array([[1, 1, 1], [-1, 1, 1], [0, -2, 1]], float)/np.array([math.sqrt(2), math.sqrt(6), math.sqrt(3)])
def solve(m, signs, k):
    tgt = np.sort(np.array(m)*np.array(signs))
    build = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]])) - k*p[0]*P0
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(build(p)); return O[:, np.argsort(np.abs(w))]
    return None
pr(f"Main tilt limits: in step {math.sqrt(md[0]/md[1]) - math.sqrt(mu[0]/mu[1]):.4f}, quarter turn {math.hypot(math.sqrt(md[0]/md[1]), math.sqrt(mu[0]/mu[1])):.4f}, opposite {math.sqrt(md[0]/md[1]) + math.sqrt(mu[0]/mu[1]):.4f} (measured 0.2250)")
UP = [U for U in (solve(mu, [*s, 1], 0.0) for s in SG) if U is not None]
ks = np.round(np.arange(-0.30, 0.301, 0.02), 3); angs = np.radians(np.arange(0, 360, 15))
def scan(basis):
    res = {}
    for k in ks:
        Ds = [D for D in (solve(md, [*s, 1], k) for s in SG) if D is not None]
        for a in angs:
            for b in angs:
                Pd = np.diag([np.exp(1j*a), np.exp(1j*b), 1.0]); Pm = Pd if basis == "dem" else F @ Pd @ F.T
                best = np.inf
                for D in Ds:
                    for U in UP:
                        A = np.abs(U.T @ Pm @ D); t = np.array([A[0,1], A[1,2], A[0,2], A[2,0]])
                        best = min(best, math.exp(max(abs(np.log(t/Vm)))))
                res[(k, round(math.degrees(a)), round(math.degrees(b)))] = best
    return res
for basis in ("hier", "dem"):
    res = scan(basis); vals = np.array(list(res.values())); kb = min(res, key=res.get)
    pr(f"\nTiming placed in the {'hierarchical (light-corner)' if basis == 'hier' else 'democratic (placement)'} basis:")
    pr(f"  best overall: k {kb[0]:+.2f}, a {kb[1]} deg, b {kb[2]} deg, worst {res[kb]:.3f}x")
    for name, sel in (("a = 90 (one i), b free, k free", lambda key: key[1] == 90), ("a = 90, b = 0, k free", lambda key: key[1] == 90 and key[2] == 0),
                      ("a = 90, b = 0, no fold", lambda key: key[1] == 90 and key[2] == 0 and key[0] == 0)):
        sub = {kk: v for kk, v in res.items() if sel(kk)}; kk = min(sub, key=sub.get); w = sub[kk]
        pr(f"  {name:32s}: k {kk[0]:+.2f}, b {kk[2]}, worst {w:.3f}x; share of all (k, a, b) as good {100*np.mean(vals <= w*1.003):.1f}%")
    pr(f"  share of grid within 1.10x: {100*np.mean(vals <= 1.10):.1f}%, within 1.05x: {100*np.mean(vals <= 1.05):.2f}%")
