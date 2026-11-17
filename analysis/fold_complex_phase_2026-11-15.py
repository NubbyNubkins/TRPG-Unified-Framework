"""Follow-up: the real fold gives a flat triangle (J = 0) and |Vtd| 1.6-1.8x too big. Give the squeeze a phase: u = (-1, -1, 2e^{i phi})/sqrt6,
M = p0 (J + diag) - k p0 u u^dagger (Hermitian). Two numbers (k, phi) against four measured: |Vus|, |Vcb|, |Vub|, |Vtd|; J reported."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
J1 = np.ones((3, 3)); mu, md = (1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857])
SG = list(itertools.product((1, -1), repeat=2))
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
            w = math.exp(max(abs(np.log(t/Vm)))); Jc = abs((V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1])).imag)
            if best is None or w < best[0]: best = (w, t, Jc)
    return best
res = []
for k in np.arange(-0.30, 0.301, 0.02):
    if abs(k) < 1e-9: continue
    for phi in np.radians(np.arange(0, 180, 15)):
        b = at(k, phi)
        if b: res.append((b[0], k, math.degrees(phi), b[1], b[2]))
res.sort(key=lambda r: r[0])
pr("Best five (k, phase): worst over Vus, Vcb, Vub, Vtd; J (measured 3.0e-5)")
for w, k, ph, t, Jc in res[:5]:
    pr(f"  k {k:+.2f}, phase {ph:5.1f} deg: Vus {t[0]:.4f} Vcb {t[1]:.5f} Vub {t[2]:.5f} Vtd {t[3]:.5f}  worst {w:.2f}x  J {Jc:.1e}")
pr(f"Real fold (phase 0) best: {min((r for r in res if r[2] == 0), key=lambda r: r[0])[0]:.2f}x")
pr(f"Share of the (k, phase) grid within 1.3x on all four: {100*np.mean([r[0] < 1.3 for r in res]):.0f}%  (grid {len(res)} points)")
