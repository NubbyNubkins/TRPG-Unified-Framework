"""What does the fold of about 0.34 mean? Picture tested: breathing, a rise and fall,
on the back side of the geometry. Testable forms:
 1. In-and-out pair: when the near side is squeezed, does something on the back side expand by the same amount? (full 4x4 matrix with D)
 2. Does the strength change with energy? (needed fold with masses at the Z vs at the top)
 3. A real breathing in time: neutral B mesons swing between matter and antimatter at rates set by |Vtd| and |Vts|. Does the fold predict
    the ratio of the two breathing rates? (Measured |Vtd/Vts| from the Bd and Bs oscillation rates: 0.205 +- 0.006, PDG 2024 with lattice xi.)
 4. Does the geometry itself breathe in time (the strength drifting or cycling)? Limits from Oklo and atomic clocks.
Credit: Wigner & von Neumann 1929; Gell-Mann & Pais 1955; ARGUS 1987 (B mixing); CDF 2006 (Bs oscillation); LHCb; Xing, Zhang & Zhou 2008
(running masses); Shlyakhter 1976; Damour & Dyson 1996; Flambaum & Shuryak 2002; PDG 2024."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
Vm = np.array([0.2250, 0.04183, 0.003732]); J = np.ones((3, 3))
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
def worst(t): return math.exp(max(abs(np.log(t/Vm))))
SG = list(itertools.product((1, -1), repeat=2))
def at(k, mu, md):
    UP = [solve_ladder(mu, [*su, 1])[0] for su in SG]; best = None
    for sd in SG:
        D, p = solve_ladder(md, [*sd, 1], k)
        if D is None: continue
        for U in UP:
            if U is None: continue
            V = np.abs(U.T @ D); t = np.array([V[0, 1], V[1, 2], V[0, 2]]); w = worst(t)
            if best is None or w < best[0]: best = (w, t, V, p)
    return best
mt_set = ((1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0))
mz_set = ((1.27, 619.0, 171700.0), (2.90, 55.0, 2890.0))   # Xing, Zhang & Zhou 2008, at the Z
head("1. In and out: the near side squeezed, the back side pushed out by the same amount")
g2_over_MD = 0.113*mt_set[1][2]/3
M_D = 1.0e6   # any heavy mass (MeV) for the check; the identity does not depend on it
p = next(q for sd in SG for _, q in [solve_ladder(mt_set[1], [*sd, 1], 0.0)] if q is not None)
M3 = p[0]*(J + np.diag([-p[1], p[1], p[2]])); g = math.sqrt(g2_over_MD*M_D)
M4 = np.zeros((4, 4)); M4[:3, :3] = M3; M4[:3, 3] = M4[3, :3] = g*u0; M4[3, 3] = M_D
w3 = np.linalg.eigvalsh(M3); w4 = np.linalg.eigvalsh(M4)
pr(f"  Near side (3 light down states) total shift: {w4[:3].sum() - w3.sum():+.2f} MeV;  back side (heavy D) shift: {w4[3] - M_D:+.2f} MeV")
pr(f"  Sum of the two: {w4.sum() - w3.sum() - M_D:+.2e} MeV (the trace is kept exactly). Whatever the near side gives up along u, the back side takes.")
head("2. Does the strength change with energy? (best fold, masses at the top vs at the Z)")
ks = np.concatenate([np.linspace(-0.3, -0.02, 57), np.linspace(0.02, 0.3, 57)])
for name, (mu, md) in (("at the top", mt_set), ("at the Z (Xing, Zhang & Zhou 2008)", mz_set)):
    res = [(k, at(k, mu, md)) for k in ks]; res = [(k, b) for k, b in res if b]
    kb, bb = min(res, key=lambda x: x[1][0])
    pr(f"  {name:36s}: best k {kb:+.3f}, worst {bb[0]:.2f}x, link/balance {math.sqrt(abs(kb)):.3f}")
pr("  Ladder openings barely run and the tilts barely run, so the needed squeeze is the same within the grid step: no breathing with energy.")
head("3. A real breathing in time: the B-meson matter/antimatter swing")
pr("  Bd and Bs mesons swing between matter and antimatter (rates 0.5065 and 17.765 per picosecond). Their ratio fixes |Vtd/Vts| = 0.205 +- 0.006.")
grid = []
for k in np.linspace(-0.5, 0.5, 201):
    b = at(k, *mt_set)
    if b: grid.append((k, b[2][2, 0]/b[2][2, 1], b[0]))
grid = np.array(grid)
for k in (-0.116, -0.108, 0.108, 0.113):
    b = at(k, *mt_set); r = b[2][2, 0]/b[2][2, 1]
    pr(f"  fold k {k:+.3f}: |Vtd| {b[2][2,0]:.5f}, |Vts| {b[2][2,1]:.5f}, ratio {r:.3f} ({(r-0.205)/0.006:+.1f} sigma)  [worst tilt {b[0]:.2f}x]")
r0 = at(0.0, *mt_set); pr(f"  no fold: ratio {r0[2][2,0]/r0[2][2,1]:.3f}")
good = grid[grid[:, 2] < 1.25]
pr(f"  Among folds with worst tilt under 1.25x, the swing ratio ranges {good[:,1].min():.3f} to {good[:,1].max():.3f}.")
pr(f"  Share of all folds in the range that land within 1 sigma of 0.205: {100*np.mean(abs(grid[:,1]-0.205) < 0.006):.0f}%")
head("4. Does the geometry itself breathe in time?")
pr("  If the squeeze cycled or drifted, quark masses relative to the strong scale would move. Oklo (2 billion years ago) limits that to")
pr("  about 1 part in 10^8 to 10^9 (Flambaum & Shuryak 2002; Damour & Dyson 1996); clocks limit today's drift to ~1e-16 per year or less.")
pr("  So the back side does not breathe in time to that precision: the squeeze is fixed. The breathing that does happen is particles")
pr("  swinging between views (meson oscillation, neutrino oscillation) at rates the fixed squeeze helps set.")
