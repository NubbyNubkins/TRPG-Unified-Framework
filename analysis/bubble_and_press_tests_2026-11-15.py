"""(A) the bubble picture against the bag model (TRPG's own B) and nuclear data.
(B) First open item: can the down ladder be the up ladder's SAME bubble shape plus only the far-side press?"""
import numpy as np, math, itertools
from scipy.optimize import fsolve, least_squares
hc = 197.327
print("=== (A) bubbles: MIT bag (Chodos, Jaffe, Johnson, Thorn & Weisskopf 1974) with TRPG's B ===")
for Bf in (57, 60, 64):
    B4 = (Bf*hc**3)**0.25
    for N, name in ((3, "proton-like (3 quarks)"), (2, "meson-like (2)")):
        x = 2.04*N
        R = (x/(4*math.pi*Bf/hc))**0.25          # fm, minimising E = 4/3 pi R^3 B + x hc / R
        E = 4/3*math.pi*R**3*Bf + x*hc/R
        print(f"B = {Bf} MeV/fm^3 (B^1/4 = {B4:.1f} MeV): {name}: edge radius {R:.3f} fm, energy {E:.0f} MeV")
print("edge moves with outside pressure P (R ~ (B+P)^-1/4):")
Bf = 60
for P in (0, 6, 30, 60, 180):
    R = (6.12/(4*math.pi*(Bf + P)/hc))**0.25; print(f"  extra outside pressure {P:4d} MeV/fm^3: radius {R:.3f} fm ({100*(R/((6.12/(4*math.pi*Bf/hc))**0.25) - 1):+.1f}%)")
print("liquid-drop surface energy a_s ~ 17.8 MeV * A^(2/3) (Weizsacker 1935): surface tension",
      f"{17.8/(4*math.pi*1.2**2):.2f} MeV/fm^2 (r0 = 1.2 fm)")
print()
print("=== (B) same bubble shape, press on down only ===")
mu = np.array([1.242, 627.8, 162500.]); md = np.array([2.678, 53.36, 2873.])
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857])
J1 = np.ones((3, 3)); SG = [(a, b) for a in (1, -1) for b in (1, -1)]
def fit_up(signs):
    tgt = np.sort(mu*np.array([*signs, 1]))
    build = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]]))
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [mu[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*mu[2]: return build(p)
    return None
def tilts(Mu, Md):
    wu, Ou = np.linalg.eigh(Mu); Ou = Ou[:, np.argsort(np.abs(wu))]
    wd, Od = np.linalg.eigh(Md); Od = Od[:, np.argsort(np.abs(wd))]
    V = Ou.conj().T @ Od; A = np.abs(V); return np.array([A[0, 1], A[1, 2], A[0, 2], A[2, 0]]), np.sort(np.abs(wd))
def score(t): return math.exp(max(abs(np.log(t/Vm))))
ups = [(s, fit_up(s)) for s in SG]; ups = [(s, M) for s, M in ups if M is not None]
print(f"up fits found for {len(ups)} sign choices")
# B1: fold direction only (democratic basis u = (-1,-1,2 e^{i phi})/sqrt6): free c, k, phi -> 3 masses, then 4 tilts predicted
best = None
for s, Mu in ups:
    for phi in np.radians(np.arange(0, 181, 5)):
        u = np.array([-1, -1, 2*np.exp(1j*phi)])/math.sqrt(6); P = np.outer(u, u.conj())
        def res(x):
            c, k = x; w = np.sort(np.abs(np.linalg.eigvalsh(c*Mu - k*162500*P)))
            return np.log(w[:2]/w[2]) - np.log(md[:2]/md[2])
        for k0 in (-0.5, -0.1, 0.05, 0.2, 0.6, 1.2):
            r = least_squares(res, [md[2]/mu[2], k0])
            if np.max(np.abs(r.fun)) < 1e-3:
                t, w = tilts(Mu, r.x[0]*Mu - r.x[1]*162500*P)
                sc = score(t)
                if best is None or sc < best[0]: best = (sc, s, math.degrees(phi), r.x[1], t)
if best: print(f"B1 fold-direction press (3 free: size, strength, timing) reproducing the down masses: best tilts {np.round(best[4], 5)} -> worst {best[0]:.2f}x (signs {best[1]}, timing {best[2]:.0f} deg, strength {best[3]:+.3f})")
else: print("B1: no press along the fold direction turns the up shape into the down masses")
# B2: generous: any complex direction (4 more numbers) -> 6 free for 7 measured
rng = np.random.default_rng(1); best2 = None
for s, Mu in ups:
    for trial in range(400):
        def res2(x):
            c, k, a, b, g, h = x; u = np.array([math.cos(a)*math.cos(b), math.cos(a)*math.sin(b)*np.exp(1j*g), math.sin(a)*np.exp(1j*h)])
            Md = c*Mu - k*162500*np.outer(u, u.conj()); t, w = tilts(Mu, Md)
            return np.concatenate([np.log(w[:2]/w[2]) - np.log(md[:2]/md[2]), np.log(t/Vm)])
        x0 = [md[2]/mu[2], rng.uniform(-1, 1), *rng.uniform(0, math.pi, 4)]
        try: r = least_squares(res2, x0, max_nfev=400)
        except Exception: continue
        sc = math.exp(np.max(np.abs(r.fun)))
        if best2 is None or sc < best2[0]: best2 = (sc, s, r.x)
print(f"B2 any press direction (6 free numbers for 7 measured): worst miss over masses + tilts {best2[0]:.3f}x (signs {best2[1]})")
c, k, a, b, g, h = best2[2]; u = np.array([math.cos(a)*math.cos(b), math.cos(a)*math.sin(b)*np.exp(1j*g), math.sin(a)*np.exp(1j*h)])
print("  best direction |u| components:", np.round(np.abs(u), 3), " strength", round(k, 3))

print()
print("=== (C) the two base shapes: split d and break e for each ladder (no press) ===")
def fit(m, signs):
    tgt = np.sort(m*np.array([*signs, 1])); build = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]]))
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p))) - tgt
    for d0 in (0.0001, 0.001, 0.01, 0.05):
        for e0 in (-0.5, -0.2, 0.01, 0.2, 0.5):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]: return p
    return None
for s in SG:
    pu, pd = fit(mu, s), fit(md, s)
    if pu is None or pd is None: print(f"signs {s}: no fit (up {pu is not None}, down {pd is not None})"); continue
    print(f"signs {s}: up split {pu[1]:.3e} break {pu[2]:+.4f} | down split {pd[1]:.3e} break {pd[2]:+.4f} | ratio split d/u {pd[1]/pu[1]:.2f}, break d/u {pd[2]/pu[2]:.2f}")

print()
print("=== (D) do the 24 and 4.5 depend on the footing? (an earlier test: a quark rule must name its footing) ===")
for name, uu, dd in (("at the top's energy (used above)", mu, md),
                     ("as weighed (2 GeV / own mass, PDG 2024)", np.array([2.16, 1273., 172570.]), np.array([4.70, 93.5, 4183.]))):
    for s in [(1, -1), (-1, 1)]:
        pu, pd = fit(uu, s), fit(dd, s)
        if pu is None or pd is None: print(f"  {name} {s}: no fit"); continue
        print(f"  {name}, signs {s}: split ratio d/u {pd[1]/pu[1]:.2f}, break ratio d/u {pd[2]/pu[2]:.2f}")
