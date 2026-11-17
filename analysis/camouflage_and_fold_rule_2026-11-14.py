"""Question: is dark matter camouflaged rather than dark, that is, overlapping particle waves that show as one wave?
Then: tests on the fold.
Credit: Bose 1924; Einstein 1925; Pauli 1925; Hu, Barkana & Gruzinov 2000; Schive, Chiueh & Broadhurst 2014; Hui, Ostriker, Tremaine & Witten 2017;
Rogers & Peiris 2021 (Lyman-alpha bound); Chikashige, Mohapatra & Peccei 1981 (majoron); Berezinsky & Valle 1993; Garcia-Cely & Heeck 2017;
Tremaine & Gunn 1979; Koide 1982; Schur 1917; Feshbach 1958."""
import math, itertools, numpy as np
from scipy.optimize import fsolve, brentq
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
hbar_c = 1.97327e-5   # eV cm
head("1. When do identical particles show as ONE wave? (the picture, made exact)")
pr("  Identical particles that can share a state (bosons: Bose 1924, Einstein 1925) pile into one wave; identical particles that cannot (fermions,")
pr("  Pauli 1925) always stay separate. So 'two identical waves look like one wave' is exact for bosons. They merge when more than one particle")
pr("  sits inside one wavelength-sized box: occupation N = (density / mass) x (wavelength)^3 > 1.")
rho = 0.4e9   # eV / cm^3, dark matter near the Sun
v = 220e5/2.998e10
for m in (1e-22, 1e-20, 1e-6, 1.0, 10.0, 30.0, 1e3, 1e9):
    lam = 2*math.pi*hbar_c/(m*v); N = rho/m*lam**3
    pr(f"   mass {m:8.0e} eV: wavelength {lam:9.2e} cm, particles per wavelength-box {N:9.2e}  -> {'ONE WAVE' if N > 1 else 'separate particles'}")
mc = brentq(lambda m: math.log(rho/m*(2*math.pi*hbar_c/(m*v))**3), 1e-3, 1e4)
pr(f"  Any dark matter BOSON lighter than about {mc:.0f} eV is camouflaged in this sense: its particles overlap into one wave everywhere in a galaxy.")
pr("  This is a real, studied class ('wave dark matter': Hu, Barkana & Gruzinov 2000; Hui et al. 2017). The famous axion is one (1e-6 eV).")
pr("  Checks: (a) it is still dark: overlap hides the particles' individuality, not their mass, so the pull is unchanged; it must also carry no charge,")
pr("  or it would shine. (b) Fermions cannot do this at all, and fermion dark matter must be heavier than ~100 eV to fit into small galaxies")
pr("  (Tremaine & Gunn 1979). (c) The lightest end (1e-22 eV, waves the size of a galaxy's core) is now disfavoured by the Lyman-alpha forest:")
pr("  mass must exceed ~2e-20 eV (Rogers & Peiris 2021). Between ~1e-20 eV and ~30 eV the idea is fully open. HOLDS as a physical possibility.")

head("2. Is there a wave-type (spin 0, chargeless, long-lived) piece already in the main line?")
pr("  an earlier test: the hidden neutrino partner's own mass is the one mass that needs no Higgs; the piece that settles it carries matter number 2.")
pr("  Reading (ii) of an earlier test (matter number a real force above ~1e11 GeV): that piece's phase is eaten by the force's carrier. Nothing left.")
pr("  Reading (i) (the spot is kept exactly, matter number never a force): the same piece's phase is LEFT OVER as a spin-0, chargeless, near-massless")
pr("  wave: the majoron (Chikashige, Mohapatra & Peccei 1981). It is not an extra state; it is part of a piece the main line already has.")
f = 6.3e10*1e9   # eV, settling scale ~ heaviest hidden partner
sum_mnu2 = (0.05)**2 + (0.0086)**2
for mJ in (1e-20, 1e-6, 1.0, 10.0, 1e3):
    Gam = mJ*sum_mnu2/(16*math.pi*f**2)    # eV, decay to neutrinos
    tau = 6.582e-16/Gam
    pr(f"   majoron mass {mJ:7.0e} eV: lifetime ~{tau:9.1e} s (universe 4.4e17 s) -> {'stable enough' if tau > 1e19 else 'too short'}; {'one wave' if mJ < mc else 'separate particles'}")
pr("  So reading (i) gives a camouflaged dark matter candidate that rule 6 allows (an earlier test said rule 6 forbade an extra state; this is not extra).")
pr("  NOT derived: its mass (set by tiny symmetry-breaking effects, e.g. gravity) and how much of it the early universe made. Those decide whether it")
pr("  is ALL the dark matter. Studied in the literature (Berezinsky & Valle 1993; Garcia-Cely & Heeck 2017). Label: CANDIDATE, open.")
pr("  Price: choosing reading (i) gives up the matter-number force (reading ii), which nothing measured needs. That is a choice.")

head("3. A rule for the fold, with nothing fitted (the earlier lead)")
pr("  Rule: the down ladder IS a lepton-like 45-degree cone, disturbed by folding in its heavy charge-1/3 partner.")
pr("  So the fold strength is not fitted: it is whatever makes the UNFOLDED down ladder sit exactly at 45 degrees. Direction: from the cube.")
md = (2.678, 53.36, 2873.0); mu = (1.242, 627.8, 162500.0); Vm = np.array([0.2250, 0.04183, 0.003732]); J = np.ones((3, 3))
def opening(ev):
    x = np.sqrt(np.sort(np.abs(ev))); a = x.mean(); c = x/a - 1; r = math.sqrt((2/3)*(c*c).sum()); return math.degrees(math.atan(r/math.sqrt(2)))
def solve(m, signs, u=None, k=0.0):
    tgt = np.sort(np.array(m)*np.array(signs))
    def build(p):
        M = p[0]*(J + np.diag([-p[1], p[1], p[2]]))
        return M - (k*p[0]*np.outer(u, u) if u is not None else 0), M
    f = lambda p: np.sort(np.linalg.eigvalsh(build(p)[0])) - tgt
    for d0 in (0.001, 0.01, 0.05):
        for e0 in (-0.2, 0.01, 0.2):
            p, info, ok, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                Meff, M0 = build(p); w, O = np.linalg.eigh(Meff); return O[:, np.argsort(np.abs(w))], np.linalg.eigvalsh(M0)
    return None
Us = {su: solve(mu, [*su, 1]) for su in itertools.product((1, -1), repeat=2)}
pr(f"  (down ladder as measured: opening {opening(md):.2f} deg)")
def rule_result(u):
    out = []
    for sd in itertools.product((1, -1), repeat=2):
        ks = np.concatenate([np.linspace(-0.6, -0.005, 120), np.linspace(0.005, 0.6, 120)]); prev = None
        for k in ks:
            r = solve(md, [*sd, 1], u, k)
            if r is None: prev = None; continue
            g = opening(r[1]) - 45.0
            if prev is not None and prev[1]*g < 0:
                k0 = prev[0] + (k - prev[0])*(-prev[1])/(g - prev[1]); r0 = solve(md, [*sd, 1], u, k0)
                if r0 is not None:
                    for su, U in Us.items():
                        if U is None: continue
                        V = np.abs(U[0].T @ r0[0]); out.append((k0, sd, su, np.array([V[0, 1], V[1, 2], V[0, 2]])))
            prev = (k, g)
    return out
dirs = {"face 12+": (1, 1, 0), "body ++-": (1, 1, -1), "axis 3": (0, 0, 1), "body +++": (1, 1, 1), "face 23+": (0, 1, 1), "face 13+": (1, 0, 1),
        "axis 1": (1, 0, 0), "axis 2": (0, 1, 0), "face 12-": (1, -1, 0), "face 23-": (0, 1, -1), "face 13-": (1, 0, -1), "body +-+": (1, -1, 1), "body -++": (-1, 1, 1)}
best_all = []
for name, u in dirs.items():
    u = np.array(u, float); u /= np.linalg.norm(u); res = rule_result(u)
    if not res: pr(f"  {name:9s}: no strength makes the unfolded ladder 45 deg"); continue
    b = min(res, key=lambda r: max(abs(np.log(r[3]/Vm)))); wm = math.exp(max(abs(np.log(b[3]/Vm))))
    best_all.append((wm, name)); pr(f"  {name:9s}: strength {b[0]:+.4f} -> turns {b[3][0]:.4f} {b[3][1]:.4f} {b[3][2]:.5f}   worst miss {wm:.2f}x")
bw = min(best_all)
pr(f"  best cube direction: {bw[1]} (worst turn {bw[0]:.2f}x). No cube direction brings all three turns within 25%, so the random control is not needed:")
pr("  the 45-degree rule for the fold FAILS. The fold still closes the gap when its strength is fitted, but this rule does not fix it.")
