"""Tests run during the full review of everything (the original statements and the record).
T1. the G2 combining rule (band-stepping = Schur/Feshbach folding, M_eff = M_LL - M_LH M_HH^-1 M_HL): fold ONE heavy state into the down ladder of the shared-state rule. The far side's coloured parts carry charge 1/3, so they can
    only touch the down ladder. The fold direction is not fitted: it must be one of the natural directions of the 'cube of growing cones'
    (3 axes, 6 face diagonals, 4 body diagonals). One number (the fold's strength) against three measured turns.
T2. The 'missing variable' (relativity, the observer) read as the energy at which masses are compared: does the shared-state
    ceiling move with energy? (an earlier test said 2.6-2.8x 'depending on energy' without a run; this checks it.)
T3. the bands (J1) pair strange with muon and bottom with tau: the Georgi-Jarlskog rule (1979), energy-free form: strange/down = (muon/electron)/9.
Credit: Schur 1917; Feshbach 1958; Harari, Haut & Weyers 1978; Fritzsch & Plankl 1990; Georgi & Jarlskog 1979; Gursey, Ramond & Sikivie 1976 (E6);
Barger et al. 1986 (exotic down quarks); FLAG 2024; PDG 2024."""
import math, itertools, numpy as np
from scipy.optimize import fsolve
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
sets = {"top's energy (162.5 GeV)": ((1.242, 627.8, 162500.0), (2.678, 53.36, 2873.0)),
        "Z energy (91 GeV), alternative inputs": ((7.04e-6*174104, 3.56e-3*174104, 0.967*174104), (1.54e-5*174104, 3.06e-4*174104, 1.630e-2*174104))}
head("T2. Does the shared-state ceiling move with energy?")
for name, (mu, md) in sets.items():
    best = None
    for su in itertools.product((1, -1), repeat=2):
        for sd in itertools.product((1, -1), repeat=2):
            U = solve_ladder(mu, [*su, 1]); D = solve_ladder(md, [*sd, 1])
            if U is None or D is None: continue
            t = turns(U, D); c = float(np.sum(np.log(t/Vm)**2))
            if best is None or c < best[0]: best = (c, t)
    pr(f"  {name}: turns {np.round(best[1], 5)}  -> 2-3 turn short by {Vm[1]/best[1][1]:.2f}x")
pr("  The ceiling barely moves with energy (mass ratios change only a little between these energies). The energy is not the missing variable.")

head("T1. the band-stepping (one heavy state folded into the down ladder), fold direction from the cube of cones")
mu, md = sets["top's energy (162.5 GeV)"]
dirs = {"axis 1": (1, 0, 0), "axis 2": (0, 1, 0), "axis 3": (0, 0, 1), "face 23+": (0, 1, 1), "face 23-": (0, 1, -1), "face 13+": (1, 0, 1),
        "face 13-": (1, 0, -1), "face 12+": (1, 1, 0), "face 12-": (1, -1, 0), "body +++": (1, 1, 1), "body ++-": (1, 1, -1), "body +-+": (1, -1, 1),
        "body -++": (-1, 1, 1)}
results = []
Us = {su: solve_ladder(mu, [*su, 1]) for su in itertools.product((1, -1), repeat=2)}
for name, u in dirs.items():
    u = np.array(u, float); u /= np.linalg.norm(u); best = None
    for sd in itertools.product((1, -1), repeat=2):
        for k in np.concatenate([np.linspace(-0.6, -0.01, 60), np.linspace(0.01, 0.6, 60)]):
            D = solve_ladder(md, [*sd, 1], (u, k))
            if D is None: continue
            for su, U in Us.items():
                if U is None: continue
                t = turns(U, D); c = float(np.sum(np.log(t/Vm)**2))
                if best is None or c < best[0]: best = (c, k, t)
    if best: results.append((name, *best))
results.sort(key=lambda r: r[1])
pr(f"  {'direction':10s} {'strength':>9s}   turns (measured 0.2250, 0.04183, 0.003732)          worst miss")
for name, c, k, t in results:
    pr(f"  {name:10s} {k:+9.3f}   {t[0]:.4f} ({t[0]/Vm[0]:.2f}x)  {t[1]:.4f} ({t[1]/Vm[1]:.2f}x)  {t[2]:.5f} ({t[2]/Vm[2]:.2f}x)   {math.exp(max(abs(np.log(t/Vm)))):.2f}x")
b = results[0]
ok = math.exp(max(abs(np.log(b[3]/Vm)))) < 1.10
pr("  " + (f"CANDIDATE: '{b[0]}' brings all three turns within 10% with one number. Attack next." if ok else
           f"Best is '{b[0]}' (worst turn off by {math.exp(max(abs(np.log(b[3]/Vm)))):.2f}x). FAIL as a one-number fix."))

head("T3. the bands pair strange with muon, bottom with tau (Georgi & Jarlskog 1979)")
r_pred = (105.6583755/0.51099895)/9; r_meas = 20.0; pr(f"  strange/down predicted (muon/electron)/9 = {r_pred:.2f}; measured 20.0 +- 0.4 (FLAG 2024): {100*(r_pred/r_meas - 1):+.0f}%, about 7 sigma. FAIL.")
pr("  (The pairing itself is real: in the bands and in many unified models. The factor 3 rule is what fails.)")

head("T1 ATTACK. (a) exact fit of the 2-3 turn; (b) every sign choice; (c) random directions as a control")
def scan(u, ks, sd_list=None, su_list=None):
    out = []
    for sd in (sd_list or list(itertools.product((1, -1), repeat=2))):
        for k in ks:
            D = solve_ladder(md, [*sd, 1], (u, k))
            if D is None: continue
            for su in (su_list or list(Us.keys())):
                U = Us[su]
                if U is None: continue
                out.append((k, sd, su, turns(U, D)))
    return out
for name in ("face 12+", "body ++-"):
    u = np.array(dirs[name], float); u /= np.linalg.norm(u)
    rows = scan(u, np.linspace(0.02, 0.3, 281))
    # for each sign pair, find k where V_cb crosses the measured value
    for sd in itertools.product((1, -1), repeat=2):
        for su in Us:
            rr = [r for r in rows if r[1] == sd and r[2] == su]
            for a, b2 in zip(rr, rr[1:]):
                if (a[3][1] - Vm[1])*(b2[3][1] - Vm[1]) < 0:
                    w = (Vm[1] - a[3][1])/(b2[3][1] - a[3][1]); t = a[3] + w*(b2[3] - a[3]); k = a[0] + w*(b2[0] - a[0])
                    pr(f"  {name}, down signs {sd}, up signs {su}: 2-3 turn set to 0.04183 at strength {k:.4f} -> 1-3 turn {t[2]:.5f} ({t[2]/Vm[2]:.2f}x), 1-2 turn {t[0]:.4f} ({t[0]/Vm[0]:.2f}x)")
rng = np.random.default_rng(1); hits10 = hits5 = 0; Nr = 300; hits_cb_ub = 0
for i in range(Nr):
    u = rng.normal(size=3); u /= np.linalg.norm(u)
    best = None
    for k in np.concatenate([np.linspace(-0.4, -0.01, 25), np.linspace(0.01, 0.4, 25)]):
        for sd in itertools.product((1, -1), repeat=2):
            D = solve_ladder(md, [*sd, 1], (u, k))
            if D is None: continue
            for su, U in Us.items():
                if U is None: continue
                t = turns(U, D); wm = max(abs(np.log(t[1:]/Vm[1:])))
                if best is None or wm < best: best = wm
    if best is not None and best < math.log(1.10): hits_cb_ub += 1
pr(f"  control: {Nr} random fold directions, one strength each: {hits_cb_ub} ({100*hits_cb_ub/Nr:.0f}%) bring BOTH the 2-3 and 1-3 turns within 10%.")
pr("  If many random directions do it, the cube direction is not special and the fit says little.")
