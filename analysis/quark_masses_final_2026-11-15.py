"""Open item 2 (heavy quark masses and the quark tilt), final pass.
(a) the pair/single difference as the cause of the two bubble shapes; (b) ONE combined candidate rule (democratic base, real;
fold through the far-side D in Nelson-Barr form with the reverse block zero; rule 4's one i on the link's light-corner part);
(c) attempts to derive the shape numbers and the fold strength; (d) summary numbers.
FOOTING (fixed before any test): every mass at the top's energy (MS-bar, mu = m_t), MeV, as in earlier tests A rule that claims an
E6 / far-side origin is also re-tested at the far-side scale 2.1e8 GeV, one-loop running. 'As weighed' is a mixed footing
and is not used. Errors: PDG 2024 relative errors with 1% added for the running. Tilts: PDG 2024 global fit with its own errors."""
import math, itertools, numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
np.set_printoptions(precision=4, suppress=True)
J1 = np.ones((3, 3)); mu = np.array([1.242, 627.8, 162500.0]); md = np.array([2.678, 53.36, 2873.0])
rel_u = np.array([0.034, 0.010, 0.006]); rel_d = np.array([0.018, 0.014, 0.005])     # PDG 2024 errors (+1% running)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857]); Vs = np.array([0.00068, 0.00074, 0.000088, 0.00018])   # |Vus| |Vcb| |Vub| |Vtd|
F = np.array([[1, 1, 1], [-1, 1, 1], [0, -2, 1]], float)/np.array([math.sqrt(2), math.sqrt(6), math.sqrt(3)])  # hierarchical basis
SG = [(1, -1), (-1, 1)]; TOP = "top's energy"
base = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]]))
def fit_base(m, signs, start=None):
    tgt = np.sort(np.array(m)*np.array([*signs, 1])); f = lambda p: np.sort(np.linalg.eigvalsh(base(p))) - tgt
    starts = [start] if start is not None else [[m[2]/3, d0, e0] for d0 in (0.0001, 0.001, 0.01, 0.05) for e0 in (-0.5, -0.2, 0.01, 0.2, 0.5)]
    for s in starts:
        p, info, ok, _ = fsolve(f, s, full_output=True)
        if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]: return p
    return None
def Hd(p, kap, u):
    m = base(p); return m @ (np.eye(3) - kap*np.outer(u, u.conj())) @ m
def obs(V):
    A = np.abs(V); t = np.array([A[0,1], A[1,2], A[0,2], A[2,0]])
    Jc = (V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1])).imag
    al = np.degrees(np.angle(-V[2,0]*np.conj(V[2,2])/(V[0,0]*np.conj(V[0,2]))))
    be = np.degrees(np.angle(-V[1,0]*np.conj(V[1,2])/(V[2,0]*np.conj(V[2,2]))))
    ga = np.degrees(np.angle(-V[0,0]*np.conj(V[0,2])/(V[1,0]*np.conj(V[1,2]))))
    return t, Jc, al, be, ga
def ckm(Uu, H): return Uu.conj().T @ np.linalg.eigh(H)[1]
def up_U(signs, m=mu):
    p = fit_base(m, signs); w, O = np.linalg.eigh(base(p)); return O[:, np.argsort(np.abs(w))], p
worst = lambda t, V=Vm: math.exp(max(abs(np.log(t/V))))

pr("="*100); pr("0. FOOTING AND INPUTS"); pr("="*100)
pr(f"masses at the top's energy (MeV): up ladder {mu}, down ladder {md}")
pr(f"relative errors used: up {rel_u}, down {rel_d}; tilts {Vm} +- {Vs}")
# one-loop running from m_t to the far-side scale (Machacek & Vaughn 1984; Arason and others 1992)
def rge(t, y):
    g1, g2, g3, yt = y; k = 1/(16*math.pi**2)
    return [k*41/10*g1**3, -k*19/6*g2**3, -k*7*g3**3, k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - 0.85*g1**2), k*yt**2]
y0 = [math.sqrt(4*math.pi/58.9), math.sqrt(4*math.pi/29.6), math.sqrt(4*math.pi*0.108), 162.5/174.1]
tfar = math.log(2.1e8/162.5); sol = solve_ivp(lambda t, y: rge(t, y[:4]), [0, tfar], y0 + [0.0], rtol=1e-9, atol=1e-12)
I_far = sol.y[4, -1]; ef = math.exp(1.5*I_far)
pr(f"running integral I = int y_t^2/(16 pi^2) dln(mu) from m_t to 2.1e8 GeV = {I_far:.4f}; so m_c/m_t falls by {1/ef:.3f}x, m_s/m_b rises by {ef:.3f}x,")
pr(f"  and the small tilts |Vcb|, |Vub|, |Vtd| fall by {1/ef:.3f}x; ratios inside the light pairs (u/c, d/s) and |Vus| do not move.")
mu_far = mu*np.array([1, 1, ef]); md_far = md*np.array([1, 1, 1/ef]); Vm_far = Vm*np.array([1, 1/ef, 1/ef, 1/ef]); Vs_far = Vs*np.array([1, 1/ef, 1/ef, 1/ef])

pr(); pr("="*100); pr("1. (a) PAIR versus SINGLE as the cause of the two bubble shapes"); pr("="*100)
shape = {}
for name, uu, dd in (("top's energy", mu, md), ("far-side 2.1e8 GeV", mu_far, md_far)):
    for s in SG:
        pu, pd = fit_base(uu, s), fit_base(dd, s); shape[(name, s)] = (pu, pd)
        pr(f"{name:20s} signs {s}: up split {pu[1]:.4e} break {pu[2]:+.5f} | down split {pd[1]:.4e} break {pd[2]:+.5f} | split ratio {abs(pd[1]/pu[1]):.2f}, break ratio {abs(pd[2]/pu[2]):.3f}")
rng = np.random.default_rng(2026); mc_ratio = {}
for s in SG:
    pu0, pd0 = shape[("top's energy", s)]; sr, br, xs = [], [], []
    for i in range(400):
        uu = mu*(1 + rel_u*rng.normal(size=3)); dd = md*(1 + rel_d*rng.normal(size=3))
        pu, pd = fit_base(uu, s, pu0), fit_base(dd, s, pd0)
        if pu is None or pd is None: continue
        sr.append(abs(pd[1]/pu[1])); br.append(abs(pd[2]/pu[2])); xs.append((pu[1]/pu[2]**2, abs(pd[1])/pd[2]**2, pu[1], pu[2], pd[1], pd[2]))
    mc_ratio[s] = (np.mean(sr), np.std(sr), np.mean(br), np.std(br), np.array(xs))
    pr(f"with the mass errors (400 draws), signs {s}: split ratio {np.mean(sr):.2f} +- {np.std(sr):.2f}, break ratio {np.mean(br):.3f} +- {np.std(br):.3f}")
r23u, r23d, r12u, r12d = mu[1]/mu[2], md[1]/md[2], mu[0]/mu[1], md[0]/md[1]
pr(f"\nSame facts in plain mass ratios (top's energy): m_c/m_t {r23u:.5f}, m_s/m_b {r23d:.5f} (ratio {r23d/r23u:.2f}); m_u/m_c {r12u:.6f}, m_d/m_s {r12d:.5f} (ratio {r12d/r12u:.1f})")
pr("\nRule PS1 (Clebsch counting, one Higgs link; Gursey, Ramond & Sikivie 1976; Georgi & Jarlskog 1979):")
pr("  in E6 (27x27x27) and SO(10) (16x16x10) a single link gives up matrix = down matrix times one number. The 1-versus-3 Clebsch factor")
pr("  of Georgi & Jarlskog acts between quarks and leptons, never between up and down. Prediction: split ratio 1, break ratio 1.")
pr(f"  measured {mc_ratio[SG[0]][0]:.1f} and {mc_ratio[SG[0]][2]:.2f}: FAIL (by 24x and 4.5x).")
pr("Rule PS2 (the pair is counted twice: up mass = pair x pair, down mass = pair x single; Georgi & Glashow 1974; Froggatt & Nielsen 1979;")
pr("  Bando & Kugo 1999). If the family steps live in the pairs, each up step is the SQUARE of the down step.")
s12 = math.sqrt((2*math.hypot(rel_d[0], rel_d[1]))**2 + rel_u[0]**2 + rel_u[1]**2); s23 = math.sqrt((2*math.hypot(rel_d[1], rel_d[2]))**2 + rel_u[1]**2 + rel_u[2]**2)
pr(f"  light step: m_u/m_c = {r12u:.6f} against (m_d/m_s)^2 = {r12d**2:.6f}: {r12d**2/r12u:.2f}x, {math.log(r12d**2/r12u)/s12:.1f} sigma (same on every footing); exponent {math.log(r12u)/math.log(r12d):.2f} (rule: 2)")
pr(f"  heavy step: m_c/m_t = {r23u:.6f} against (m_s/m_b)^2 = {r23d**2:.6f}: {r23u/r23d**2:.1f}x off, {math.log(r23u/r23d**2)/s23:.0f} sigma; exponent {math.log(r23u)/math.log(r23d):.2f}")
pr(f"  heavy step at the far-side scale: {r23u/ef/(r23d*ef)**2:.1f}x off; exponent {math.log(r23u/ef)/math.log(r23d*ef):.2f}")
pr("  FAIL as a rule (4 sigma and 11x). Holds in kind only: the up ladder is steeper, with exponents 2.1 and 1.4 instead of 2 and 2.")
pr("Rule PS3/PS4 (counting directions: a pair is 2 colour directions, a single is 1; right-side hypercharge 2/3 against 1/3; electric pull 4/9 against 1/9):")
pr("  the only numbers this gives are 2 and 4. Rule PS5 (colour count, Georgi-Jarlskog's 3): 3 and 9.")
for s in SG:
    a, sa, b, sb, _ = mc_ratio[s]
    pr(f"  signs {s}: break ratio {b:.3f}+-{sb:.3f}: vs 2: {(b-2)/sb:+.0f} sigma, vs 3: {(b-3)/sb:+.0f}, vs 4: {(b-4)/sb:+.1f}, vs 9: {(b-9)/sb:+.0f} | split ratio {a:.2f}+-{sa:.2f}: vs 4: {(a-4)/sa:+.0f}, vs 9: {(a-9)/sa:+.0f}, vs 16: {(a-16)/sa:+.0f}, vs 27: {(a-27)/sa:+.1f}")
pr("  FAIL: no pre-stated number lands within 2 sigma on either sign choice.")
# luck check for after-the-fact matches among simple numbers 2^a 3^b
S = sorted({2.0**a*3.0**b for a in range(-2, 6) for b in range(0, 4) if 1 < 2.0**a*3.0**b < 100})
pr(f"\nLuck check for matches found AFTER looking (for example 9/2 for the break, 24 for the split). Simple numbers 2^a x 3^b between 1 and 100: {len(S)} of them.")
hits = []
for (name, s), (pu, pd) in shape.items():
    for lab, val, sg in (("split", abs(pd[1]/pu[1]), mc_ratio[s][1]/mc_ratio[s][0]), ("break", abs(pd[2]/pu[2]), mc_ratio[s][3]/mc_ratio[s][2])):
        near = min(S, key=lambda q: abs(math.log(q/val))); z = math.log(val/near)/sg; hits.append(abs(z) < 2)
        pr(f"  {name:20s} signs {s} {lab}: {val:.3f}, nearest simple number {near:g}, {z:+.1f} sigma")
lr = rng.uniform(0, math.log(100), size=(20000, 2)); sgs = (mc_ratio[SG[0]][1]/mc_ratio[SG[0]][0], mc_ratio[SG[0]][3]/mc_ratio[SG[0]][2])
lS = np.log(S); hit1 = np.array([[np.min(np.abs(lS - x)) < 2*sg for x, sg in zip(row, sgs)] for row in lr])
p_both = np.mean(hit1[:, 0] & hit1[:, 1]); p_one = np.mean(hit1, axis=0)
pr(f"  a random number matches some simple number within 2 sigma: split {100*p_one[0]:.0f}%, break {100*p_one[1]:.0f}%; both at once {100*p_both:.1f}%;")
pr(f"  with 4 chances (2 footings x 2 sign choices) at least one double match: {100*(1-(1-p_both)**4):.0f}%. Seen: {'a' if any(hits[i] and hits[i+1] for i in range(0, 8, 2)) else 'no'} double match.")

pr(); pr("="*100); pr("2. (a, c) WHAT E6 ITSELF SAYS ABOUT PAIRS AND SINGLES: the singles have a far-side twin, the pairs do not"); pr("="*100)
pr("In the 27 the right-side down partner (a single) has a twin of the same charges in the 10 (D's partner); the right-side up partner")
pr("(a pair) has none. With the reverse block empty (Nelson 1984; Barr 1984; Bento, Branco & Parada 1991) the heavy D mixes only the")
pr("right-side singles, and the light down matrix becomes   H_d = m (1 - kappa u u^dagger) m^T,   kappa = |F|^2/(|F|^2 + M_D^2),")
pr("with m the real base, F the far-side link (direction u) and M_D the D's own weight. This is a press from the RIGHT side.")
# 2.1 one link matrix for everything -> no tilt at all
lam = base([1.0, 0.03, 0.2]); lam3 = lam.copy()
Fm, Mm = 3.0*lam, 1.7*lam; K = Fm.conj().T @ np.linalg.inv(Fm @ Fm.conj().T + Mm @ Mm.conj().T) @ Fm
Hd1 = lam @ (np.eye(3) - K) @ lam.T; V1 = np.linalg.eigh(lam @ lam.T)[1].T @ np.linalg.eigh(Hd1)[1]
pr(f"\n2.1 One 27x27x27 link matrix for the base, the far-side link and D's weight (three D's, one per family): the press is {K[0,0]:.4f} x the unit matrix,")
pr(f"    largest off-diagonal tilt {np.max(np.abs(np.abs(V1) - np.eye(3))):.1e}. Exactly no tilt and identical ladder shapes: FAIL. A second, independent link matrix is needed,")
pr("    and its entries are new free numbers (this is Bento, Branco & Parada's minimal model). E6 does not supply them.")
# 2.2 same bubble (down base = c x up matrix) + one D pressing from the right: closed form
pr("\n2.2 Same bubble for both ladders (down base = c x up matrix), one D pressing the singles. In the up ladder's own basis")
pr("    H_d = c^2 D_u (1 - kappa w w^dagger) D_u = c^2 D_u^2 - v v^dagger: a rank-one change. Its eigenvectors are v_i/(lam_i - lam'_j)")
pr("    (Loewner 1934; Golub 1973), so every phase of w drops out: J = 0 EXACTLY, whatever the timing. The masses fix |w| for each c,")
pr("    leaving ONE free number (c) for the four tilts.")
def cauchy(c):
    lam_ = (c*mu)**2; lp = md**2
    v2 = np.array([-np.prod(lp - lam_[i])/np.prod([lam_[k] - lam_[i] for k in range(3) if k != i]) for i in range(3)])
    if np.any(v2 < 0): return None
    V = np.array([[math.sqrt(v2[i])/abs(lam_[i] - lp[j]) for j in range(3)] for i in range(3)]); V /= np.linalg.norm(V, axis=0)
    return V, 1 - np.prod(lp)/np.prod(lam_)
lo, hi = md[0]/mu[0], md[2]/mu[1]
pr(f"    interlacing (Cauchy) allows c only between m_d/m_u = {lo:.3f} and m_b/m_c = {hi:.3f}.")
bestc = None
for c in np.linspace(lo*1.0001, hi*0.9999, 4000):
    r = cauchy(c)
    if r is None: continue
    V, kap = r; t = np.array([V[0,1], V[1,2], V[0,2], V[2,0]]); w = worst(t)
    if bestc is None or w < bestc[0]: bestc = (w, c, t, kap)
pr(f"    best c = {bestc[1]:.3f}: tilts {np.round(bestc[2], 5)} (measured {Vm}), worst {bestc[0]:.1f}x, 1 - kappa = {1-bestc[3]:.1e}, J = 0.")
pr("    FAIL: the two bubbles are different shapes before the press (as 81.5 found for the other form of the press).")
pr("    With all three D's pressing, 1 - K is any 3x3 matrix (9 numbers for 7 measured): it fits everything and predicts nothing.")

pr(); pr("="*100); pr("3. (b) THE ONE COMBINED CANDIDATE RULE"); pr("="*100)
pr("Pieces: (i) democratic base p0 (J + diag(-d, d, e)), REAL, for both ladders (CP exact before the facing is chosen);")
pr("(ii) the fold through D in Nelson-Barr form, reverse block zero: H_d = m (1 - kappa u u^dagger) m, on the down ladder only;")
pr("(iii) rule 4's one i: the far-side link's LIGHT-CORNER part is a quarter turn out of step: u = (i s1, s2, s3) in the hierarchical basis.")
pr("Free numbers: up p0, d, e (3); down p0, d, e (3); press kappa (1); link direction (2 angles). Total 9 for 10 measured (6 masses,")
pr("3 tilts, 1 twist): ONE prediction. The six masses are inputs, so the test is: 3 numbers (kappa, 2 angles) against 4 tilt sizes.")
ups = {s: up_U(s) for s in SG}
def mk_u(basis, ph, a, b):
    s = np.array([math.cos(a)*math.cos(b)*np.exp(1j*ph[0]), math.cos(a)*math.sin(b)*np.exp(1j*ph[1]), math.sin(a)*np.exp(1j*ph[2])])
    return F @ s if basis == "hier" else s
def fit_nb(basis, ph, su, n, seed, kapfix=None, mdn=md, Uu=None, Vt=Vm, Vsg=Vs, allmin=False, dirmap=None):
    """3 free numbers (kappa, 2 angles) + the down base (3, pinned by the 3 down masses). Staged least squares from n random starts
    for each of the two down-base sign branches. Returns (chi2 on the four tilts, x, (kappa, u, H_d))."""
    r_ = np.random.default_rng(seed); Uu = ups[su][0] if Uu is None else Uu; best = None; mins = []
    def parts(x):
        kap = kapfix if kapfix is not None else 1/(1 + math.exp(-x[3])); a_, b_ = (x[4], x[5]) if dirmap is None else dirmap(x[4]); u = mk_u(basis, ph, a_, b_); return kap, u, Hd(x[:3], kap, u)
    def res(x, wt):
        kap, u, H = parts(x); w = np.sqrt(np.abs(np.linalg.eigvalsh(H))); o = obs(ckm(Uu, H))
        return np.concatenate([wt*(np.log(w) - np.log(mdn)), (o[0] - Vt)/Vsg])
    for sd in SG:
        p0 = fit_base(mdn, sd)
        for trial in range(n):
            x0 = [*p0*(1 + 0.3*r_.normal(size=3)), r_.uniform(2, 8), r_.uniform(0.9, 1.5)*r_.choice([-1, 1]), r_.uniform(0, 2*math.pi)]
            try: r = least_squares(res, x0, args=(30.0,), max_nfev=300); r = least_squares(res, r.x, args=(1e3,), max_nfev=300)
            except Exception: continue
            c = float(np.sum(r.fun**2)); mins.append(round(c, 1))
            if best is None or c < best[0]: best = (c, r.x, parts(r.x))
    return (best, sorted(set(mins))) if allmin else best
QT = (math.pi/2, 0.0, 0.0)
pr("\n3.1 The candidate, for each sign choice of the up base (the down base's numbers are free, so its sign branch is found by the fit):")
cand = {}
for su in SG:
    b, mins = fit_nb("hier", QT, su, 60, 11, allmin=True); cand[su] = b; chi, x, (kap, u, H) = b; t, Jc, al, be, ga = obs(ckm(ups[su][0], H))
    pr(f"  up {su}: best chi2 {chi:6.2f} (1 degree of freedom) | 1-kappa {1-kap:.2e} | tilts {np.round(t, 5)} | worst {worst(t):.3f}x | |J| {abs(Jc):.2e} | alpha {abs(al):.1f} beta {abs(be):.1f} gamma {abs(ga):.1f}")
    pr(f"      other local best points found (chi2): {mins[1:6]}")
su0 = min(cand, key=lambda q: cand[q][0]); chi0, x0, (kap0, u0, H0) = cand[su0]; sd0 = (1, -1) if x0[2] < 0 else (-1, 1)
t0, J0, al0, be0, ga0 = obs(ckm(ups[su0][0], H0)); s0 = F.T @ u0; s0 = s0*np.exp(-1j*np.angle(s0[2]))
pr(f"\n3.2 Best: up signs {su0}, down branch {sd0}. The ten numbers:")
wd = np.sqrt(np.abs(np.linalg.eigvalsh(H0)))
pr(f"  six masses: up {mu} (inputs, exact), down {np.round(wd, 3)} (inputs, exact)")
for nm, v, m_, s_ in zip(("|Vus|", "|Vcb|", "|Vub|", "|Vtd|"), t0, Vm, Vs): pr(f"  {nm}: {v:.5f}  measured {m_:.5f} +- {s_:.5f}  pull {(v - m_)/s_:+.2f} sigma")
from math import erfc
pr(f"  chi2 = {chi0:.2f} for 1 degree of freedom: {math.sqrt(chi0):.1f} sigma, chance {100*erfc(math.sqrt(chi0/2)):.0f}% -> the single prediction is NOT contradicted")
pr(f"  twist |J| = {abs(J0):.3e} (measured 3.12e-5 +- 0.13e-5: {(abs(J0) - 3.12e-5)/0.13e-5:+.1f} sigma). Its SIGN is not predicted (i or -i).")
pr(f"  triangle: alpha {abs(al0):.1f} deg (measured 91.6 +- 1.5 fit; 85.2 +- 4.6 direct), beta {abs(be0):.1f} (22.6 +- 0.4 from sin 2beta = 0.709 +- 0.011; model sin 2beta = {math.sin(math.radians(2*abs(be0))):.3f}, {(math.sin(math.radians(2*abs(be0))) - 0.709)/0.011:+.1f} sigma),")
pr(f"            gamma {abs(ga0):.1f} (measured 65.9 +- 3.0 direct: {(abs(ga0) - 65.9)/3.0:+.1f} sigma). |Vtd/Vts| = {t0[3]/np.abs(ckm(ups[su0][0], H0))[2,1]:.4f} (measured 0.205 +- 0.006)")
pr(f"  fitted numbers: down base split {x0[1]:+.5f}, break {x0[2]:+.5f} (the unpressed fit had {shape[(TOP, sd0)][1][1]:+.5f}, {shape[(TOP, sd0)][1][2]:+.5f});")
wb = np.sort(np.abs(np.linalg.eigvalsh(base(x0[:3]))))
pr(f"    the down base BEFORE the press would weigh {np.round(wb, 1)} MeV (span {wb[2]/wb[0]:.0f}); after the press {md} (span {md[2]/md[0]:.0f}); up span {mu[2]/mu[0]:.0f}.")
pr(f"    So the press lowers the three down masses by {wb[0]/md[0]:.1f}x, {wb[1]/md[1]:.1f}x, {wb[2]/md[2]:.1f}x and makes the ladder steeper, not flatter (answers 79.1's open test).")
pr(f"    press kappa = {kap0:.5f}, so D's own weight is {math.sqrt((1-kap0)/kap0):.3f} of its link (the link is {math.sqrt(kap0/(1-kap0)):.0f}x D's own weight);")
pr(f"    link direction in the hierarchical basis: |s| = {np.round(np.abs(s0), 4)}, phases {np.round(np.degrees(np.angle(s0)), 1)} deg;")
ud = u0*np.exp(-1j*np.angle(u0[2])); pr(f"    in the placement basis: {np.round(ud/abs(ud[0].real), 3)}  = (1, 1, {ud[2].real/ud[0].real:.3f}) + {abs(ud[0].imag/ud[0].real):.3f} i (1, -1, 0):")
pr("    the same three-part shape as a bubble (equal + third-placement break + split of the first two), with the split a quarter turn off.")
pr(f"    angle between the link and the old fold direction (-1,-1,2)/sqrt6: {math.degrees(math.acos(abs(s0[1]))):.0f} deg; and the all-equal direction: {math.degrees(math.acos(abs(s0[2]))):.0f} deg.")

pr("\n3.3 Where the one i sits: all six single placements (3 numbers against 4 tilts each, best sign choice)")
place = {}
for basis in ("hier", "dem"):
    for ip in range(3):
        ph = [0.0, 0.0, 0.0]; ph[ip] = math.pi/2
        b = min(((*fit_nb(basis, ph, su, 20, 21 + ip), su) for su in SG), key=lambda z: z[0]); t = obs(ckm(ups[b[3]][0], b[2][2]))[0]; place[(basis, ip)] = b[0]
        pr(f"  i on part {ip+1} of the {'hierarchical' if basis == 'hier' else 'placement'} basis: chi2 {b[0]:8.1f}, worst {worst(t):.3f}x, tilts {np.round(t, 5)}")
pr(f"  No phase at all (real link): J = 0, FAIL. Only the light-corner placement (hierarchical part 1) passes; the next best is {math.sqrt(sorted(place.values())[1]):.1f} sigma off.")

pr(); pr("="*100); pr("4. ATTACKS ON THE COMBINED CANDIDATE"); pr("="*100)
pr("4.1 Luck: put RANDOM timings on the link's first two parts (instead of a quarter turn and zero) and refit the same 3 numbers.")
r4 = np.random.default_rng(404); lucky = []
for i in range(80):
    ph = (*r4.uniform(0, 2*math.pi, 2), 0.0); lucky.append(min(fit_nb("hier", ph, su, 5, 1000 + i)[0] for su in SG))
lucky = np.array(lucky); share = np.mean(lucky <= chi0 + 0.05)
pr(f"  80 random timing pairs: {100*share:.0f}% do as well as the quarter turn (chi2 <= {chi0:.2f}); {100*np.mean(lucky <= 4):.0f}% reach chi2 <= 4; median chi2 {np.median(lucky):.0f}.")
pr(f"  (the random fits use fewer starts than the candidate's, so these shares are, if anything, too LOW.)")
pr(f"  The quarter turn on the light corner was named in an earlier test before this test, so the odds for it alone are about {100*share:.0f}%. Counting the 6 placements")
pr(f"  of 3.3 as 6 tries: {100*(1 - (1 - share)**6):.0f}%. A general timing pair (5 numbers for 4 tilts) always fits exactly, so this is the size of the 'good' band.")
pr("\n4.2 Data with the mass errors included (all six masses allowed to move inside their errors): 9 numbers, 10 measured.")
pu0 = ups[su0][1]
def res_all(x):
    wu = np.sort(np.abs(np.linalg.eigvalsh(base(x[:3])))); kap = 1/(1 + math.exp(-x[6])); u = mk_u("hier", QT, x[7], x[8]); H = Hd(x[3:6], kap, u)
    wdd = np.sqrt(np.abs(np.linalg.eigvalsh(H))); Ou = np.linalg.eigh(base(x[:3])); Uu = Ou[1][:, np.argsort(np.abs(Ou[0]))]
    return np.concatenate([(wu/mu - 1)/rel_u, (wdd/md - 1)/rel_d, (obs(ckm(Uu, H))[0] - Vm)/Vs])
ra = least_squares(res_all, [*pu0, *x0], x_scale=np.abs([*pu0, *x0]) + 1e-6, max_nfev=3000); chi_all = float(np.sum(ra.fun**2))
pr(f"  chi2 = {chi_all:.2f} for 1 degree of freedom ({math.sqrt(chi_all):.1f} sigma, chance {100*erfc(math.sqrt(chi_all/2)):.0f}%). Pulls: masses {np.round(ra.fun[:6], 2)}, tilts {np.round(ra.fun[6:], 2)}")
pr("\n4.3 Footing: the same candidate with masses and tilts carried to the far-side scale (2.1e8 GeV), where the D lives.")
far = {}
for su in SG:
    Uf = up_U(su, mu_far)[0]; far[su] = fit_nb("hier", QT, su, 40, 77, mdn=md_far, Uu=Uf, Vt=Vm_far, Vsg=Vs_far)
suf = min(far, key=lambda q: far[q][0]); cf, xf, (kf, uf, Hf) = far[suf]; tf = obs(ckm(up_U(suf, mu_far)[0], Hf))
pr(f"  best chi2 {cf:.2f} (up signs {suf}); 1-kappa {1-kf:.2e}; tilts {np.round(tf[0], 5)} against {np.round(Vm_far, 5)}; alpha {abs(tf[2]):.1f}, gamma {abs(tf[4]):.1f}.")
pr(f"  {'The candidate survives the change of footing.' if cf < 6 else 'The candidate does NOT survive the change of footing (chi2 above 6).'}")
pr("\n4.4 The theorem's conditions (Nelson 1984; Barr 1984; Bento, Branco & Parada 1991), checked on the full 4x4 down matrix [[m, h], [F, M]].")
m0 = base(x0[:3]); detu = np.linalg.det(base(pu0)); detm = np.linalg.det(m0)
for Mh in (2.4e8, 2.1e11):
    Fv = math.sqrt(kap0)*Mh*u0.conj(); Mo = math.sqrt(1 - kap0)*Mh
    full = np.zeros((4, 4), complex); full[:3, :3] = m0; full[3, :3] = Fv; full[3, 3] = Mo
    UL, S, _ = np.linalg.svd(full); order = np.argsort(S); S = S[order]; UL = UL[:, order]
    V4 = ups[su0][0].T @ UL[:3, :]; t4 = np.abs(np.array([V4[0,1], V4[1,2], V4[0,2], V4[2,0]]))
    leak = abs(V4[0, 3])**2
    pr(f"  heavy partner at {Mh/1e3:.1e} GeV: light masses {np.round(S[:3], 3)} MeV, tilts {np.round(t4, 5)} (same as the 3x3 formula: {np.allclose(t4, t0, rtol=1e-4)});")
    pr(f"    first-row leak into the heavy partner |V_uD|^2 = {leak:.1e} (the measured first-row shortfall is 1.6e-3: the leak is {1.6e-3/leak:.0e} times too small to be it)")
pr(f"  tree-level strong angle: the timing sits only in F, which does not enter det = det(m) x M. sign det(up) = {np.sign(detu):+.0f}, sign det(m) = {np.sign(detm):+.0f}:")
pr(f"    theta-bar = {'0' if detu*detm > 0 else 'pi unless the sign of the real number M (D own weight) is taken negative, a free discrete choice; then 0'}, with J = {abs(J0):.2e} not zero. HOLDS as algebra.")
hdir = np.ones(3)/math.sqrt(3); slope = abs((math.sqrt(kap0/(1 - kap0))*u0.conj() @ np.linalg.inv(m0) @ hdir).imag)
slope_w = np.linalg.norm((math.sqrt(kap0/(1 - kap0))*u0.conj() @ np.linalg.inv(m0)).imag)
pr(f"  reverse block switched on: theta-bar = {slope:.1e} x h/MeV for h along the all-equal direction, {slope_w:.1e} x h/MeV for the worst direction.")
pr(f"    For theta-bar < 1e-10 (Abel and others 2020) h must stay below {1e-10/slope:.0e} MeV ({1e-10/slope_w:.0e} MeV worst case), that is {1e-10/slope/174100:.0e} to {1e-10/slope_w/174100:.0e}")
pr(f"    of the Higgs value. So the block must be empty by an exact rule, not by smallness.")
pr("  CP exact in the base: the candidate uses real bases, and the only timing is the link's quarter turn. Consistent.")
pr("  NOT consistent with earlier tests: there the fold was written M_d - k p0 u u^dagger, which needs BOTH blocks (it is h M^-1 F). With the")
pr("    reverse block empty the press acts from the right side only. The old fold numbers (k = 0.11 to 0.28, 'link = 0.34 of the balance point',")
pr("    direction (-1,-1,2)) therefore do not carry over: see 3.2 for the new ones.")
pr("\n4.5 What the candidate says about the triangle (the tilt sizes are refitted with alpha pinned):")
def res_pin(x, a0):
    kap = 1/(1 + math.exp(-x[3])); u = mk_u("hier", QT, x[4], x[5]); H = Hd(x[:3], kap, u); o = obs(ckm(ups[su0][0], H))
    return np.concatenate([1e3*(np.log(np.sqrt(np.abs(np.linalg.eigvalsh(H)))) - np.log(md)), (o[0] - Vm)/Vs, [(abs(o[2]) - a0)/0.02]])
prof = {}
for direction in (-1, 1):
    x = x0.copy()
    for a0 in np.arange(round(abs(al0)), round(abs(al0)) + direction*8, direction):
        r = least_squares(res_pin, x, args=(a0,), max_nfev=400); x = r.x
        kap = 1/(1 + math.exp(-x[3])); o = obs(ckm(ups[su0][0], Hd(x[:3], kap, mk_u("hier", QT, x[4], x[5])))); prof[int(a0)] = (float(np.sum(r.fun[3:7]**2)), abs(o[3]), abs(o[4]))
for a0 in sorted(prof): pr(f"  alpha {a0}: chi2 {prof[a0][0]:.2f}, beta {prof[a0][1]:.2f}, gamma {prof[a0][2]:.2f}")
ok = [a for a in prof if prof[a][0] <= chi0 + 1]; pr(f"  within 1 unit of the best: alpha {min(ok)} to {max(ok)} deg, gamma {prof[max(ok)][2]:.1f} to {prof[min(ok)][2]:.1f} deg. The quarter turn in the link gives 'alpha near 90',")
pr(f"  not exactly 90. The sharp number gamma = 67.4 +- 0.5 belongs to the separate statement 'alpha is exactly 90' (78.3); at alpha = 90 the candidate has chi2 {prof[90][0]:.2f}, gamma {prof[90][2]:.1f}.")

pr(); pr("="*100); pr("5. (c) CAN STRUCTURE GIVE THE NUMBERS?"); pr("="*100)
pr("5.1 The press kappa from the balance rule and other fixed values (2 numbers left for 4 tilts):")
for nm, kf, n_ in (("1/2  (link = D's own weight: the balance point)", 0.5, 8), ("1/3  (one placement of three)", 1/3, 8), ("2/3", 2/3, 8), ("7/9  (1 - the 2/9 twist)", 7/9, 8), ("8/9", 8/9, 8), ("26/27", 26/27, 8), ("0.99", 0.99, 8), ("0.9985", 0.9985, 25), ("0.9990", 0.9990, 25), ("0.9994", 0.9994, 25), ("0.9997", 0.9997, 25)):
    b = min((fit_nb("hier", QT, su, n_, 55, kapfix=kf) for su in SG), key=lambda z: z[0])
    pr(f"  kappa = {nm:48s}: chi2 {b[0]:8.1f}")
pr(f"  FAIL for every structural value. The data want 1 - kappa = {1-kap0:.1e} (D's own weight about {100*math.sqrt((1-kap0)/kap0):.0f}% of its link). No rule gives this number.")
pr("5.1b Link shapes with one number fewer (2 numbers, kappa and one angle, for 4 tilts):")
for nm, dm in (("no part along the old fold direction (s2 = 0): all-equal + i x split", lambda t: (t, 0.0)),
               ("no all-equal part (s3 = 0)", lambda t: (0.0, t)),
               ("no link to the third placement (C = 0)", lambda t: (math.asin(math.cos(t)*math.sqrt(2/3)), math.atan2(math.cos(t)/math.sqrt(3), math.sin(t))))):
    b = min((fit_nb("hier", QT, su, 15, 66, dirmap=dm) for su in SG), key=lambda z: z[0])
    pr(f"  {nm:70s}: chi2 {b[0]:8.1f}")
pr("  FAIL: the link needs all three parts.")
pr("5.2 The old (two-block) fold with link = balance point, k = 1 (the earlier identity k = (g/sqrt(p0 M_D))^2):")
def solve_add(m, signs, k, u):
    tgt = np.sort(np.array(m)*np.array(signs)); Pu = np.outer(u, u.conj()); bld = lambda p: base(p) - k*p[0]*Pu
    f = lambda p: np.sort(np.linalg.eigvalsh(bld(p))) - tgt
    for d0 in (0.001, 0.01, 0.05, 0.2):
        for e0 in (-0.5, -0.2, 0.01, 0.2, 0.5, 1.0, 2.0):
            p, info, okk, _ = fsolve(f, [m[2]/3, d0, e0], full_output=True)
            if okk == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]:
                w, O = np.linalg.eigh(bld(p)); return O[:, np.argsort(np.abs(w))]
    return None
for k in (1.0, 0.5):
    bw = None
    for phi in np.radians(np.arange(0, 180, 15)):
        u = np.array([-1, -1, 2*np.exp(1j*phi)])/math.sqrt(6)
        for sd in SG:
            D = solve_add(md, [*sd, 1], k, u)
            if D is None: continue
            for su in SG:
                t = obs(ups[su][0].conj().T @ D)[0]; w = worst(t)
                if bw is None or w < bw[0]: bw = (w, math.degrees(phi), t)
    pr(f"  k = {k}: " + (f"best worst-tilt {bw[0]:.2f}x (timing {bw[1]:.0f} deg), tilts {np.round(bw[2], 5)}: FAIL" if bw else "no solution reproduces the down masses: FAIL"))
pr("5.3 The four shape numbers against a list of structure numbers fixed beforehand (E6 and TRPG counts):")
Lst = {"1/3": 1/3, "2/9": 2/9, "1/9": 1/9, "1/8": 1/8, "1/10": 1/10, "1/12": 1/12, "1/16": 1/16, "1/24": 1/24, "1/27": 1/27, "2/27": 2/27, "1/36": 1/36, "1/45": 1/45, "1/72": 1/72, "1/78": 1/78, "1/248": 1/248}
lL = np.log(list(Lst.values())); r5 = np.random.default_rng(5); nhit = 0; exp_hits = 0.0
for s in SG:
    xs = mc_ratio[s][4]
    for lab, col in (("up split", 2), ("up break", 3), ("down split", 4), ("down break", 5)):
        val = abs(np.mean(xs[:, col])); sg = np.std(xs[:, col])/val; nm = min(Lst, key=lambda q: abs(math.log(Lst[q]/val))); z = math.log(val/Lst[nm])/sg
        draws = math.log(val) + r5.uniform(-math.log(3), math.log(3), 4000); ph = np.mean([np.min(np.abs(lL - x)) < 2*sg for x in draws]); exp_hits += ph; nhit += abs(z) < 2
        pr(f"  signs {s} {lab:10s}: {val:.5f} +- {100*sg:.1f}%; nearest {nm:5s} ({z:+.1f} sigma); a random value within 3x of it matches some list number {100*ph:.0f}% of the time")
pk = 1 - sum(math.exp(-exp_hits)*exp_hits**k/math.factorial(k) for k in range(nhit))
pr(f"  matches within 2 sigma: {nhit} (1/45 on one sign choice, 1/12 on the other: they cannot both hold); expected by luck: {exp_hits:.1f}; chance of {nhit} or more: {100*pk:.0f}%.")
pr("  Not evidence. And in the combined candidate the down base's own numbers change (3.2), so these two matches do not carry over.")
pr("5.4 A relation between the split and the break: split = 3 x break^2 (found after looking; 'three placements')")
for s in SG:
    xs = mc_ratio[s][4]; a, b = xs[:, 0], xs[:, 1]
    pr(f"  top's energy, signs {s}: up split/break^2 = {np.mean(a):.3f} +- {np.std(a):.3f} ({(np.mean(a)-3)/np.std(a):+.1f} sigma from 3); down {np.mean(b):.3f} +- {np.std(b):.3f} ({(np.mean(b)-3)/np.std(b):+.1f} sigma from 3)")
    pu, pd = shape[("far-side 2.1e8 GeV", s)]; pr(f"  far-side scale, signs {s}: up {pu[1]/pu[2]**2:.3f}, down {abs(pd[1])/pd[2]**2:.3f}")
Xu, Xd = mu[0]*mu[2]**2/mu[1]**3, md[0]*md[2]**2/md[1]**3
sXu = math.sqrt(rel_u[0]**2 + 4*rel_u[2]**2 + 9*rel_u[1]**2); sXd = math.sqrt(rel_d[0]**2 + 4*rel_d[2]**2 + 9*rel_d[1]**2)
pr(f"  in plain masses: m1 m3^2/m2^3 = {Xu:.1f} +- {Xu*sXu:.1f} (up), {Xd:.1f} +- {Xd*sXd:.1f} (down): {abs(Xu-Xd)/math.hypot(Xu*sXu, Xd*sXd):.1f} sigma apart at the top's energy;")
pr(f"    at the far-side scale {Xu*ef**2:.1f} and {Xd/ef**2:.1f} ({abs(Xu*ef**2-Xd/ef**2)/math.hypot(Xu*ef**2*sXu, Xd/ef**2*sXd):.1f} sigma apart); as weighed (mixed footing) 31 and 101.")
r6 = np.random.default_rng(6); st = r6.uniform(math.log(5), math.log(1000), size=(200000, 4))   # log steps m2/m1, m3/m2 for two random ladders
lx = (-st[:, 0] + 2*st[:, 1]) - (-st[:, 2] + 2*st[:, 3]); pl = np.mean(np.abs(lx) <= abs(math.log(Xu/Xd)))
pr(f"  luck: two random ladders (each step between 5x and 1000x) agree this well in m1 m3^2/m2^3 {100*pl:.1f}% of the time; with the 6 simplest such")
pr(f"    combinations one could have tried, {100*(1-(1-pl)**6):.0f}%. As a rule 'split = 3 x break^2' it FAILS on the down ladder (2.5 or 7.8 sigma, by sign choice)")
pr("    and moves with the footing. For the up ladder alone it is a coincidence to note (0.1 and -1.0 sigma), not evidence.")
pr("5.5 The octonion 3x3 (exceptional Jordan algebra; Jordan, von Neumann & Wigner 1934; Freudenthal 1954): its three invariants")
pr("  (trace, pair-sum, determinant) are exactly the three masses' three symmetric sums, so the algebra by itself leaves all three free.")
for nm, m_, s in (("up", mu, SG[0]), ("down", md, SG[0])):
    w = np.sort(m_*np.array([*s, 1])); tr, s2, dt = np.sum(w), w[0]*w[1] + w[1]*w[2] + w[0]*w[2], np.prod(w)
    Q = np.sum(m_)/np.sum(np.sqrt(m_))**2
    pr(f"  {nm:4s}: pair-sum/trace^2 = {s2/tr**2:+.5f}, determinant/trace^3 = {dt/tr**3:+.2e}; the all-equal shape (a primitive idempotent) has 0 and 0. Cone number (Koide 1981) {Q:.4f} (leptons: 2/3 = 0.6667)")
pr(f"  The natural added conditions fail: determinant = 0 needs a massless up quark (measured {mu[0]} +- {mu[0]*rel_u[0]:.3f} MeV, {1/rel_u[0]:.0f} sigma from 0);")
pr("  the 45-degree cone holds for leptons only. FAIL: the octonion 3x3 does not fix the split or the break.")

pr(); pr("="*100); pr("6. SUMMARY NUMBERS"); pr("="*100)
pr(f"(a) pair/single: no counting rule gives the split ratio {mc_ratio[SG[0]][0]:.1f} / {mc_ratio[SG[1]][0]:.1f} or the break ratio {mc_ratio[SG[0]][2]:.2f} / {mc_ratio[SG[1]][2]:.2f}; E6's own pair/single")
pr("    statement (singles have a far-side twin) with one bubble shape gives J = 0 exactly and tilts off by more than 100x. FAIL.")
pr(f"(b) combined candidate: 9 free numbers for 10 measured; its one test gives chi2 {chi0:.2f} ({math.sqrt(chi0):.1f} sigma), {chi_all:.2f} with mass errors, {cf:.2f} at the far-side footing;")
pr(f"    luck {100*share:.0f}% (alone) to {100*(1 - (1 - share)**6):.0f}% (six placements). gamma {prof[max(ok)][2]:.1f} to {prof[min(ok)][2]:.1f} deg (best {abs(ga0):.1f}); exactly right-angled: 67.4. LEAD, not a derivation.")
pr("(c) no structural value for the press, the link shape, or the four shape numbers survives. The quark numbers stay measured inputs.")
