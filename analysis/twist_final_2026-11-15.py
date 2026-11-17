#!/usr/bin/env python3
"""
twist_final_2026-11-15.py   OPEN ITEM 3: the twist (the same number as the electron mass).

Every known or plausible route to fixing the lepton twist delta in
    sqrt(m_i) = A (1 + sqrt2 cos(delta + 2 pi i/3))
is tested with high-precision arithmetic and a luck (look-elsewhere) check.

Footing: POLE (directly measured) masses unless a section says otherwise.
Inputs: PDG 2024 / CODATA 2022.
Credits: Koide 1981-83, 1990; Foot 1994; Brannen 2006; Rivero & Gsponer 2005; Sumino 2009; Xing & Zhang 2006;
Xing, Zhang & Zhou 2008; Huang & Zhou 2021; Gray, Broadhurst, Grafe & Schilcher 1990 (two-loop pole mass);
Tarrach 1981; Lindemann 1882 & Weierstrass 1885; Niven 1956; Landau, Abrikosov & Khalatnikov 1954; Shulga 2026;
Machacek & Vaughn 1983-85 (through the earlier numbers); Particle Data Group 2024; CODATA 2022.
Runtime: about a minute. Deterministic (fixed seeds).
"""
import math, sys, itertools
import numpy as np
import mpmath as mp
mp.mp.dps = 40
OUT = open(__file__.replace(".py", ".txt"), "w")
def pr(s=""):
    print(s); OUT.write(s + "\n")
def head(s):
    pr("\n" + "=" * 110); pr(s); pr("=" * 110)

# ---------------------------------------------------------------- inputs
me, sme = mp.mpf("0.51099895069"), mp.mpf("1.6e-10")        # MeV, CODATA 2022
mmu, smmu = mp.mpf("105.6583755"), mp.mpf("2.3e-6")         # MeV, PDG 2024
mtau, smtau = mp.mpf("1776.93"), mp.mpf("0.09")             # MeV, PDG 2024
mu_over_e, s_mu_over_e = mp.mpf("206.7682827"), mp.mpf("4.6e-6")   # CODATA 2022 (22 ppb): the sharpest lepton number
alpha = 1/mp.mpf("137.035999177")                           # CODATA 2022
MPl, sMPl_rel = mp.mpf("1.220890e22"), mp.mpf("1.1e-5")     # MeV, CODATA 2022 (G known to 2.2e-5)
MPl_red = MPl/mp.sqrt(8*mp.pi)
r2 = mp.sqrt(2); pi = mp.pi
ppm = mp.mpf("1e-6")

def fam(delta, r=r2):
    """the three x_k = 1 + r cos(delta + 2 pi k/3), sorted light to heavy; masses are x^2 (times a scale)"""
    return sorted([(1 + r*mp.cos(delta + 2*pi*k/3)) for k in range(3)])
def ratios(delta, r=r2):
    x = fam(delta, r); return (x[0]/x[1])**2, (x[2]/x[1])**2          # e/mu, tau/mu
def params(m):
    """scale M0, lean r, twist delta from three masses"""
    s = [mp.sqrt(v) for v in m]; A = sum(s)/3; x = [v/A - 1 for v in s]
    r = mp.sqrt(2*sum(v*v for v in x)/3); c3 = 4*x[0]*x[1]*x[2]/r**3
    return A*A, r, mp.acos(c3)/3
def Q(m): return sum(m)/sum(mp.sqrt(v) for v in m)**2
f = lambda v, n=10: mp.nstr(v, n)

# ======================================================================================================
head("SECTION 1. THE NUMBER ITSELF (pole masses; tau error propagated)")
M0, r_meas, d_meas = params([me, mmu, mtau])
dQ = (Q([me, mmu, mtau + smtau]) - Q([me, mmu, mtau - smtau]))/2
dd = (params([me, mmu, mtau + smtau])[2] - params([me, mmu, mtau - smtau])[2])/2
dr = (params([me, mmu, mtau + smtau])[1] - params([me, mmu, mtau - smtau])[1])/2
pr(f"1a. All three masses measured, nothing assumed:")
pr(f"    Q = {f(Q([me,mmu,mtau]),9)} +- {f(abs(dQ),2)}   (2/3: {f((Q([me,mmu,mtau])-mp.mpf(2)/3)/abs(dQ),2)} sigma)")
pr(f"    lean r = {f(r_meas,9)} +- {f(abs(dr),2)}  (sqrt2 = {f(r2,9)}: {f((r_meas-r2)/abs(dr),2)} sigma);  self/relationship weight = {f(2/r_meas**2,8)}")
pr(f"    twist delta = {f(d_meas,9)} +- {f(abs(dd),2)} rad  (2/9 = {f(mp.mpf(2)/9,9)}: {f((d_meas-mp.mpf(2)/9)/abs(dd),2)} sigma);  scale M0 = {f(M0,8)} MeV")
SIG2 = abs(dd)                                    # tier-2 error on delta: nothing assumed (tau error)

emu = 1/mu_over_e; s_emu = emu*s_mu_over_e/mu_over_e
d45 = mp.findroot(lambda d: ratios(d)[0] - emu, mp.mpf("0.2222"))
d45_hi = mp.findroot(lambda d: ratios(d)[0] - (emu + s_emu), mp.mpf("0.2222"))
SIG1 = abs(d45_hi - d45)                          # tier-1 error on delta: 45 degrees exact (e/mu error only)
tau45 = mmu*ratios(d45)[1]; tau45_err = mp.sqrt((smmu/mmu)**2 + (abs(mmu*ratios(d45_hi)[1]-tau45)/tau45)**2)*tau45
eps = 1 - d45/(mp.mpf(2)/9)
pr(f"\n1b. If the 45-degree rule (rule 7) is exact, the electron/muon ratio (known to 22 parts per billion) fixes the twist:")
pr(f"    delta = {f(d45,13)} +- {f(SIG1,2)} rad  = {f(d45*180/pi,9)} degrees")
pr(f"    2/9 - delta = {f(mp.mpf(2)/9-d45,4)} = {f(eps/ppm,5)} ppm of 2/9  ({f((mp.mpf(2)/9-d45)/SIG1,4)} sigma).   Call this shortfall eps = {f(eps,5)} +- {f(SIG1/(mp.mpf(2)/9),2)}")
pr(f"    tau mass then = {f(tau45,9)} +- {f(tau45_err,2)} MeV (Koide 1981);  measured {f(mtau,6)} +- {f(smtau,2)}: {f((tau45-mtau)/smtau,3)} sigma")
gap = pi/12 - d45
pr(f"    massless electron at pi/12 = {f(pi/12,10)}; gap = {f(gap,8)} rad = {f(gap*180/pi,6)} degrees; delta/(pi/12) = {f(d45/(pi/12),9)}")
x = fam(d45); pr(f"    electron's root-mass share x_e = sqrt(m_e/M0) = {f(x[0],9)};  m_e/M0 = {f(x[0]**2,9)};  m_e/m_tau = {f((x[0]/x[2])**2,9)}")

e29, t29 = ratios(mp.mpf(2)/9)
pr(f"\n1c. Both exact (45 degrees AND 2/9):  e/mu predicted {f(e29,12)} vs measured {f(emu,12)} +- {f(s_emu,2)}:")
pr(f"    off {f((e29/emu-1)/ppm,4)} ppm = {f((e29-emu)/s_emu,4)} sigma.  FAIL (as in an earlier test, which quoted 447 sigma with a slightly different error on the ratio).   tau would be {f(mmu*t29,8)} MeV.")
r29 = mp.findroot(lambda r: ratios(mp.mpf(2)/9, r)[0] - emu, r2)
NS = f(abs((e29-emu)/s_emu),3)
pr(f"\n1d. If instead 2/9 is exact and the lean gives way: r = {f(r29,12)} (sqrt2 {f((r29/r2-1)/ppm,3)} ppm), Q = {f((1+r29**2/2)/3,12)}, self/relationship = {f(2/r29**2,10)}")
tau29 = mmu*ratios(mp.mpf(2)/9, r29)[1]
pr(f"    tau mass then = {f(tau29,9)} MeV ({f((tau29-mtau)/smtau,3)} sigma).")
# self-consistent variant 3 delta = Q
rq = mp.findroot(lambda r: ratios((1+r*r/2)/9, r)[0] - emu, r2); dq = (1+rq*rq/2)/9
tauq = mmu*ratios(dq, rq)[1]
pr(f"1e. Variant 'three times the twist equals the Koide value' (3 delta = Q, one free number): r = {f(rq,12)}, delta = {f(dq,12)}, tau = {f(tauq,9)} MeV.")
pr(f"\n1f. What the tau mass can and cannot decide:")
pr(f"    45 exact:        {f(tau45,9)} MeV")
pr(f"    2/9 exact:       {f(tau29,9)} MeV   (differs by {f((tau29-tau45)*1000,3)} keV)")
pr(f"    3 delta = Q:     {f(tauq,9)} MeV   (differs by {f((tauq-tau45)*1000,3)} keV)")
pr(f"    measured error today: 90 keV. So the tau mass tests the 45-degree rule, and CANNOT tell which of these is right until it is known to about 1 keV.")
sens_me = (mp.findroot(lambda d: ratios(d)[0] - emu*(1+mp.mpf("1e-6")), d45) - d45)/mp.mpf("1e-6")
pr(f"    sensitivity: d(delta)/d ln(m_e) = {f(sens_me,6)} rad  (1 ppm in the electron mass moves the twist by {f(abs(sens_me)*1e-6/d45/ppm,4)} ppm)")

# ======================================================================================================
head("SECTION 2. ROUTE 1: 2/9 PLUS A RADIATIVE (QED) CORRECTION")
a_pi = alpha/pi
z2, z3 = mp.zeta(2), mp.zeta(3)
pr("2a. Footing A: running masses at ONE common energy mu. One-loop QED (Tarrach 1981): m(mu) = M [1 - (alpha/pi)(1 + (3/2) ln(mu/M))].")
for name, mu in [("m_tau", mtau), ("M_Z", mp.mpf(91187.6)), ("1e16 GeV", mp.mpf("1e19"))]:
    run = [M*(1 - a_pi*(1 + mp.mpf(3)/2*mp.log(mu/M))) for M in (me, mmu, mtau)]
    M0r, rr, dr_ = params(run)
    pr(f"    mu = {name:9s}: Q = {f(Q(run),8)} (2/3 {f((Q(run)*3/2-1)*100,3)}%), lean off sqrt2 by {f((rr/r2-1)*100,3)}%, twist = {f(dr_,8)} (2/9 {f((dr_*9/2-1)*100,3)}%)")
hz = [mp.mpf("0.48307"), mp.mpf("101.766"), mp.mpf("1728.56")]
M0r, rr, dr_ = params(hz)
pr(f"    full Standard Model running to M_Z (Huang & Zhou 2021): Q = {f(Q(hz),8)}, twist = {f(dr_,8)} (2/9 {f((dr_*9/2-1)*100,3)}%)")
shift1 = mp.mpf(3)/2*a_pi*mp.log(mmu/me)
pr(f"    The mass-dependent part is mu-independent at this order: e/mu shifts by (3 alpha/2 pi) ln(m_mu/m_e) = {f(shift1*100,4)}% = {f(shift1/ppm,6)} ppm.")
pr(f"    Needed to turn exact (45, 2/9) into the data: +{f((emu/e29-1)/ppm,4)} ppm.  The one-loop shift is {f(shift1/(emu/e29-1),5)} times too big.  FAIL on any common-energy footing.")
pr("    (This is Sumino's 2009 point: the rules hold for the directly weighed masses, where QED says they should not, unless something cancels the log.)")

pr("\n2b. Footing B: each mass at its OWN scale, m_i(m_i), to two loops (Gray, Broadhurst, Grafe & Schilcher 1990, turned into QED: C_F = 1, C_A = 0, T = 1).")
d1 = mp.mpf(7)/128 - mp.mpf(15)/8*z2 - mp.mpf(3)/4*z3 + 3*z2*mp.log(2)
d3 = -mp.mpf(71)/96 - z2/2
d4 = -mp.mpf(143)/96 + z2
Dl = lambda r: pi**2/8*r - mp.mpf("0.597")*r**2 + mp.mpf("0.230")*r**3
pr(f"    two-loop numbers: d1 = {f(d1,6)}, per massless lighter lepton d3 = {f(d3,6)}, own loop d4 = {f(d4,6)}")
leps = [("e", me), ("mu", mmu), ("tau", mtau)]
R = {}
for i, (n, M) in enumerate(leps):
    lighter = leps[:i]
    a_i = alpha/(1 - (2*alpha/(3*pi))*sum(mp.log(M/ml) for _, ml in lighter))/pi     # MS-bar alpha at mu = M, leptons up to this one
    K = d1 + d4 + sum(d3 + 2*Dl(ml/M) for _, ml in lighter)
    R[n] = 1/(1 + a_i + a_i**2*K)
    pr(f"    {n:3s}: alpha(M)/pi = {f(a_i,9)} (1/alpha = {f(1/(a_i*pi),7)}), K = {f(K,5)}, m(M)/M = {f(R[n],12)}")
sh_emu = (R["mu"]/R["e"] - 1); sh_tmu = (R["mu"]/R["tau"] - 1)
need = emu/e29 - 1
pr(f"    If exact (45, 2/9) held for m_i(m_i), the pole ratio e/mu would be the family ratio times R_mu/R_e = 1 {f(sh_emu/ppm,4)} ppm.")
pr(f"    Needed: +{f(need/ppm,4)} ppm.  Got: {f(sh_emu/ppm,4)} ppm.  Same size (alpha^2 x log), OPPOSITE sign. Miss = {f((sh_emu-need)/ppm,4)} ppm = {f((sh_emu-need)*emu/s_emu,4)} sigma. FAIL.")
piece1 = -(alpha/(1 - (2*alpha/(3*pi))*mp.log(mmu/me)) - alpha)/pi; piece2 = sh_emu - piece1
pr(f"    The two pieces: alpha grows between the two scales: {f(piece1/ppm,4)} ppm (certain); the electron loop inside the muon: +{f(piece2/ppm,3)} ppm (two-loop number).")
pr(f"    Robustness: even if the two-loop number were dropped or doubled the shift stays negative ({f(piece1/ppm,4)} to {f((piece1+2*piece2)/ppm,3)} ppm). The sign cannot come out right.")
own = [me*R["e"], mmu*R["mu"], mtau*R["tau"]]
M0o, ro, do = params(own)
pr(f"    (own-scale masses: twist = {f(do,10)} with tau measured; and with 45 exact the e/mu ratio needs twist = {f(mp.findroot(lambda d: ratios(d)[0]-emu*R['e']/R['mu'], d45),10)}: further from 2/9, not closer.)")
pr("    Not included: hadron loops in alpha at the tau (they touch only the tau, which is not the sharp test), three-loop terms (below 1 ppm).")

pr("\n2c. Reverse question: what correction WOULD do it?  Suppose exact (45, 2/9) holds for m~_i = M_i (1 + c ln M_i).")
c_need = need/mp.log(mmu/me)
pr(f"    needed c = {f(c_need,5)} per unit of log  = {f(c_need/a_pi**2,5)} x (alpha/pi)^2  = {f(c_need/alpha**2,5)} x alpha^2  = {f(c_need/(1.5*a_pi),4)} x the one-loop QED log.")
tau_c = mmu*t29*(1 + c_need*mp.log(mtau/mmu))**-1
pr(f"    tau mass in that case: {f(tau_c,9)} MeV (differs from the 45-exact value by {f((tau_c-tau45)*1000,3)} keV).")
pr("    The size is that of a two-loop QED log, so a correction of this kind is not absurd. But c is one free number fitted to one datum:")
pr("    it predicts nothing until some rule gives 0.34 (alpha/pi)^2, and it also needs something to cancel the 1900-times-larger one-loop log.")
pr("    Sumino's 2009 cancellation uses a family gauge force with strength alpha/4. TRPG has no family force: the three generations are placements,")
pr("    the carriers are exactly 8|3|1, and rule 6 is 'nothing extra'. So that mechanism is not available on the main line.")

# ======================================================================================================
head("SECTION 3. ROUTES 2 AND 3: OCTONION / G2 / F4 / E6 PHASES, AND THREE-PLACEMENT (CIRCULANT) PHASES")
pr("3a. A theorem first (Lindemann 1882; Weierstrass 1885). If an angle in radians is a non-zero ALGEBRAIC number (any fraction is one),")
pr("    its cosine, sine and e^(i angle) are TRANSCENDENTAL. And any structure built from finite algebra (octonion structure constants,")
pr("    group characters, Clebsch-Gordan numbers, a circulant with algebraic entries) has algebraic entries and algebraic phases e^(i delta).")
pr("    So: NO finite algebraic structure can have a phase of exactly 2/9 radian. Group theory hands out rational fractions of a TURN, not of a radian.")
pr("    A rational radian can only come from a length-over-length or an integral (an area, a Berry-type phase, a loop correction).")
pr("    Check by computer: is cos(2/3) (the quantity the masses use, cos 3 delta) a root of a small polynomial?")
for target, nm in [(mp.cos(mp.mpf(2)/3), "cos(2/3)"), (mp.cos(2*pi/9), "cos(2 pi/9) [control: algebraic]")]:
    pol = mp.findpoly(target, 6, maxcoeff=2000, tol=mp.mpf("1e-30"))
    pr(f"      {nm:34s}: degree <= 6, coefficients <= 2000 -> {pol}")

pr("\n3b. Does the measured twist sit at a rational fraction of a turn (what group theory could give)?")
best = []
for q in range(1, 361):
    p = int(mp.nint(d45/(2*pi)*q))
    if p > 0 and math.gcd(p, q) == 1: best.append((abs(p*2*pi/q - d45), p, q))
best.sort()
for dist, p, q in best[:4]:
    pr(f"      {p}/{q} of a turn = {f(2*pi*p/q,10)} rad: off by {f(dist/d45*100,3)}%  = {f(dist/SIG2,4)} sigma (nothing assumed) / {f(dist/SIG1,3)} sigma (45 exact)")
nq = sum(1 for q in range(1, 361) for p in range(1, q) if math.gcd(p, q) == 1 and 2*math.pi*p/q < math.pi/12)
pr(f"    {nq} fractions of a turn with denominator <= 360 lie below pi/12; their mean spacing is {math.pi/12/nq:.2e} rad, the tau-error window (+-3 sigma) is {float(6*SIG2):.1e} rad.")
pr(f"    Expected hits by luck: {nq*float(6*SIG2)/(math.pi/12):.2f}. Found within 3 sigma: {sum(1 for d_,p,q in best if d_ < 3*SIG2)} (11/311, a denominator no structure offers; it fails by 32,000 sigma if 45 is exact).")
pr("    Smallest denominators: 1/28 of a turn is %.2f%% off, 1/29 is %.2f%% off. FAIL." % (100*(2*math.pi/28/float(d45)-1), 100*(2*math.pi/29/float(d45)-1)))

pr("\n3c. Group-theory numbers of the 3x3 octonion matrix that EQUAL 2/9 (as pure numbers):")
S = {"dim G2 (14)": 14, "dim F4 (52)": 52, "dim E6 (78)": 78, "dim SO(10) (45)": 45, "27 of E6": 27, "26 of F4": 26, "16": 16, "10": 10,
     "octonion (8)": 8, "imaginary units (7)": 7, "colours/placements (3)": 3, "h(G2)=6": 6, "h(F4)=h(E6)=12": 12, "dual h(G2)=4": 4,
     "dual h(F4)=9": 9, "two (2)": 2, "one (1)": 1, "rank E6 (6)": 6, "24 off-diagonal": 24, "dim SU(3) (8)": 8, "dim SO(8) (28)": 28}
vals = sorted(set(S.values()))
hits = [(a, b) for a in vals for b in vals if a < b and abs(a/b - 2/9) < 1e-12]
pr(f"      ratios a/b of these {len(vals)} integers equal to 2/9: {hits}   e.g. 2/(dual Coxeter number of F4) and 6/27.")
rat = sorted(set(a/b for a in vals for b in vals if a < b and a/b < math.pi/12))
allr = sorted(set(a/b for a in vals for b in vals if a < b))
pr(f"      but the same integers make {len(allr)} different ratios below 1, {len(rat)} of them below pi/12: {', '.join(f'{v:.4f}' for v in rat)}")
pr("      So some structural ratio lands on any simple fraction. And by 3a such a ratio is a number, not an angle: it needs the extra rule")
pr("      'this number is an angle in radians', which is the same missing unit that sank the ring reading (34.5). FAIL as a derivation.")

pr("\n3d. Is ANY simple algebraic value hit by the invariants the three placements actually use? (twist with 45 exact, +-3 sigma of each tier)")
def target_fns():
    return {"cos(3 delta)": lambda d: mp.cos(3*d), "sin(3 delta)": lambda d: mp.sin(3*d), "tan(3 delta)": lambda d: mp.tan(3*d),
            "cos(delta)": lambda d: mp.cos(d), "sin(delta)": lambda d: mp.sin(d), "tan(delta)": lambda d: mp.tan(d),
            "det of root-matrix e3 = cos(3d)/sqrt2 - 1/2": lambda d: mp.cos(3*d)/r2 - mp.mpf(1)/2,
            "x_e = sqrt(m_e/M0)": lambda d: fam(d)[0], "m_e/M0": lambda d: fam(d)[0]**2,
            "m_e/m_tau": lambda d: (fam(d)[0]/fam(d)[2])**2, "m_e/m_mu": lambda d: ratios(d)[0], "m_mu/m_tau": lambda d: 1/ratios(d)[1],
            "delta/(pi/12)": lambda d: d/(pi/12), "delta (radians)": lambda d: d, "delta/pi": lambda d: d/pi}
famA = set()
for p in range(1, 61):
    for q in range(1, 61):
        famA.add(p/q); famA.add(math.sqrt(p/q))
for c in (2, 3, 5, 6, 7):
    for a in range(-9, 10):
        for b in range(-9, 10):
            if b == 0: continue
            for d_ in range(1, 13):
                v = (a + b*math.sqrt(c))/d_
                if v > 0: famA.add(v)
famA = np.array(sorted(famA)); lf = np.log(famA)
pr(f"      family: p/q and sqrt(p/q) with p, q <= 60, and (a + b sqrt c)/d with |a|,|b| <= 9, c in 2,3,5,6,7, d <= 12: {len(famA)} distinct values")
rng = np.random.default_rng(20261115)
def count(v, tol):                       # members within relative tol of v
    lo, hi = np.searchsorted(lf, [math.log(v) - tol, math.log(v) + tol]); return hi - lo, famA[lo:hi]
tot1 = tot2 = 0; exp1 = exp2 = 0.0
pr(f"      {'invariant':46s} {'value':>14s} {'hits 45-exact':>14s} {'luck':>8s} {'hits tau-error':>15s} {'luck':>8s}  nearest simple value (off)")
for nm, fn in target_fns().items():
    v = fn(d45); s1 = abs(fn(d45 + SIG1) - v)/abs(v); s2 = abs(fn(d45 + SIG2) - v)/abs(v)
    v = float(v); t1, t2 = 3*float(s1), 3*float(s2)
    h1, _ = count(v, t1); h2, m2 = count(v, t2)
    dec = v*np.exp(rng.uniform(-0.3, 0.3, 4000))
    l1 = np.mean([count(u, t1)[0] for u in dec]); l2 = np.mean([count(u, t2)[0] for u in dec])
    near = famA[np.argmin(abs(lf - math.log(v)))]
    simple = min((abs(p/q - v)/v, f"{p}/{q}") for p in range(1, 13) for q in range(1, 13))
    tot1 += h1; tot2 += h2; exp1 += l1; exp2 += l2
    pr(f"      {nm:46s} {v:14.10f} {h1:14d} {l1:8.3f} {h2:15d} {l2:8.2f}  {simple[1]} ({simple[0]*1e6:.3g} ppm)")
pr(f"      totals: 45-exact tier {tot1} hits against {exp1:.2f} expected by luck; tau-error tier {tot2} hits against {exp2:.1f} expected by luck.")
pr("      Reading: the only small-number hit anywhere is delta = 2/9 itself (off 0.79 ppm). No algebraic invariant of the three placements is simple.")

pr("\n3e. Is the nearness to 2/9 itself luck?")
fr = sorted(set(p/q for q in range(1, 13) for p in range(1, q) if p/q < math.pi/12))
for label, dev in [("45 exact: |delta - 2/9| = %.2e" % float(mp.mpf(2)/9-d45), float(mp.mpf(2)/9 - d45)),
                   ("nothing assumed: tau error %.1e" % float(SIG2), float(max(abs(d_meas - mp.mpf(2)/9), SIG2)))]:
    p_luck = len(fr)*2*dev/(math.pi/12)
    pr(f"      {label}: {len(fr)} fractions with denominator <= 12 lie below pi/12; chance a random twist is this close to one: {p_luck:.1e}")
    pr(f"         allowing 10 ways of writing the angle (delta, 3 delta, per pi, per turn, tan, ...): {10*p_luck:.1e}")
pr(f"      So the nearness to 2/9 is NOT luck (about 1 in 10,000 even after generous trials). And yet exact 2/9 with exact 45 fails by {NS} sigma.")
pr("      Both statements are true together: the twist is 2/9 less a small, real, unexplained shortfall of 0.79 ppm.")

# ======================================================================================================
head("SECTION 4. '2/9 TIMES (1 - small correction)': IS ANY SIMPLE CORRECTION PICKED OUT?")
e0 = float(eps); se = float(SIG1/(mp.mpf(2)/9))
pr(f"    shortfall eps = {e0:.5e} +- {se:.1e} (45 exact).  Family: eps = k alpha^n pi^m, k = p/q with p, q <= 12, n = 1..4, m = -4..4.")
al = float(alpha); cands = {}
ks = sorted(set(p/q for p in range(1, 13) for q in range(1, 13)))
for n in range(1, 5):
    for m in range(-4, 5):
        for p in range(1, 13):
            for q in range(1, 13):
                if math.gcd(p, q) != 1: continue
                cands[p/q*al**n*math.pi**m] = f"({p}/{q}) alpha^{n} pi^{m}"
cv = np.array(sorted(cands)); lcv = np.log(cv)
def cnt(v, nsig):
    lo, hi = np.searchsorted(lcv, [math.log(v) - nsig*se/e0, math.log(v) + nsig*se/e0]); return hi - lo, cv[lo:hi]
for ns in (1, 3):
    h, mem = cnt(e0, ns)
    dec = e0*np.exp(rng.uniform(-1.5, 1.5, 20000)); luck = np.mean([cnt(u, ns)[0] for u in dec])
    pr(f"    within {ns} sigma: {h} member(s) {[cands[v] + ' = %.4e' % v for v in mem]};  expected by luck for a random eps of this size: {luck:.2f}")
for name, v in [("2 alpha^3", 2*al**3), ("(alpha/pi)^2 / 7", (al/math.pi)**2/7), ("alpha^2/(8 pi^2)", al**2/(8*math.pi**2)), ("alpha^2/(22 pi)", al**2/(22*math.pi))]:
    pr(f"      named near-misses: {name:18s} = {v:.4e}: off {100*(v/e0-1):+.2f}% = {(v-e0)/se:+.1f} sigma")
pr(f"    The family has {len(cv)} members; any eps of this size meets some member by luck. A hit here would not be evidence, and the simplest-looking")
pr("    one (2 alpha^3) misses by 6 sigma. FAIL: no correction is picked out.")

# ======================================================================================================
head("SECTION 5. ROUTE 4: A SECOND BALANCE ('the link sits at a balance point')")
pr("    Rule 7 is one balance (self weight = relationship weight in the ROOT-mass matrix a I + b P + b* P^2, |b| = a/sqrt2). The twist is the link's angle.")
pr("    A second balance would be one more equation, fixing that angle. Every natural one is tried:")
def bal_mass_matrix(d):      # relationship/self weight of the MASS matrix (the square of the root matrix)
    return 6*abs(r2*mp.e**(1j*d) + mp.mpf(1)/2*mp.e**(-2j*d))**2/12
rows = [
 ("B1 mass matrix itself balanced (its links weigh as much as its diagonal)", "cos 3d = -1/(4 sqrt2)", mp.acos(-1/(4*r2))/3),
 ("B2 link half real, half imaginary (push = pull in the link)", "delta = pi/4", pi/4),
 ("B3 lightest member's self exactly cancelled by its links", "delta = pi/12", pi/12),
 ("B4 link purely real (no twist; two members equal)", "delta = 0", mp.mpf(0)),
 ("B7 middle member sits at the family average (x_mu = 1)", "cos(d + 4pi/3) = 0", pi/6),
 ("B8 geometric balance: muon is the mean of electron and tau (x_mu^2 = x_e x_tau)", "", mp.findroot(lambda d: fam(d)[1]**2 - fam(d)[0]*fam(d)[2], mp.mpf("0.1"))),
 ("B9 three-way product equals the two-way sum (det = e2 share): e3 = 1/18", "cos 3d = (1/2+1/18) sqrt2", mp.acos((mp.mpf(1)/2+mp.mpf(1)/18)*r2)/3),
]
for nm, cond, dv in rows:
    e_, t_ = ratios(dv)
    pr(f"      {nm:82s} twist {f(dv,7):>10s}  ({f((dv/d45-1)*100,3):>7s}% off)  -> m_e/m_mu x 1e3 = {f(e_*1000,5)} (measured {f(emu*1000,5)})")
pr(f"    measured value of B1's ratio (mass-matrix relationship/self): {f(bal_mass_matrix(d45),8)};  the tau's share of the total mass = {f(fam(d45)[2]**2/6,8)}")
b9 = rows[-1][2]
pr(f"    The closest (B9) is off by {f((b9/d45-1)*100,3)}% in the twist: {f((b9-d_meas)/SIG2,3)} sigma with nothing assumed, {f((b9-d45)/SIG1,3)} sigma with 45 exact, and it was picked after looking.")
pr("    A balance with virial splitting was already ruled out (39.3: Q <= 0.5).")
pr("    FAIL: no second balance of the three placements gives the twist. 'The electron stops short of massless' (B3) is the right picture and gives no number.")
xe2 = float(fam(d45)[0]**2)
pr(f"\n    A curiosity, recorded so it is not rediscovered: m_e/M0 = {xe2:.6e} = {xe2/al:.5f} x alpha, and 2/9 = 0.22222 (off {100*(xe2/al/(2/9)-1):+.2f}%).")
pr(f"    With M0 from mu and tau that would 'predict' m_e = {float(2/9*alpha*M0):.5f} MeV against 0.51100: off {100*(float(2/9*alpha*M0/me)-1):+.2f}%, i.e. millions of sigma. FAIL.")

# ======================================================================================================
head("SECTION 6. ROUTE 5: A GRAVITY TIE (the electron mass from the Planck mass)")
L = mp.log(MPl/me); Lr = mp.log(MPl_red/me)
pr(f"6a. m_e/M_Planck = {f(me/MPl,6)};  ln(M_Planck/m_e) = {f(L,9)} +- {f(sMPl_rel,2)};  with the reduced Planck mass: {f(Lr,9)}")
pr(f"    As m_e = M_Planck exp(-c/alpha): c = {f(L*alpha,8)} (reduced: {f(Lr*alpha,8)}).   3/8 = 0.375 gives m_e = {f(MPl*mp.e**(-mp.mpf(3)/8/alpha),5)} MeV ({f((MPl*mp.e**(-mp.mpf(3)/8/alpha)/me-1)*100,3)}% off).")
# family
aMZ = 1/127.951
facs = {"1": 1.0, "pi": math.pi, "1/pi": 1/math.pi, "2pi": 2*math.pi, "pi/2": math.pi/2, "4pi": 4*math.pi, "pi^2": math.pi**2, "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3)}
G = {}
for (mn, Mv) in (("M_Pl", float(MPl)), ("reduced M_Pl", float(MPl_red))):
    for (an, av) in (("alpha(0)", al), ("alpha(M_Z)", aMZ)):
        for fn_, fv in facs.items():
            for p in range(1, 17):
                for q in range(1, 17):
                    if math.gcd(p, q) != 1: continue
                    for n in range(-2, 3):
                        lnm = math.log(Mv) + n*math.log(av) - p/q*fv/av
                        G[lnm] = f"{mn} x {an}^{n} x exp(-({p}/{q}) {fn_}/{an})"
gv = np.array(sorted(G)); target = math.log(float(me))
def gcount(t, tol): lo, hi = np.searchsorted(gv, [t - tol, t + tol]); return hi - lo, gv[lo:hi]
pr(f"    family: M x alpha^n x exp(-k f/alpha), M in (Planck, reduced), alpha at 0 or at the Z, k = p/q <= 16/16, f in {list(facs)}, n = -2..2: {len(gv)} formulas")
for tolname, tol in (("the Planck mass's own error (11 ppm)", 1.1e-5), ("1 per cent", 0.01), ("10 per cent", 0.1)):
    h, mem = gcount(target, tol)
    decoy = target + rng.uniform(-3, 3, 20000); luck = np.mean([gcount(u, tol)[0] for u in decoy])
    ex = "; ".join(G[v] + f" ({100*(math.exp(v-target)-1):+.3f}%)" for v in mem[:3])
    pr(f"      within {tolname:38s}: {h:3d} formula(s); expected by luck {luck:7.3f}.  {ex}")
pr("    Reading: formulas of this kind land within 1% about as often as luck says, and NONE lands within the Planck mass's error. FAIL / not evidence.")

pr("\n6b. The best any gravity formula could ever do.  The Planck mass is known to 11 ppm (G to 22 ppm, CODATA 2022).")
dd_g = abs(sens_me)*sMPl_rel
pr(f"    An exact formula m_e = (number) x M_Planck would fix the electron mass to 11 ppm, hence the twist to +-{f(dd_g,3)} rad = +-{f(dd_g/d45/ppm,3)} ppm.")
pr(f"    The twist is MEASURED to +-{f(SIG1/d45*1e9,3)} parts per billion (45 exact). So measurement is {f(dd_g/SIG1,4)} times sharper than any gravity tie can be tested.")
pr(f"    (The 2/9 shortfall, 0.79 ppm, is {f(eps*d45/dd_g*1,3)} of that gravity error bar: a gravity tie could not even see it.)")
pr("    Also: the twist is m_e relative to the LEPTON SCALE M0 (313.86 MeV), not m_e alone. A gravity tie for m_e fixes the twist only if a second rule")
pr("    gives M0 (equivalently m_tau). So the gravity route needs TWO rules, not one.")

pr("\n6c. The electron-gravity balance of an earlier test read backwards (far-side mass M_far = sqrt(m_e M_Planck), so m_e = M_far^2/M_Planck).")
tab = [(1e7, 0.23557), (3e7, 0.23401), (1e8, 0.23229), (1e9, 0.22902)]          # the earlier two-loop table (electron_gravity_balance_2loop .txt)
slope = (tab[3][1] - tab[2][1])/(math.log(tab[3][0]) - math.log(tab[2][0]))     # d sin^2 / d ln M
Mclose = 2.127e8*1e3  # MeV
me_back = Mclose**2/float(MPl)
pr(f"    The far-side mass that closes the strengths exactly (two loops): 2.127e8 GeV  ->  m_e = {me_back:.3f} MeV: {me_back/float(me):.1f} times the real electron mass.")
pr(f"    d sin^2/d ln(M_far) = {slope:.5f}. The +-0.3% left in the meeting itself is +-{0.003*0.23122/abs(slope):.2f} in ln M_far, i.e. a factor {math.exp(2*0.003*0.23122/abs(slope)):.1f} either way in m_e.")
me_need = 1e-6*float(d45/abs(sens_me))          # relative accuracy on m_e that gives the twist to 1 ppm
need_s = me_need*abs(slope)/2/0.23122
pr(f"    To give the twist to 1 ppm the electron mass is needed to about {me_need*1e6:.0f} ppm, so sin^2 would have to be CALCULATED to {need_s:.1e} relative;")
pr(f"    it is measured to 1.7e-4 and calculated to about 3e-3. Short by a factor of {3e-3/need_s:.0e}.")
pr("    FAIL as a derivation of the twist. If item 1 adopts the balance, it would explain why the electron is light to within a factor of a few, no better.")
pr(f"    In twist language: a factor 7.3 in m_e is a twist of {f(mp.findroot(lambda d: ratios(d)[0]-emu*me_back/float(me), (mp.mpf('0.05'), mp.mpf('0.25')), solver='bisect'),5)} instead of {f(d45,5)}.")

# ======================================================================================================
head("SECTION 7. ROUTE 6: EARLIER IDEAS, RE-CHECKED AT TODAY'S PRECISION")
ring = sum((1 - mp.cos(2*pi*n/3))/n**2 for n in range(1, 200001))/pi**2
pr(f"    M1 ring of three seats, two paths: (1/3)(2/3) = 2/9; wave sum (200,000 terms) = {f(ring,7)}. As a fraction of a turn it is {f(mp.mpf(2)/9*2*pi,6)} rad, not 0.2222.")
pr(f"       Needs 'one ring = one radian' (no rule gives it) and then still misses by {NS} sigma with 45 exact. WITHDRAWN (34.5); unchanged.")
pr(f"    M2 twist = Koide value / 3: with 45 exact this is 2/9 again ({NS} sigma). As a free-lean rule it is variant 1e: one free number, tau {f(tauq,8)} MeV, no test.")
pr(f"    M3 twist ladder 2/9 x (1, 1/2, 1/3) for leptons, down, up: dropped in an earlier test (predicts m_u 2.83 vs 2.13). Unchanged.")
pr(f"    M4 twist as polarity: electron and positron masses agree to 8 parts per billion, so the twist is the same for both. FAIL (72.5). Unchanged.")
pr(f"    M5 'a limit approached and never reached': with 45 exact the electron is massless at pi/12 = 15 degrees and sits at {f(d45*180/pi,6)} degrees.")
pr(f"       HOLDS as the picture (Landau, Abrikosov & Khalatnikov 1954: a charged particle cannot be massless). It gives the reason, not the number.")
pr(f"    M6 the two ends (top link about 1, electron link about 3e-6) under one rule: tested as 6c above. FAIL as a rule; hint only.")
pr(f"    M7 mass as the balance point: tested as Section 5. FAIL for the twist; rule 7 is still the balance and still holds:")
pr(f"       self/relationship = {f(2/r_meas**2,8)} +- {f(abs(4*dr/r_meas**3),2)}.")

# ======================================================================================================
head("SECTION 8. SUMMARY TABLE AND THE SHARP PREDICTION")
pr(f"    {'route':58s} {'result':60s} label")
tabl = [
 ("1 QED correction, masses at one common energy", f"shift {f(shift1/ppm,5)} ppm vs +{f(need/ppm,3)} needed (x{f(shift1/need,4)})", "FAIL"),
 ("1 QED correction, each mass at its own scale (2 loops)", f"{f(sh_emu/ppm,4)} ppm vs +{f(need/ppm,3)} needed: wrong sign", "FAIL"),
 ("1 a two-loop-size log with a free coefficient", f"needs c = {f(c_need/a_pi**2,3)} (alpha/pi)^2; one number for one datum", "not a prediction"),
 ("2 octonion / G2 / F4 / E6 phases", "Lindemann: no algebraic phase is a rational radian; no turn-fraction fits", "FAIL"),
 ("3 three-placement circulant invariants", f"{tot1} simple values hit in {len(famA)} tried x 15 invariants", "FAIL"),
 ("4 a second balance", "nearest natural balance 0.05% off (18 sigma), others far", "FAIL"),
 ("5 gravity formulas M_Pl exp(-k/alpha)", "hits no better than luck; none within G's error", "FAIL"),
 ("5 electron-gravity balance backwards", f"m_e within a factor {me_back/float(me):.0f}; can never reach ppm", "FAIL (hint only)"),
 ("6 The ring, Q/3, ladder, polarity", "as recorded in earlier tests", "FAIL / withdrawn"),
 ("6 'approached, never reached'", "right picture, no number", "HOLDS as picture"),
 ("- nearness to 2/9", f"luck about 1 in 10,000; shortfall 0.79 ppm is real ({NS} sigma)", "real, unexplained"),
]
for a, b, c in tabl: pr(f"    {a:58s} {b:60s} {c}")
pr(f"\n    THE TWIST, as a measured input (pole masses, 45 exact):  delta = {f(d45,12)} +- {f(SIG1,2)} rad  = 2/9 x (1 - {f(eps,5)})")
pr(f"    SHARP PREDICTION (from rule 7 alone, twist read from e/mu):  m_tau = {f(tau45,8)} +- {f(tau45_err,1)} MeV.   PDG 2024: {f(mtau,6)} +- {f(smtau,2)} ({f((tau45-mtau)/smtau,3)} sigma).")
pr(f"    Rule 7 is broken if the tau mass settles outside {f(tau45,7)} by more than its error; e.g. a final 1776.93 +- 0.01 would be {f((tau45-mtau)/mp.mpf('0.01'),3)} sigma.")
pr(f"    The twist itself has no tau test: 45-exact, 2/9-exact and the log-corrected version differ by {f(abs(tau29-tau45)*1000,2)} to {f(abs(tau_c-tau45)*1000,2)} keV.")
OUT.close()
