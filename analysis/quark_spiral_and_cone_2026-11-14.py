"""The quark masses as another spiral, or a cone. Both pictures are written as exact statements FIRST (section 0) and only then tested, each at one common energy with error bars.
  cone   = the three rungs of a ladder sit on a circle around one axis: fixed opening, the rungs 120 degrees apart, one twist (Koide 1982; Foot 1994).
  spiral = each turn also climbs: every rung multiplies the one before by a set factor (a ladder of equal steps), possibly with a pitch tied to the tilt.
Masses: lattice ratios (FLAG 2024; HPQCD; Fermilab-MILC-TUMQCD), PDG 2024, two-loop running to the top's own energy, as in an earlier test.
Credit: Koide 1982; Foot 1994; Brannen 2006; Rivero & Gsponer 2005; Froggatt & Nielsen 1979; Wolfenstein 1983; Gatto, Sartori & Tonin 1968;
Fritzsch 1977; Harari, Haut & Weyers 1978; Rodejohann & Zhang 2011; Kartavtsev 2011; Cabibbo 1963; Kobayashi & Maskawa 1973;
Gross & Wilczek 1973; Politzer 1973; Tarasov, Vermaseren & Larin 1980; FLAG 2024; PDG 2024."""
import math, numpy as np
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)

MZ = 91.1876
def b(nf): return (11 - 2*nf/3)/4, (102 - 38*nf/3)/16
def g(nf): return 1.0, (202/3 - 20*nf/9)/16
def run(a0, mu0, mu1, th, lnm=0.0):
    pts = sorted([mu0, mu1] + [t for t in th if min(mu0, mu1) < t < max(mu0, mu1)], reverse=(mu1 < mu0)); a_, l = a0, lnm
    for x, y in zip(pts, pts[1:]):
        mid = math.sqrt(x*y); nf = 6 if mid > th[2] else 5 if mid > th[1] else 4 if mid > th[0] else 3; b0, b1 = b(nf); g0, g1 = g(nf)
        sol = solve_ivp(lambda t, z: [-b0*z[0]**2 - b1*z[0]**3, -g0*z[0] - g1*z[0]**2], [math.log(x*x), math.log(y*y)], [a_, l], rtol=1e-9, atol=1e-11).y[:, -1]; a_, l = sol
    return a_, l
def cone(m):
    """Cone coordinates of three masses: size a^2 (MeV), opening (degrees from the axis), twist (degrees), and relationship/self = r^2/2."""
    x = np.sqrt(np.sort(np.asarray(m, float))); a = x.mean(); c = x/a - 1; r2 = (2/3)*(c*c).sum(); r = math.sqrt(r2)
    d = math.acos(max(-1.0, min(1.0, (4/3)*(c**3).sum()/r**3)))/3
    return a*a, math.degrees(math.atan(r/math.sqrt(2))), math.degrees(d), r2/2
def ms(v): v = np.asarray(v, float); return v.mean(), v.std()

head("0. The two pictures as exact statements, written before any number is looked at")
pr("  CONE.    sqrt(mass of rung k) = size x (1 + opening x cos(twist + k x 120 degrees)).  Three rungs, three numbers: size, opening, twist.")
pr("           C1  each ladder is a cone                       (any three masses are: no content by itself)")
pr("           C2  all ladders share the leptons' opening, 45 degrees")
pr("           C3  the two quark ladders share ONE opening (whatever it is)")
pr("           C4  all six quarks sit on one 45-degree cone   (Kartavtsev 2011)")
pr("           C5  both quark ladders are cones about the SAME axis -> then the weak force's tilt between them is exactly zero")
pr("  SPIRAL.  each turn climbs: rung k+1 = rung k x one set factor.")
pr("           S1  each ladder climbs by equal steps           (middle^2 = bottom x top)")
pr("           S2  the up ladder is the down ladder with every step squared (one spiral, two pitches)")
pr("           S3  the step is set by the tilt: down ladder steps = leak^2, up ladder steps = leak^4   (Froggatt & Nielsen 1979; Wolfenstein 1983)")
pr("           S4  the two ladders are ONE spiral seen from its two ends: same opening, twist mirrored (the spiral of an earlier test)")
pr("           S5  the leptons' exact cone, wound into a spiral: same opening and twist, plus one climb factor per rung")
pr("  9 statements. A pass needs the measured value inside its own error bar, at one energy, for BOTH ladders.")

head("1. The masses at one energy (the top's own, 162.5 GeV), with errors")
rng = np.random.default_rng(47); N = 600; D, Uq, six = [], [], []
for _ in range(N):
    sud = rng.normal(27.23, 0.10); ud = rng.normal(0.465, 0.024); cs = rng.normal(11.766, 0.030); bc = rng.normal(4.577, 0.008)
    aZ = rng.normal(0.1180, 0.0009)/math.pi; mc = rng.normal(1.273, 0.0046); mt = rng.normal(162.5, 0.8); th = [1.273, 4.183, 162.5]
    mcT = mc*math.exp(run(run(aZ, MZ, mc, th)[0], mc, mt, th)[1])*1e3; mtT = mt*1e3          # MeV at the top's energy
    msT = mcT/cs; mdT = msT/(sud*(1 + ud)/2); muT = mdT*ud; mbT = bc*mcT
    D.append((mdT, msT, mbT)); Uq.append((muT, mcT, mtT))
D, Uq = np.array(D), np.array(Uq); L = np.array([0.51099895, 105.6583755, 1776.93])
Ls = np.column_stack([np.full(4000, L[0]), np.full(4000, L[1]), rng.normal(1776.93, 0.09, 4000)])
for nme, A in (("down, strange, bottom", D), ("up, charm, top", Uq)):
    pr(f"  {nme:<22}: " + ", ".join(f"{A[:, i].mean():.4g} +/- {A[:, i].std():.2g}" for i in range(3)) + " MeV")
pr("  electron, muon, tau   : 0.5110, 105.66, 1776.93 +/- 0.09 MeV (their own masses; these barely move with energy)")
pr("  Steps up each ladder (second rung / first, third / second):")
for nme, A in (("leptons", Ls), ("down ladder", D), ("up ladder", Uq)):
    s1, s2 = A[:, 1]/A[:, 0], A[:, 2]/A[:, 1]; pr(f"     {nme:<12} {s1.mean():8.1f} +/- {s1.std():.1f}    {s2.mean():8.1f} +/- {s2.std():.1f}")

head("2. The cone: where each ladder sits (size, opening, twist)")
co = {}
for nme, A in (("leptons", Ls), ("down ladder", D), ("up ladder", Uq)):
    C = np.array([cone(r_) for r_ in A]); co[nme] = C
    pr(f"  {nme:<12} size {C[:, 0].mean():9.2f} +/- {C[:, 0].std():.2f} MeV   opening {C[:, 1].mean():.3f} +/- {C[:, 1].std():.3f} deg   twist {C[:, 2].mean():.3f} +/- {C[:, 2].std():.3f} deg   relationship/self {C[:, 3].mean():.4f} +/- {C[:, 3].std():.4f}")
pr("  C1  each ladder is a cone: TRUE for any three masses, so it says nothing by itself. The content is in the three numbers above.")
oL, oD, oU = co["leptons"][:, 1], co["down ladder"][:, 1], co["up ladder"][:, 1]
pr(f"  C2  all share the leptons' 45 degrees: leptons {oL.mean():.4f} +/- {oL.std():.4f} (holds). Down {oD.mean():.2f} ({(oD.mean() - 45)/oD.std():+.0f} sigma). Up {oU.mean():.2f} ({(oU.mean() - 45)/oU.std():+.0f} sigma). FAIL.")
pr(f"  C3  the two quark ladders share one opening: {oD.mean():.2f} against {oU.mean():.2f}: apart by {(oU.mean() - oD.mean())/math.hypot(oD.std(), oU.std()):.0f} sigma. FAIL.")
allsix = np.concatenate([D, Uq], axis=1); Q6 = allsix.sum(1)/np.sqrt(allsix).sum(1)**2
pr(f"  C4  all six quarks on one 45-degree cone: sum of masses / (sum of square roots)^2 = {Q6.mean():.4f} +/- {Q6.std():.4f}; a 45-degree cone needs 0.6667: {(Q6.mean() - 2/3)/Q6.std():+.0f} sigma. FAIL.")
pr("  C5  both ladders are cones about the SAME axis: then one turn of 120 degrees maps each ladder to itself, the two mass patterns can be")
pr("      sorted by the same three directions, and the weak force changes up-type into down-type with NO leak between families.")
pr("      (an earlier test found the same.) Measured leaks: 0.2243, 0.041, 0.0037. Zero is the right first guess (the leaks are small) and FAILS strictly")
pr("      (the first by more than 100 sigma).")
pr("      So the tilt between the ladders IS the amount by which the two cones do not share an axis. A cone alone cannot give it.")
tw = {k_: v_[:, 2] for k_, v_ in co.items()}
pr(f"  Seen in the numbers, going leptons -> down -> up: the opening widens ({oL.mean():.1f}, {oD.mean():.1f}, {oU.mean():.1f} deg) and the twist unwinds")
pr(f"      ({tw['leptons'].mean():.2f}, {tw['down ladder'].mean():.2f}, {tw['up ladder'].mean():.2f} deg). Two orderings of three things line up like this by chance about 1 time in 3. A trend, not a law.")
pr("  -> The cone is exact for the leptons and for nothing else. For quarks it is a way of writing three masses, with openings that match no simple value")
pr("     (an earlier test: 5/4 and 5/3 fail) and that differ between the ladders.")

head("3. The spiral: does each turn climb by a set factor?")
for nme, A in (("leptons", Ls), ("down ladder", D), ("up ladder", Uq)):
    e = A[:, 1]**2/(A[:, 0]*A[:, 2]); l = np.log(A[:, 2]/A[:, 1])/np.log(A[:, 1]/A[:, 0])
    pr(f"  S1  {nme:<12} middle^2 / (bottom x top) = {e.mean():8.3f} +/- {e.std():.3f}  (equal steps need 1)     second step / first step, in powers: {l.mean():.3f} +/- {l.std():.3f}  (need 1)")
pr("      FAIL for all three as an exact rule. But look at the last column: the quark ladders are 1.3 and 0.9, the leptons 0.5.")
pr("      The quarks are much NEARER a spiral of equal steps than the leptons are. The leptons are the cone; the quarks lean toward the spiral.")
p1 = np.log(Uq[:, 1]/Uq[:, 0])/np.log(D[:, 1]/D[:, 0]); p2 = np.log(Uq[:, 2]/Uq[:, 1])/np.log(D[:, 2]/D[:, 1])
pr(f"  S2  up ladder = down ladder with every step squared: power needed on the first step {p1.mean():.3f} +/- {p1.std():.3f}, on the second {p2.mean():.3f} +/- {p2.std():.3f}; need 2 and 2.")
pr(f"      One power for both steps would also do: they differ by {(p1.mean() - p2.mean())/math.hypot(p1.std(), p2.std()):.0f} sigma. FAIL.")
lam = 0.2243
pr("  S3  the step is set by the tilt (leak = 0.2243). Steps DOWN each ladder, measured against the rule:")
for nme, r_, pw in (("down / strange", D[:, 0]/D[:, 1], 2), ("strange / bottom", D[:, 1]/D[:, 2], 2), ("up / charm", Uq[:, 0]/Uq[:, 1], 4), ("charm / top", Uq[:, 1]/Uq[:, 2], 4)):
    pr(f"         {nme:<17} {r_.mean():.5f} +/- {r_.std():.5f}    leak^{pw} = {lam**pw:.5f}    measured / rule = {r_.mean()/lam**pw:.2f}    power of the leak that fits: {math.log(r_.mean())/math.log(lam):.2f}")
pr("      One of four lands (down / strange, the Cabibbo link that already passes: Gatto, Sartori & Tonin 1968). The other three are off by")
pr("      factors of 0.4 to 1.5. The fitted powers (2.0, 2.7, 4.2, 3.7) are near whole numbers and not on them. An old, real, rough pattern. As a law it FAILS.")
pr(f"  S4  one spiral seen from its two ends: needs the same opening (C3: FAIL) and a mirrored twist. Twists: down {tw['down ladder'].mean():.2f} +/- {tw['down ladder'].std():.2f}, up {tw['up ladder'].mean():.2f} +/- {tw['up ladder'].std():.2f} deg.")
pr(f"      Mirror through the axis would need them to add to 60: they add to {tw['down ladder'].mean() + tw['up ladder'].mean():.1f}. Mirror in the turning sense leaves a cone's three masses")
pr("      unchanged, so the two ladders would be EQUAL (up = down, charm = strange, top = bottom). FAIL both ways.")
pr("      What survives of the earlier spiral here: the two ladders ARE the two ends of one state (the quark pair joins the up partner one way and")
pr("      the down partner the other way, 45.3). That fixes which links exist. It does not fix their sizes.")
cl = np.sqrt(L/L.mean()); x = np.sqrt(L); cl = x/x.mean()
for nme, A in (("down ladder", D), ("up ladder", Uq)):
    r1 = np.sqrt(A[:, 1]/A[:, 0])/(cl[1]/cl[0]); r2 = np.sqrt(A[:, 2]/A[:, 1])/(cl[2]/cl[1])
    pr(f"  S5  leptons' cone wound into a spiral, {nme}: climb factor needed on the first step {r1.mean():.3f} +/- {r1.std():.3f}, on the second {r2.mean():.3f} +/- {r2.std():.3f}; must be equal. FAIL.")

head("4. Attacks on my own tests")
pr("  (a) Did the choice of energy decide the outcome? The down ladder's cone and its steps do not change with energy under the strong force (all")
pr("      three masses run alike), so C2, C3, S1 and S3 for the down ladder are the same at every energy below the top. The up ladder was taken")
pr("      at the top's own energy; an earlier test followed both to 10^16 GeV and they drift in opposite directions. No energy rescues C2, C3 or C4.")
pr("  (b) Was a signed root missed? With the lightest root taken negative the opening only widens, so C2 to C4 still fail:")
for nme, A in (("leptons", Ls), ("down ladder", D), ("up ladder", Uq)):
    xs = np.sqrt(A); xs[:, 0] *= -1; Qn = (xs**2).sum(1)/xs.sum(1)**2; pr(f"      {nme}: with a negative lightest root, opening {math.degrees(math.atan(math.sqrt(3*Qn.mean() - 1))):.2f} deg.")
pr("  (c) Nine statements were tested. One sub-test of S3 landed, and it was already known (the Cabibbo link). Nothing new passed, so there is")
pr("      no question of a lucky hit to discount.")

head("5. Summary")
pr("  CONE   : exact for the leptons (45.000 degrees, a part in a million). For quarks it is only a way of writing three masses. FAIL as a rule.")
pr("  SPIRAL : the right rough shape for quarks (each rung climbs by a large factor; quarks are nearer equal steps than leptons; the steps are")
pr("           roughly whole powers of the leak). FAILS as an exact rule in every form tested.")
pr("  LEARNED: (1) leptons are the cone, quarks lean to the spiral: the two are different shapes, which is why the lepton rule never carried over.")
pr("           (2) The tilt between the ladders is exactly the failure of the two cones to share an axis. (3) The down ladder's shape does not")
pr("           depend on the energy at which it is looked at, so it is a definite object to explain: opening 48.15 degrees.")
pr("  NOT MOVED: bottom, charm, strange and top stay measured. Scorecard unchanged.")
