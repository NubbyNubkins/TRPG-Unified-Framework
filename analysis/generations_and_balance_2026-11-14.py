"""Routes B and C (14 Nov 2026; earlier tests).
B. Can three generations sit beside the weak force of an earlier test without holes?
C. The lepton 45 degree rule as a balance: is there a derivation? And do the quark families show a balance of their own, tested strictly?
Uses the matrix tools of forces_from_entanglement_2026-11-14.py and the 2-loop running of quark_common_scale_2026-11-14.py (same folder).
Credit: Jordan, von Neumann & Wigner 1934; Kobayashi & Maskawa 1973; Cabibbo 1963; Koide 1981, 1990; Foot 1994; Brannen 2006; Rivero & Gsponer 2005;
Boyle 2020; Dubois-Violette 2016; Furey 2014, 2018; Gillard & Gresnigt 2019; Manogue & Dray 1999; Tarrach 1981; Chetyrkin 1997; Vermaseren, Larin & van Ritbergen 1997;
FLAG 2024 (lattice quark-mass ratios); HPQCD; Fermilab-MILC-TUMQCD; PDG 2024; ATLAS & CMS (Higgs rates); LEP (Z width); Eberhardt et al. 2012."""
import math, os, numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "forces_from_entanglement_2026-11-14.py")).read()
exec(src[src.index("import math, itertools"):src.index('head("1. Entanglement')])
exec(src[src.index("def jordan_algebra(n):"):src.index("names = {2:")])

head("B1. Three placements of the same relationship matrix (three seats), instead of three slots of one matrix")
n = 8; dim, T = jordan_algebra(n); der = derivations(dim, T)
W = [0, 1, 2] + [3 + b*n + c for b in range(3) for c in (0, 1)]; notW = [i for i in range(dim) if i not in W]
face = restrict(der, lambda g: g[np.ix_(notW, W)].ravel()); both = restrict(face, lambda g: g[:, 2])
def three(a, b, c):
    Z = np.zeros((3*dim, 3*dim)); Z[:dim, :dim] = a; Z[dim:2*dim, dim:2*dim] = b; Z[2*dim:, 2*dim:] = c; return Z
O = np.zeros((dim, dim))
allturns = [three(g, O, O) for g in der] + [three(O, g, O) for g in der] + [three(O, O, g) for g in der]
P = np.zeros((3*dim, 3*dim))
for s in range(3): P[((s+1) % 3)*dim:((s+1) % 3 + 1)*dim, s*dim:(s+1)*dim] = np.eye(dim)       # move every seat one place round
blind = restrict(allturns, lambda g: (g@P - P@g).ravel())
seen = [three(g, O, O) for g in both] + [three(O, g, O) for g in both] + [three(O, O, g) for g in both]
blind_seen = restrict(seen, lambda g: (g@P - P@g).ravel())
pr(f"  Three copies of the 3x3 octonion matrix, one per seat. No number systems are joined, and nothing squares to nothing (a sum of three such matrices).")
pr(f"  All turns: {len(allturns)}. Turns that cannot tell the seats apart (commute with moving every seat one place round): {len(blind)} (rank {rank_of(blind)}).")
pr(f"  With the observer's spot and facing direction: {len(seen)} carriers if each seat had its own forces; {len(blind_seen)} (rank {rank_of(blind_seen)}) if the forces are blind to the seat.")
pr(f"  Link states: 3 x 16 = {3*16} components = three generations of weak pairs. One weak force, the same for all three.")
pr("  -> A hole-free arrangement EXISTS: the same shape in three placements (the earlier own wording), with forces that are blind to the placement.")
pr("     The seat space carries the 3x3 mass matrices of rule 5; the octonion matrix carries the forces. They are two different matrices.")
pr("  Measured: the weak force is the same for all three generations to 0.1-0.2% (tau, muon and pion decays); no force that tells generations apart is seen.")

head("B2. What that costs, and what can replace it")
pr("  'At most 3 generations because octonion matrices stop at 3x3' used the octonion matrix's slots as the generations. With seats")
pr("  instead of slots, that argument no longer applies: nothing limits the number of seats.")
pr("  What still fixes 3:")
pr("   - at least 3: an irremovable twist (matter-antimatter difference) needs 3 (Kobayashi & Maskawa 1973; an earlier test). A theorem.")
pr("   - at most 3: only rule 6 ('nothing extra'), which holds again in the earlier picture: the minimum the twist needs, and no more. A rule, not a theorem.")
pr("  Data: a fourth generation of the same kind is excluded (Higgs production would be about 9 times larger; Z width gives 2.9963 +/- 0.0074 light neutrinos).")
pr("  Other published routes to 3 inside one structure, checked against rule 1:")
pr("   - doubling the octonion matrix with complex numbers (Boyle 2020): complex x octonions has holes.")
pr("   - the next number system, 16 parts (sedenions; Gillard & Gresnigt 2019): has holes (it fails rule 1 by Hurwitz's theorem).")
pr("   - the three four-part systems containing the facing direction (Manogue & Dray 1999): that '3' is already the three colours in TRPG.")
pr("  -> None passes rule 1. FAIL for 'one structure gives both the forces and exactly three generations'.")

head("C1. The lepton 45 degree rule as a balance of two pushes")
m = np.array([0.51099895, 105.6583755, 1776.93]); s = np.sqrt(m); a = s.mean(); rho = ((s - a)**2).sum()/(3*a*a)
pr(f"  relationship weight / self weight = {rho:.6f}. In the tractor picture: at balance the two pushes are EQUAL (Newton's third law). That is the 45 degree rule.")
pr("  A published energy function gives it (Koide 1990): energy = (terms in the total) - k x (self weight) x (relationship weight). With the total kept")
pr("  (rule 2), the product is largest when the two weights are equal. It works BECAUSE the function treats self and relationship alike. That is put in.")
pr("  Wave test: three seats on a ring joined only by their links (no outside anchor: the full glass), with a twist phi through the ring:")
best = 0
for phi in np.linspace(0, 2*math.pi, 2001):
    w = 2*np.abs(np.sin(phi/2 + math.pi*np.arange(3)/3)); Q = (w**2).sum()/w.sum()**2 if w.sum() > 0 else 0; best = max(best, Q)
pr(f"   each seat's energy IS half own motion, half links (the virial theorem), but the three frequencies give a ratio of at most Q = {best:.3f}; leptons have 0.667.")
pr("  -> The balance reading restates the rule; the simplest wave model with that balance FAILS to give the lepton masses. The 45 degree rule stays measured.")

head("C2. The twist and the electron's mass are the same number")
M0 = a; d0 = math.pi/12
r = math.sqrt(2*rho); k = int(np.argmin(s))
delta = None
from scipy.optimize import brentq
f = lambda d: sorted(M0*(1 + r*np.cos(d + 2*math.pi*np.arange(3)/3)))[0] - s.min()
delta = brentq(f, 0.15, d0)
pr(f"  With exact 45 degree balance, the lightest member is massless when the twist is exactly pi/12 = 15 degrees = {d0:.6f}.")
pr(f"  Measured twist: {delta:.6f} ({math.degrees(delta):.3f} degrees). The gap, {math.degrees(d0 - delta):.3f} degrees, IS the electron's mass: sqrt(m_e)/average = {s.min()/M0:.5f}; gap in radians = {d0 - delta:.5f}.")
pr("  -> 'Why is the twist 2/9?' and 'why is the electron so light but not massless?' are one question, not two. It does not answer either.")

head("C3. Do the quark families have a balance of their own? Tested at one common energy, with errors")
MZ = 91.1876
def b(nf): return (11 - 2*nf/3)/4, (102 - 38*nf/3)/16
def g(nf): return 1.0, (202/3 - 20*nf/9)/16
def run(a0, mu0, mu1, th, lnm=0.0):
    pts = sorted([mu0, mu1] + [t for t in th if min(mu0, mu1) < t < max(mu0, mu1)], reverse=(mu1 < mu0)); a_, l = a0, lnm
    for x, y in zip(pts, pts[1:]):
        mid = math.sqrt(x*y); nf = 6 if mid > th[2] else 5 if mid > th[1] else 4 if mid > th[0] else 3; b0, b1 = b(nf); g0, g1 = g(nf)
        sol = solve_ivp(lambda t, z: [-b0*z[0]**2 - b1*z[0]**3, -g0*z[0] - g1*z[0]**2], [math.log(x*x), math.log(y*y)], [a_, l], rtol=1e-9, atol=1e-11).y[:, -1]; a_, l = sol
    return a_, l
def Q3(x): return sum(x)/sum(math.sqrt(v) for v in x)**2
rng = np.random.default_rng(39); N = 4000
sets = {"FLAG 2024, four sea quarks": dict(sud=(27.23, 0.10), ud=(0.465, 0.024)), "FLAG 2024, three sea quarks": dict(sud=(27.42, 0.12), ud=(0.485, 0.019))}
for nm, S in sets.items():
    Qd, Qu = [], []
    for _ in range(N):
        sud = rng.normal(*S["sud"]); ud = rng.normal(*S["ud"]); cs = rng.normal(11.766, 0.030); bc = rng.normal(4.577, 0.008)
        ms_md = sud*(1 + ud)/2; mb_ms = bc*cs
        Qd.append(Q3([1/(mb_ms*ms_md), 1/mb_ms, 1.0]))
        mc_mu = cs*sud*(1 + ud)/(2*ud)
        if len(Qu) < 600:                                  # running is slow; 600 samples are enough for the error
            aZ = rng.normal(0.1180, 0.0009)/math.pi; mc = rng.normal(1.273, 0.0046); mt = rng.normal(162.5, 0.8); th = [1.273, 4.183, 162.5]
            a_mc = run(aZ, MZ, mc, th)[0]; mc_at_mt = mc*math.exp(run(a_mc, mc, mt, th)[1])
            Qu.append(Q3([mc_at_mt/mc_mu, mc_at_mt, mt]))
    Qd, Qu = np.array(Qd), np.array(Qu)
    pr(f"  {nm}:")
    pr(f"    down family (d, s, b): Q = {Qd.mean():.5f} +/- {Qd.std():.5f}  -> relationship/self = {3*Qd.mean()-1:.4f} +/- {3*Qd.std():.4f};  against 3/4: {(Qd.mean()-0.75)/Qd.std():+.1f} sigma")
    pr(f"    up family   (u, c, t): Q = {Qu.mean():.5f} +/- {Qu.std():.5f}  -> relationship/self = {3*Qu.mean()-1:.4f} +/- {3*Qu.std():.4f};  against 8/9: {(Qu.mean()-8/9)/Qu.std():+.1f} sigma")
mtau = rng.normal(1776.93, 0.09, 20000); Ql = (0.51099895 + 105.6583755 + mtau)/(math.sqrt(0.51099895) + math.sqrt(105.6583755) + np.sqrt(mtau))**2
pr(f"  For comparison, leptons: Q = {Ql.mean():.6f} +/- {Ql.std():.6f} -> relationship/self = {3*Ql.mean()-1:.6f} +/- {3*Ql.std():.6f}; against 2/3: {(Ql.mean()-2/3)/Ql.std():+.1f} sigma.")
pr("  The down family's ratio does not change with energy under the strong force (all three masses run alike); the up family's is taken at the top's own scale.")
Qw = []
for _ in range(400):                                       # the same test with the wider top-mass error (PDG's value from production rates)
    sud = rng.normal(27.23, 0.10); ud = rng.normal(0.465, 0.024); cs = rng.normal(11.766, 0.030)
    aZ = rng.normal(0.1180, 0.0009)/math.pi; mc = rng.normal(1.273, 0.0046); mt = rng.normal(162.5, 1.8); th = [1.273, 4.183, 162.5]
    a_mc = run(aZ, MZ, mc, th)[0]; mc_at_mt = mc*math.exp(run(a_mc, mc, mt, th)[1]); Qw.append(Q3([mc_at_mt/(cs*sud*(1 + ud)/(2*ud)), mc_at_mt, mt]))
Qw = np.array(Qw)
pr(f"  Up family again with the wider top-mass error (162.5 +/- 1.8): Q = {Qw.mean():.5f} +/- {Qw.std():.5f}; against 8/9: {(Qw.mean()-8/9)/Qw.std():+.1f} sigma")
pr("  Attack: is there another energy at which both hold? Above the top's scale the top's own link pulls the ratios (one-loop running, Standard Model):")
def drift(tmax):
    def rhs(t, z):
        g1, g2, g3, yt, L = z; k = 1/(16*math.pi**2)
        return [k*41/10*g1**3, -k*19/6*g2**3, -k*7*g3**3, k*yt*(4.5*yt**2 - 17/20*g1**2 - 9/4*g2**2 - 8*g3**2), k*1.5*yt**2]
    return solve_ivp(rhs, [0, tmax], [0.462, 0.648, 1.166, 0.937, 0.0], rtol=1e-8, dense_output=True).sol
sol = drift(math.log(1e16/162.5)); S = sets["FLAG 2024, four sea quarks"]; sud, ud, cs, bc = S["sud"][0], S["ud"][0], 11.766, 4.577
aZ = 0.1180/math.pi; th = [1.273, 4.183, 162.5]; mc_mt = 1.273*math.exp(run(run(aZ, MZ, 1.273, th)[0], 1.273, 162.5, th)[1])
def QQ(logmu):
    L = sol(math.log(10**logmu/162.5))[4]; tc = 162.5/mc_mt*math.exp(L); bs = bc*cs*math.exp(-L); mcu = cs*sud*(1 + ud)/(2*ud); msd = sud*(1 + ud)/2
    return Q3([1/(tc*mcu), 1/tc, 1.0]), Q3([1/(bs*msd), 1/bs, 1.0])
for lm in (2.21, 4, 6, 8, 10, 12, 14, 16): qu, qd = QQ(lm); pr(f"    at 10^{lm:<5} GeV: up family Q = {qu:.4f} ({(qu/(8/9)-1)*100:+.2f}% from 8/9), down family Q = {qd:.4f} ({(qd/0.75-1)*100:+.2f}% from 3/4)")
lx = brentq(lambda l: QQ(l)[0] - 8/9, 2.3, 16) if (QQ(2.3)[0] - 8/9)*(QQ(16)[0] - 8/9) < 0 else None
if lx: pr(f"    The up family reaches 8/9 near 10^{lx:.1f} GeV, where the down family is {(QQ(lx)[1]/0.75-1)*100:+.2f}% from 3/4. The two move in OPPOSITE directions.")
pr("    -> No single energy gives both. Choosing a different energy for each family would be a free choice, which is fitting, not deriving.")
pr("  -> At one common energy the families sit NEAR relationship/self = 1 (leptons), 5/4 (down) and 5/3 (up); self shares near 1/2, 4/9, 3/8.")
pr("     Strictly: the down family is 0.1-0.2% from 3/4 but 2 to 3 sigma low. The up family is 0.4% from 8/9 and 4 to 5 sigma low.")
pr("     Both FAIL as exact rules. There is also no derivation of 3/4 or 8/9, and simple fractions are easy to find near any number.")
pr("     NOT adopted. Only the leptons have a balance rule that survives its own error bar.")
