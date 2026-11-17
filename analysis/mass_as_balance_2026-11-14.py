"""Mass as a balance point between containment and outward push ('leakage').
Part A checks the idea where physics already uses it, at four scales. Part B applies it to the open problems and tries to break it.
Credit: Poincare 1906; von Laue 1911; Chodos, Jaffe, Johnson, Thorn & Weisskopf 1974 (bag); DeGrand, Jaffe, Johnson & Kiskis 1975;
Friedberg & Lee 1977; Theberge, Thomas & Miller 1980; Burkert, Elouadrhiri & Girod 2018; Weizsacker 1935; Bohr & Wheeler 1939; Meitner & Frisch 1939;
Eddington; Chandrasekhar; Englert & Brout 1964; Higgs 1964; Nambu & Jona-Lasinio 1961; Weisskopf 1939; Pendleton & Ross 1981; Hill 1981;
Chanowitz, Furman & Hinchliffe 1978; Bagger & Naculich 1991; Degrassi et al. 2012; Buttazzo et al. 2013; Koide 1981; Foot 1994; Brannen 2006;
HotQCD 2019; PDG 2024; AME 2020 (atomic masses)."""
import math, numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
hc = 197.3269804          # MeV fm

head("A1. The proton: mass as the balance between containment (B) and the quarks' outward push")
x0 = 2.04279
pr("  Energy of three quarks held in a region of radius R:  E(R) = (outward push) 3 x0 hc / R  +  (containment) (4 pi/3) B R^3")
pr(f"  {'B (MeV/fm^3)':>14s} {'balance R (fm)':>15s} {'E at balance':>13s} {'containment share':>18s} {'push share':>11s}")
for B in (30.0, 57.1, 60.1, 63.7, 120.0):
    E = lambda R: 3*x0*hc/R + 4*math.pi/3*B*R**3
    R = minimize_scalar(E, bounds=(0.2, 5), method="bounded", options=dict(xatol=1e-10)).x; cont = 4*math.pi/3*B*R**3
    pr(f"  {B:14.1f} {R:15.3f} {E(R):13.1f} {cont/E(R):18.4f} {1-cont/E(R):11.4f}")
pr("  -> At the balance point the split is ALWAYS 1/4 containment and 3/4 outward push, whatever the strength B: mass = 4 x containment energy.")
pr("     The split is set by how each side changes with size (R^3 against 1/R), not by how strong each side is. The mass grows only as B^(1/4).")
pr("  TRPG's record already does this in full (proton_sets_B: the containment balance), with the pion cloud as the leakage:")
pr("     B = 60.1 MeV/fm^3 (57.1-63.7), balance radius 0.995 fm, leakage (cloud) energy -220 MeV, proton 938.9 and Delta 1232.0 MeV reproduced.")
pr("     Skin: record step M1 gives a wall about 0.11 fm thick at a radius near 1 fm.")
pr("  Measured (Burkert, Elouadrhiri & Girod 2018): pressure inside the proton pushes OUT near the centre and IN beyond about 0.6 fm,")
pr("     and the two cancel exactly over the whole proton (von Laue's stability condition, 1911). First measurement; model-dependent.")

head("A2. 'When the outside stops pushing in, the centre bursts out': the temperature at which containment fails")
pr("  Balance of pressures: hot quark-gluon gas pressure (g pi^2/90) T^4 = B. Above that temperature nothing holds the skin.")
for B in (57.1, 60.1, 63.7):
    B4 = B*hc**3
    Ts = [(90*B4/(g*math.pi**2))**0.25 for g in (37.0, 34.0, 47.5)]
    pr(f"  B = {B:5.1f}: T = {Ts[0]:6.1f} MeV (2 quark kinds), {Ts[1]:6.1f} (pion gas subtracted), {Ts[2]:6.1f} (3 quark kinds)")
pr("  Measured (lattice QCD, HotQCD 2019; seen at RHIC and the LHC): 156.5 +/- 1.5 MeV, about 1.8e12 kelvin.")
pr(f"  -> The simple balance gets the right kind of answer (a trillion degrees) but is {100*(1-105.4/156.5):.0f}% low. As a strict prediction it FAILS.")
pr("     Known reason: the real change is gradual, not a wall giving way all at once.")

head("A3. The nucleus: fission as a shift in the balance, not a breaking of what is inside")
av, as_, ac, aa, ap = 15.75, 17.8, 0.711, 23.7, 11.18
def terms(A, Z):
    N = A - Z; d = 0.0
    if A % 2 == 0: d = ap/math.sqrt(A) if Z % 2 == 0 else -ap/math.sqrt(A)
    return np.array([av*A, -as_*A**(2/3), -ac*Z*(Z-1)/A**(1/3), -aa*(A-2*Z)**2/A, d])
before = terms(235, 92); after = terms(141, 56) + terms(92, 36); d = after - before
lab = ["volume (bulk holding)", "surface (the skin: containment)", "electric push (outward)", "neutron-proton balance", "pairing"]
pr("  U-235 + neutron -> Ba-141 + Kr-92 + 3 neutrons. Binding energy by part (liquid-drop formula), MeV:")
for l, b, a, dd in zip(lab, before, after, d): pr(f"    {l:34s} before {b:9.1f}   after {a:9.1f}   change {dd:+8.1f}")
Qm = 40.9218 + 8.0713 - (-79.733 - 68.769 + 3*8.0713)
pr(f"    total released: {d.sum():.0f} MeV by the formula; measured {Qm:.1f} MeV (from the measured masses). {100*(d.sum()/Qm-1):+.0f}%.")
pr(f"  Protons before and after: 92 and {56+36}. Neutrons: 144 and {141-56+92-36+3}. Nothing inside is destroyed.")
pr(f"  The skin gets BIGGER (costs {-d[1]:.0f} MeV) and the outward electric push drops by much more ({d[2]:.0f} MeV released).")
pr(f"  Energy released / total mass = {Qm/(236*931.494)*100:.3f}%.  -> The energy is the change in the balance. The nucleus does end up as two pieces.")
x = (ac*92**2/236**(1/3))/(2*as_*236**(2/3))
pr(f"  How close uranium sits to the tipping point: electric push / (2 x skin energy) = {x:.2f} (1.00 = no skin can hold it; Bohr & Wheeler 1939).")

head("A4. A star: the same balance at the largest scale")
G, Ms, Rs, mp, kB = 6.6743e-11, 1.98847e30, 6.957e8, 1.67262e-27, 1.380649e-23
Tbar = 0.5*G*Ms*0.61*mp/(Rs*kB)
pr(f"  Gravity pulls in, heat pushes out. Balance requires the Sun's average temperature to be about {Tbar/1e6:.1f} million K (centre: 15.7 million; surface: 5,800).")
pr("  Take away the outward push and it collapses; take away the containment and it flies apart (a supernova is the rebound).")
pr("  TRPG's neutron-star result (2.0-2.15 solar masses) is this same balance with the SAME B as the proton.")

head("A5. The Higgs field: the balance point that sets every elementary mass")
v, mH = 246.22, 125.20
lam = mH**2/(2*v**2); mu = mH/math.sqrt(2); depth = lam*v**4/4
pr(f"  Energy of the Higgs field: -(mu^2) h^2 pushes it AWAY from zero; +lambda h^4 holds it IN. Balance point v = mu/sqrt(lambda) = {mu/math.sqrt(lam):.2f} GeV.")
pr(f"  From the measured Higgs mass: lambda = {lam:.4f}, mu = {mu:.1f} GeV, energy stored at the balance = ({depth**0.25:.0f} GeV)^4.")
pr("  Every elementary particle's mass = (its link strength) x v / sqrt(2). ONE balance point serves all of them; only the link strengths differ.")
for nm, m in (("electron", 0.00051099895), ("muon", 0.1056583755), ("tau", 1.77693), ("charm", 1.273), ("bottom", 4.183), ("top", 172.57)):
    pr(f"    {nm:9s} link strength = {m*math.sqrt(2)/v:.3e}")
pr("  -> 'Mass is a balance point' is TRUE here too. But the balance fixes v, not the link strengths. The open problems are the link strengths.")

head("B1. Does an elementary particle make its own skin? The same balance, applied to the Higgs field")
R = (x0/(4*math.pi*depth))**0.25; Eb = (4/3)*x0/R
pr(f"  A particle could push the Higgs field out of a small region and sit inside it (a bag in the Higgs field). Cost at the balance point: at least {Eb:.0f} GeV")
pr(f"  (radius {R*0.1973269804:.1e} fm), before counting the wall. A particle forms such a skin only if it is HEAVIER than that.")
pr(f"  top: {172.57/Eb*100:.0f}% of the threshold. bottom: {4.183/Eb*100:.1f}%. electron: {0.000511/Eb*100:.1e}%.")
pr("  -> FAILS for every known elementary particle: none is heavy enough to make a skin of its own. They are points riding on the common balance v.")
pr("     Consistent with the old bound that no elementary particle can be much heavier than about 500 GeV (Chanowitz, Furman & Hinchliffe 1978).")

head("B2. The electron as a balance of its own electric push against a containing skin (Poincare 1906)")
re = hc/137.035999/0.51099895
pr(f"  Electric push energy e^2/r equals the electron's mass at r = {re:.2f} fm. Measured: no structure down to about 2e-5 fm (LEP).")
pr(f"  -> A skin that small would store about {re/2e-5:.0e} times the electron's mass in electric push. FAILS by that factor.")
frac = 3/(2*math.pi*137.035999)*math.log(1.22e19/0.51099895e-3)
pr(f"  In quantum theory the electric push adds only {100*frac:.0f}% to the electron's mass even counted down to the Planck length, and it is")
pr("  proportional to the mass already there (Weisskopf 1939). It cannot create the mass.")

head("B3. Is the top quark's link strength a balance point between its own push and the strong force's pull?")
def rge(t, y):
    g1, g2, g3, yt = y; k = 1/(16*math.pi**2)
    return [k*(41/6)*g1**3, -k*(19/6)*g2**3, -k*7*g3**3, k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - (17/12)*g1**2)]
g0 = [0.35830, 0.64779, 1.1666]; mt = 172.57; yt_meas = 0.9369
pr(f"  Exact balance (Pendleton & Ross 1981): link^2 = (2/9) x strong strength^2 -> link = {math.sqrt(2/9)*g0[2]:.3f}; measured {yt_meas:.3f}. Off by {100*(math.sqrt(2/9)*g0[2]/yt_meas-1):+.0f}%.")
for top_scale in (1e16, 1.22e19):
    up = solve_ivp(rge, [0, math.log(top_scale/mt)], g0 + [yt_meas], rtol=1e-9, atol=1e-12)        # gauge strengths at the high scale
    ys = []
    for y_hi in (2.0, 5.0):
        # iterate: run gauge couplings with the large Yukawa (its effect on them is zero at one loop), then run down
        dn = solve_ivp(rge, [math.log(top_scale/mt), 0], list(up.y[:3, -1]) + [y_hi], rtol=1e-9, atol=1e-12); ys.append(dn.y[3, -1])
    pr(f"  Approached balance (Hill 1981), starting from a strong link at {top_scale:.0e} GeV: link at the top's scale = {ys[0]:.3f} to {ys[1]:.3f}; measured {yt_meas:.3f}. Off by {100*(ys[1]/yt_meas-1):+.0f}%.")
pr("  -> FAILS both ways. The measured top sits between the two balance points, on neither.")

head("B4. Is the Higgs itself sitting exactly at the edge of balance?")
mt_, as_m = 172.57, 0.1180
mcrit = 129.6 + 2.0*(mt_ - 173.34) - 0.5*(as_m - 0.1184)/0.0007
err = math.sqrt((2.0*0.5)**2 + (0.5*0.0009/0.0007)**2 + 0.3**2)
logL = 11.3 + 1.0*(mH - 125.66) - 1.2*(mt_ - 173.10) + 0.4*(as_m - 0.1184)/0.0007
pr(f"  The Higgs's containment (lambda) weakens at higher energy. Published fit (Buttazzo et al. 2013), with today's top mass {mt_} and Higgs mass {mH}:")
pr(f"    Higgs mass for containment to last exactly to the Planck scale: {mcrit:.1f} +/- {err:.1f} GeV. Measured {mH}. Gap {mcrit-mH:.1f} GeV = {(mcrit-mH)/err:.1f} sigma.")
pr(f"    Published energy where the Higgs field becomes unstable: about 10^{logL:.1f} GeV (uncertain by about a factor of 10 either way).")
k_ = 1/(16*math.pi**2)
def rge2(t, y):
    """two-loop running of the three force strengths, the top's link and the Higgs's containment lambda (Machacek & Vaughn 1983-85; Luo & Xiao 2003)"""
    g1, g2, g3, yt, l = y; k = k_
    b1 = k*(41/6)*g1**3 + k*k*g1**3*(199/18*g1**2 + 9/2*g2**2 + 44/3*g3**2 - 17/6*yt**2)
    b2 = -k*(19/6)*g2**3 + k*k*g2**3*(3/2*g1**2 + 35/6*g2**2 + 12*g3**2 - 3/2*yt**2)
    b3 = -k*7*g3**3 + k*k*g3**3*(11/6*g1**2 + 4.5*g2**2 - 26*g3**2 - 2*yt**2)
    by = k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - (17/12)*g1**2) + k*k*yt*(-12*yt**4 + 6*l*l - 12*l*yt**2 + yt**2*(131/16*g1**2 + 225/16*g2**2 + 36*g3**2)
         + 1187/216*g1**4 - 0.75*g1**2*g2**2 + 19/9*g1**2*g3**2 - 23/4*g2**4 + 9*g2**2*g3**2 - 108*g3**4)
    bl = k*(24*l*l - 6*yt**4 + 0.375*(2*g2**4 + (g2*g2 + g1*g1)**2) + l*(12*yt*yt - 9*g2*g2 - 3*g1*g1)) + k*k*(-312*l**3 + l*l*(108*g2**2 + 36*g1**2 - 144*yt**2)
         + l*yt**2*(80*g3**2 + 22.5*g2**2 + 85/6*g1**2) - 3*l*yt**4 + 30*yt**6 - 32*g3**2*yt**4 - 8/3*g1**2*yt**4 + l*(-73/8*g2**4 + 39/4*g2**2*g1**2 + 629/24*g1**4)
         + 305/16*g2**6 - 289/48*g2**4*g1**2 - 559/48*g2**2*g1**4 - 379/48*g1**6 + yt**2*(-9/4*g2**4 + 21/2*g2**2*g1**2 - 19/4*g1**4))
    return [b1, b2, b3, by, bl]
def zero_of_lambda(mh, mtop, a3):
    d3 = (a3 - 0.1184)/0.0007                      # starting values at the top's scale: Buttazzo et al. 2013
    y0 = [0.35830 + 0.00011*(mtop-173.34), 0.64779 + 0.00004*(mtop-173.34), 1.1666 + 0.00314*d3 - 0.00046*(mtop-173.34),
          0.93690 + 0.00556*(mtop-173.34) - 0.00042*d3, 0.12604 + 0.00206*(mh-125.15) - 0.00004*(mtop-173.34)]
    T = math.log(1.22e19/mtop); so = solve_ivp(rge2, [0, T], y0, rtol=1e-10, atol=1e-13, dense_output=True)
    ts = np.linspace(0, T, 8000); ls = so.sol(ts)[4]
    return math.log10(mtop*math.exp(ts[int(np.argmax(ls < 0))])), ls[-1]
z0, lpl = zero_of_lambda(mH, mt_, as_m); zlo, _ = zero_of_lambda(mH, mt_ + 0.58, as_m); zhi, _ = zero_of_lambda(mH, mt_ - 0.58, as_m)
pr(f"    My own two-loop check: the containment lambda reaches zero at 10^{z0:.1f} GeV (10^{zlo:.1f} to 10^{zhi:.1f} for the top mass within 2 sigma),")
pr(f"    and at the Planck scale it is {lpl:+.3f} (zero would be exact balance; today's value is {lam:.3f}).")
pr("  -> 'Exactly at the edge' is disfavoured, not excluded. The universe sits within about 2.5% of the edge in the Higgs mass: near balance, not on it.")
pr("     Lead, not a result: the containment fades out near 10^10-10^11 GeV, the same region as TRPG's hidden partners (10^10.2 and 10^11.0 GeV)")
pr("     and close to the right side's scale (10^11.7). That the Higgs's containment runs out near the scale of heavy neutrino partners is noted in")
pr("     published work (Elias-Miro et al. 2012). It is a coincidence within about a factor of ten, not a derivation.")

head("B5. The lepton 45 degree rule IS a balance. Does the idea explain why the balance is 1 to 1?")
m = np.array([0.51099895, 105.6583755, 1776.93]); s = np.sqrt(m); a = s.mean(); selfp = 3*a*a; rel = ((s - a)**2).sum()
m2 = m.copy(); m2[2] += 0.09; s2 = np.sqrt(m2); a2 = s2.mean(); dr = abs(((s2 - a2)**2).sum()/(3*a2*a2) - rel/selfp)
pr(f"  Leptons: (relationship part) / (self part) of the root-mass vector = {rel/selfp:.6f} +/- {dr:.6f}. One to one, to 2 parts in 100,000.")
pr("  For a two-sided balance  E = a x^p + b / x^q  the split at the balance point is q/p, whatever the strengths a and b (as in A1: 1/3).")
pr("  So a 1 to 1 split says the two sides change with size in equal and opposite ways (p = q). That is a restatement, not a derivation:")
pr("  nothing in the balance idea says what 'size' is here or why p = q. Quarks: 1.19 (down family) and 1.55 (up family), with wide errors.")
pr("  -> The 45 degree rule is the balance idea already in TRPG (rule 7). It stays a measured rule. The twist (2/9) is untouched.")

head("C1. Where the forces' SHAPES come from: the same kind of field, three different responses of space")
pr("  A charge sends out a fixed amount of field (flux). What happens to that flux decides the force law.")
pr("   (i)  Space lets the flux spread freely      -> it thins out over a sphere -> force falls as 1/r^2, unlimited reach (electric force, gravity).")
pr("   (ii) Space pushes the flux out (containment) -> it is squeezed into a tube -> force stays CONSTANT with distance (strong force between quarks).")
pr("   (iii) Space is filled with something the flux must push through -> it dies away -> short reach (weak force; the Higgs field does this).")
pr("  Case (ii) is exactly a balance: per unit length, energy = (flux^2)/(2 x tube area) [wants to spread] + B x (tube area) [containment squeezes].")
B, alpha = 60.1, 1.55                                   # TRPG's own values from the proton and the Delta (proton_sets_B)
Phi2 = 4*math.pi*alpha*(4/3)                            # (colour charge of a quark)^2
def per_len(A, B4): return Phi2/(2*A) + B4*A
for Bv in (57.1, 60.1, 63.7):
    B4 = Bv*hc**3; A = minimize_scalar(lambda A_: per_len(A_, B4), bounds=(1e-9, 1e-2), method="bounded", options=dict(xatol=1e-16)).x
    sig = per_len(A, B4); Afm = A*hc**2
    pr(f"   B = {Bv}: balance tube radius {math.sqrt(Afm/math.pi):.2f} fm; field share {Phi2/(2*A)/sig:.2f}, containment share {B4*A/sig:.2f}; "
       f"pull = {sig/1e6:.3f} GeV^2 = {sig/hc/1000:.2f} GeV per fm = {sig/hc*1.602177e-13/1e-15/9.80665/1000:.1f} tonnes")
sig0 = math.sqrt(32*math.pi*alpha*B*hc**3/3)
pr(f"  Closed form sqrt(32 pi alpha B / 3) = {sig0/1e6:.3f} GeV^2 (Johnson & Thorn 1976). With the 1975 bag numbers (B = 57.5, alpha = 2.2) it gives {math.sqrt(32*math.pi*2.2*57.5*hc**3/3)/1e6:.3f}.")
pr("  Measured pull between quarks: 0.18-0.19 GeV^2 (about 0.9-1.0 GeV per fm, 15-16 tonnes), from heavy-quark pairs, lattice QCD and the spin ladder of mesons.")
pr(f"  -> The balance gives a constant pull of the right size from the proton's own two numbers, with nothing new put in: {100*(sig0/1e6/0.185-1):+.0f}% against 0.185.")
pr("     Right kind of force, right size; as a strict prediction it is low and FAILS. (The strength alpha here is the bag model's own, and is known to be rough.)")
pr("     The tube also comes out too wide (1.4 fm in radius; lattice QCD sees roughly 0.5 fm). Note the split: exactly half field, half containment.")

head("C2. The weak force's short reach comes from the Higgs balance point")
g2 = 0.6517
MW = g2*v/2; rng_w = hc/1000/MW
pr(f"  The W's mass = (weak strength) x (Higgs balance point) / 2 = {g2} x {v} / 2 = {MW:.1f} GeV. Measured 80.37 GeV. Reach = {rng_w:.1e} fm.")
pr("  This was predicted before the W was found (Weinberg 1967; Salam 1968; found at CERN 1983). Honest note: the balance point v is itself")
pr("  read off the weak force at low energy, so this is a consistency check that passes, not an independent derivation.")

head("C3. The force between protons and neutrons IS shared leakage")
mpi = 138.04; gA_bag = (5/3)*(1 - (2*x0 - 3)/(3*(x0 - 1)))
pr(f"  The proton's pion cloud (the 'leakage' of the record's proton step) reaches a neighbour and binds it (Yukawa 1935). Reach = hbar/(m_pi c) = {hc/mpi:.2f} fm.")
for lab_, gA in (("bag balance alone", gA_bag), ("measured g_A", 1.2754)):
    g = gA*938.92/92.3; f2 = (g*mpi/(2*938.92))**2/(4*math.pi)
    pr(f"   strength from the {lab_:18s}: g_A = {gA:.3f} -> f^2/4pi = {f2:.4f} (measured 0.0778: {100*(f2/0.0778-1):+.0f}%)")
pr("  -> Right kind of force and right reach. The simple balance gives the strength about a third low; the record uses the measured value.")

head("C4. Do the forces' basic STRENGTHS come from a balance? (strength made entirely by the matter around it: Sakharov 1967 for gravity; Zel'dovich 1967)")
pr("  If a force's strength were made entirely by the particles it acts on, its inverse strength would run to zero at some high energy.")
ia = dict(hypercharge=(59.0, -41/10), weak=(29.59, 19/6), strong=(1/0.1180, 7.0))
for nm, (i0, b) in ia.items():
    if b < 0: pr(f"   {nm:11s}: inverse strength {i0:.1f} today, falling: reaches zero at 10^{math.log10(91.19) + (i0*2*math.pi/(-b))/math.log(10):.0f} GeV (far beyond the Planck energy, 10^19)")
    else: pr(f"   {nm:11s}: inverse strength {i0:.1f} today, RISING with energy: never reaches zero")
pr("  -> FAILS for the weak and strong forces with the known particles. Their strengths are not made by this kind of balance.")
pr("     One true thing nearby: the strong force's strength is not a free number. It is fixed by ONE scale, the scale of the proton's balance.")

head("C5. What the balance cannot say")
pr("  A balance needs two things to balance. It explains the shape, reach and (roughly) the size of each force once the force exists.")
pr("  It does not say WHICH forces exist (8 | 3 | 1), why the weak force is left-handed, or what the three basic strengths are.")
Nstar = (hc*1e6*1.602177e-19*1e-15/(G*mp*mp))**1.5
pr(f"  It does set the sizes of things from the strengths: e.g. a star's mass scale, where gravity balances quantum push, is (hbar c / G m_p^2)^(3/2) protons = {Nstar*mp/Ms:.2f} solar masses.")

head("B6. Summary: does the balance idea make the open problems easier?")
rows = [("bottom, charm, strange masses", "no: they are link strengths; the balance point v is common to all (A5)"),
        ("the top's link", "no: both balance points fail (B3)"),
        ("the 2/9 twist", "no: untouched (B5)"),
        ("lepton 45 degrees", "already this idea (rule 7); not derived by it (B5)"),
        ("why the forces look so different", "YES: shape and reach of each force follow from how space treats its flux (C1-C3)"),
        ("pull between quarks", "roughly: 0.155 against 0.18-0.19 GeV^2 from the proton's own numbers (C1); strict FAIL"),
        ("which forces, handedness, strengths", "no (C4, C5)"),
        ("matter excess", "no numbers; only the picture (the burst outward)"),
        ("where does mass come from at all", "YES: 99% is the proton's balance (A1), the rest the Higgs balance (A5). No magic."),
        ("one law at every scale", "YES: proton, nucleus, star, vacuum (A1-A5), with TRPG's one B for proton and neutron star")]
for a_, b_ in rows: pr(f"  {a_:36s} {b_}")
