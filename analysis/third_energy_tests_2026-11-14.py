"""Idea tested: the Higgs is not the whole thing; place the Higgs with gravity and a different energy in the missing scale.
Tested here, each one attacked:
  1. Is there an energy in the matrix that the Higgs does not give, and what allows it?
  2. Does a force go with it, and which of the two kinds is it?
  3. Can that energy tell the three families apart?
  4. What the LHC says about a family-blind Higgs plus a second source of mass.
  5. The Higgs value, dark energy and gravity as 'constants of space'.
  6. The unplaced pieces: the coloured far-side parts, and dark matter.
Uses forces_build_the_atom_2026-11-14.py, perspective_and_gravity_2026-11-14.py and three_way_relationship_2026-11-14.py (same folder).
Credit: Majorana 1937; Minkowski 1977; Yanagida 1979; Gell-Mann, Ramond & Slansky 1979; Mohapatra & Senjanovic 1980; Davidson 1979;
Marshak & Mohapatra 1980; Wilczek & Zee 1979; Georgi 1975; Fritzsch & Minkowski 1975; Chikashige, Mohapatra & Peccei 1981; Krauss & Wilczek 1989;
Martin 1992; Kadastik, Kannike & Raidal 2010; Furry 1939; Carena, Daleo, Dobrescu & Tait 2004; Wagner, Schlamminger, Gundlach & Adelberger 2012; Koide 1982; Brannen 2006;
Rivero & Gsponer 2005; ATLAS and CMS 2022-2023; Lange et al. 2021; Hofmann & Mueller 2018; DESI 2025; Planck 2018; D'Onofrio & Rummukainen 2016;
Arkani-Hamed, Hall, Kolda & Murayama 2000; Fardon, Nelson & Weiner 2004; Elias-Miro et al. 2012; Buttazzo et al. 2013; Nath & Fileviez Perez 2007;
KamLAND-Zen 2024; Cowsik & McClelland 1972; PDG 2024."""
import os, io, contextlib, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src5 = open(os.path.join(HERE, "forces_build_the_atom_2026-11-14.py")).read()
with contextlib.redirect_stdout(io.StringIO()):                       # the 12 turns, the far side, 'blind to the facing', the settled Higgs side
    exec(src5[src5.index("import os, io"):src5.index('head("3. The atom from the charges")')])

def kindval(gen):                                                     # value of a single turn on each kind of state
    out = {}
    for nme in names_:
        Pk = U[:, kinds[nme]]; out[nme] = sorted(set(np.round(np.linalg.eigvalsh(herm(Pk.conj().T@(-1j*gen)@Pk)), 4).tolist()))
    return out
def sym_count(bs):                                                    # how many of the allowed pairings are symmetric (a mass needs a symmetric one)
    if not bs: return 0
    S = np.array([((b + b.T)/2).ravel() for b in bs]); return int((np.linalg.svd(S, compute_uv=False) > 1e-8).sum())

head("1. An energy the Higgs does not give: which pairs of states can link with NOTHING else present?")
G13r = restrict(Kd, two(ray(e1), ray(e2))); G11r = restrict(Kd, two(exact(e1), exact(e2)))
G13 = [C(g) for g in G13r]; G11 = [C(g) for g in G11r]
b13 = invariant_bilinear(G13); b11 = invariant_bilinear(G11)
pr("  A 'Higgs-free' mass is a pairing of two states that every kept turn leaves alone. Count of such pairings (symmetric ones, as a mass needs):")
pr(f"     all 45 turns (no perspective chosen)                         : {sym_count(b45)}")
pr(f"     13 turns (spot and facing both kept only as directions)      : {sym_count(b13)}")
pr(f"     12 turns (rule 4: the spot kept EXACTLY, the facing a direction): {sym_count(b12)}")
pr(f"     11 turns (both kept exactly)                                 : {sym_count(b11)}")
Bm = b12[0]; Bs = (Bm + Bm.T)/2; wgt = np.array([[abs(U[:, a_]@Bs@U[:, b_])**2 for b_ in range(16)] for a_ in range(16)]); a0, b0 = np.unravel_index(wgt.argmax(), wgt.shape)
hid = kinds["hidden neutrino partner"][0]; ov = abs(np.vdot(U[:, hid], e1))**2/np.vdot(e1, e1).real
pr(f"  Under the 12, the one pairing is the hidden neutrino partner with itself: {a0 == hid and b0 == hid} (weight {wgt[a0, b0]/wgt.sum():.3f}).")
pr(f"  And that state IS the state rule 4 keeps exactly, the observer's 'spot' state of an earlier test (overlap of the two: {ov:.6f}).")
pr("     (This had to be so: a state kept exactly is untouched by all 12 turns, so it carries no charge of any kind, and the only such state is the")
pr("     hidden neutrino partner. Georgi 1975; Fritzsch & Minkowski 1975.)")
# the 13th turn: the part of the 13 that is not in the 12
V12 = np.array([g.ravel() for g in G]); q12, _ = np.linalg.qr(V12.T); V13 = np.array([g.ravel() for g in G13r])
rest = V13 - (V13@q12)@q12.T; u_, s_, v_ = np.linalg.svd(rest, full_matrices=False); n13 = int((s_ > 1e-8*s_[0]).sum())
T13 = C(v_[0].reshape(32, 32)); Ph = U[:, kinds["hidden neutrino partner"]]; sc = 5.0/(Ph.conj().T@(-1j*T13)@Ph).real[0, 0]; Xg = T13*sc; xv = kindval(Xg)
comm = max(np.abs(Xg@g - g@Xg).max() for g in gens12)
pr(f"  The 13 hold the 12 and exactly {n13} more turn. It commutes with all 12 (error {comm:.0e}), so it is a single one, like the electric-type one.")
pr(f"     Its values (scaled so the spot state has 5): { {k_: v_ for k_, v_ in xv.items()} }")
BLg = (Xg + 4*Yg)/5; blv = kindval(BLg)
pr(f"     Mixed with the earlier single one, (this + 4 x that)/5: { {k_: v_ for k_, v_ in blv.items()} }")
pr("     Read as particles (partners are written here as their antiparticles): quarks +1/3, leptons -1. That is 'matter number':")
pr("     (number of quark-triples) minus (number of leptons), B-L in the physics literature (Davidson 1979; Marshak & Mohapatra 1980).")
pr(f"  The spot state carries {blv['hidden neutrino partner'][0]:+.0f} of it, so the spot's pairing with itself carries {2*blv['hidden neutrino partner'][0]:+.0f}.")
pr("  -> The hidden partner's own mass is allowed EXACTLY BECAUSE rule 4 keeps the spot exactly: keeping the spot exactly is the same act as giving up")
pr("     the 13th turn. With no perspective (45) or with the spot kept only as a direction (13), no mass of any kind exists without the Higgs.")
pr("     So the matrix holds exactly ONE energy the Higgs does not give, and it belongs to the observer's own spot. the 'different energy'")
pr("     is there. It is the standard heavy-partner mass (Majorana 1937; Minkowski 1977; Yanagida 1979; Gell-Mann, Ramond & Slansky 1979).")
pr("  ATTACK 1: with both states kept exactly (11 turns) more pairings appear: the electron partner would get a Higgs-free mass of its own, and the")
pr(f"     single carrier is lost altogether. That pairing's electric charge is {2*kindval(Qg)['electron partner'][0]:+.1f}, not zero, so the measured world rules it out.")
pr("     'Spot exact, facing a direction' is needed. It is the same reading an earlier test already needed, so no new choice is added.")
pr("  ATTACK 2: the pairing is ALLOWED, not REQUIRED, and its SIZE is not fixed by the 12 (it is unchanged by all of them). Its size is an input.")

head("2. The turn that is given up: is it a force, and which of the two kinds?")
gens13 = gens12 + [BLg]
tests = [("strong + matter number (9)", col8 + [BLg]), ("strong + electric (9)  [an earlier test]", col8 + [Qg]), ("strong + electric + matter number (10)", col8 + [Qg, BLg]),
         ("strong + the earlier single one + matter number (10)", col8 + [Yg, BLg]), ("strong + weak + matter number (12)", col8 + wk3 + [BLg]), ("all 13", gens13)]
for nme, gs in tests: pr(f"     {nme:<54} blind to the facing: {facing_blind(gs)}")
grid = np.round(np.arange(-2, 2.0001, 0.25), 3); ok = sorted(set(float(a) for a in grid for c in (0.0, 0.7, -1.3) if facing_blind(col8 + [Yg + a*t3 + c*BLg])))
pr(f"  Scan: strong + ONE mixture (earlier single one + a x weak turn + c x matter number), a from -2 to 2, c = 0, 0.7, -1.3: blind only for a = {ok}, any c.")
pr("  -> Inside the 13, the LARGEST part that is blind to the facing has 10 turns: strong + electric + matter number. (Any 10 or more must hold")
pr("     the whole strong force plus two single mixtures; the scan shows the mixtures must be electric (a = +1 or -1) and matter number.)")
pr("     So the turn given up at the spot is of the FIRST kind, the push/pull kind, the same kind as the strong and electric forces.")
rho13 = [act10(g)[0] for g in gens13]
def kept13(h):
    A = np.array([np.concatenate([(r_@h).real, (r_@h).imag]) for r_ in rho13]).T; u2, s2, v2 = np.linalg.svd(A); return v2[(s2 > 1e-9*s2[0]).sum():]
N13 = kept13(h0); k13 = [sum(c*g for c, g in zip(cv, gens13)) for cv in N13]; blh = np.linalg.norm(rho13[-1]@h0)/np.linalg.norm(h0)
blc = np.linalg.norm(rho13[-1]@hc)/np.linalg.norm(hc)
pr(f"  The settled Higgs side (45.3) keeps {len(N13)} of the 13 (blind to the facing: {facing_blind(k13)}). Matter number of the Higgs pair: {blh:.1e}; of the coloured far-side parts: {blc:.3f}.")
pr("  -> The Higgs does NOT give up matter number. Only the spot does. Two separate settlings, each giving up one thing and giving one kind of mass:")
pr("        the far side settles  : gives up the weak turns   -> every mass that links a pair to its partner (the Higgs link)")
pr("        the spot is kept exactly: gives up matter number  -> the one mass that needs no Higgs (the hidden partner's own)")
pr("     Today both have happened: 13 -> 9 turns = strong + electric. That is why matter number is not a force we feel.")
pr("  In ordinary matter the source of this force would be: proton +1, electron -1, neutron +1. Hydrogen: 0. Any neutral atom: its number of neutrons.")
gL = 0.56; MN = [7.2e5, 1.3e9, 6.3e10]; hb = 0.19733
pr("  If it were long-reach it would pull on neutrons only and break 'everything falls alike'. Torsion balances limit such a pull to a strength")
pr("     below about 1e-24 (as a charge; Wagner, Schlamminger, Gundlach & Adelberger 2012); the forces here have about 0.5. So its carrier must be heavy.")
pr("  How heavy: the same settling gives the carrier's mass and the hidden partner's mass, so the carrier weighs about as much as the heaviest hidden")
pr(f"     partner or more (unless that partner's link to the settling is several times the force's own strength): about {MN[2]:.1e} GeV with rule 5's link. Reach {hb/MN[2]:.1e} fm,")
pr(f"     which is {1.58e-3/(hb/MN[2]):.0e} times shorter than the Higgs link's reach. Its pull at ordinary energies, against the weak force's: {(246.22/(2*MN[2]/gL))**2:.0e}.")
pr(f"     Collider limit: carrier mass / strength > about 7,000 GeV (LEP; Carena, Daleo, Dobrescu & Tait 2004). Here: {MN[2]/gL:.0e}. Passed, by a wide margin.")
pr("  ATTACK (the reading): there are two ways to hold rule 4, and no experiment today separates them.")
pr("     (i)  The spot is exact always. Then matter number is never a force: it is a turn the perspective excludes, with no carrier at any energy,")
pr("          and the hidden partner's mass is simply an allowed number.")
pr("     (ii) The spot becomes exact below some energy (as the Higgs side settles below 160 GeV). Then above that energy matter number is a real")
pr("          third push/pull force, and its carrier has the mass worked out above.")
pr("     Both give the same world at every energy reached so far. Reading (ii) is the one in which the 'different energy' is a scale.")
pr("  -> A third push/pull force IS in the matrix, with its own charge (matter number). It binds nothing and cannot be detected directly.")
pr("     ATTACK: 'cannot be detected' means this is not a confirmation. Its footprints are only: neutrino masses; neutrinos being their own")
pr("     antiparticles; and the matter excess. Not new physics: it is the B-L force of Georgi 1975 and Fritzsch & Minkowski 1975,")
pr("     and the 4x4 alternative already had it. What is new here is where it sits in the main line: it is the turn that")
pr("     keeping the spot exactly removes.")

head("3. Can that energy tell the three families apart?")
me, mmu, mtau = 0.51099895e-3, 0.1056583755, 1.77693; mnu = [0.36e-12, 8.6e-12, 50.0e-12]; vh = 174.104
pr(f"  (a) Charged particles. Higgs-free pairings that involve any charged state: {sym_count(b12) - 1}. So this energy cannot enter the electron's, muon's, tau's")
pr("      or any quark's mass directly. Indirectly (through a loop) it changes a charged lepton's mass by a fraction")
for nme, m_ in (("electron", me), ("muon", mmu), ("tau", mtau)):
    y = m_/vh; pr(f"         {nme:<9} {y*y/(16*math.pi**2)*math.log(MN[2]/91.19):.1e}   (link {y:.1e})")
pr("      of a mass that is already there. It multiplies; it cannot separate. FAIL for the charged particles.")
Ms = [m_**2/n_ for m_, n_ in zip((me, mmu, mtau), mnu)]
pr("  (b) Neutrinos. The light neutrino's mass = (its Higgs link)^2 / (the hidden partner's own mass). With rule 5's link and TRPG's neutrino")
pr(f"      masses (0.36, 8.6, 50 meV) the three values of the third energy are {Ms[0]:.1e}, {Ms[1]:.1e}, {Ms[2]:.1e} GeV.")
ang = lambda m: math.degrees(math.acos(sum(math.sqrt(x) for x in m)/math.sqrt(3*sum(m))))
pr(f"      TEST 'one value for all three families': the neutrino masses would then go as the lepton masses squared: m2/m3 = {(mmu/mtau)**2:.4f}; measured {8.6/50:.3f}.")
pr(f"      FAIL by a factor {(8.6/50)/(mmu/mtau)**2:.0f}. So this energy must itself differ between the families.")
pr(f"      TEST 'its three values follow the 45 degree balance': charged leptons {ang((me, mmu, mtau)):.2f} deg; the three hidden-partner values {ang(Ms):.1f} deg; the neutrinos {ang(mnu):.1f} deg. FAIL.")
pr("      Three values, no pattern, fixed only by measured lepton and neutrino masses. It CARRIES the family difference for the neutrinos;")
pr("      it does not EXPLAIN it.")
pr("  (c) States. One matrix has one spot, so one hidden partner and one such energy. Three need three placements, which is the reading TRPG")
pr("      already uses. The energy is a pairing of a state that is already counted; it adds no states. The count 24 against 48 (45.8) is unchanged.")
pr("  -> The 'different energy' exists, but it does not supply the missing families. FAIL as the source of the three generations.")

head("4. What the LHC says: is the Higgs blind to the family, with something else supplying the difference?")
kap = {"top": (1.01, 0.11), "bottom": (0.99, 0.16), "tau": (0.92, 0.08), "muon": (1.12, 0.21)}
pr("  Measured Higgs link / (mass / 246 GeV), which is 1 if the Higgs gives the whole mass (CMS 2022; ATLAS 2022 agrees):")
for k_, (c_, e_) in kap.items(): pr(f"     {k_:<7} {c_:.2f} +/- {e_:.2f}    largest share of the mass from another source (2 sigma): {max(0.0, 1 - (c_ - 2*e_))*100:.0f}%")
r = mtau/mmu; kt, et = kap["tau"]; km, em = kap["muon"]
pr(f"  TEST 'family-blind Higgs' (same link to the muon and the tau): the muon's number above would be {r*kt:.1f}; measured {km:.2f} +/- {em:.2f}: FAIL by {(r*kt - km)/em:.0f} sigma.")
pr(f"     In rates: Higgs -> two taus would equal Higgs -> two muons. Measured: {r*r:.0f} x ({(kt/km)**2:.2f} +/- {2*(kt/km)**2*math.hypot(et/kt, em/km):.2f}) to one, in line with the standard {r*r:.0f}.")
pr("  So the Higgs itself tells the tau from the muon, and it links to the top and the bottom in proportion to their masses. For the second and")
pr("  third families a second source of mass is limited to the shares above. FAIL for 'the Higgs is a family-blind constant'.")
ke = math.sqrt(3.0e-4/5.0e-9)
pr(f"  HONEST GAP: the first family is not measured. Electron: link < {ke:.0f} x standard (CMS 2023). Charm: < 8.5 x (ATLAS 2022). Up, down, strange: no limit that matters.")
pr("     Experiment alone would still allow a different source for the lightest masses. The matrix does not (test 3a): it predicts 1.00 for every one.")
pr("  WHAT SURVIVES: 'the Higgs sets the scale, the geometry sets the pattern' (mass = one Higgs number x a pattern from the matrix). That is exactly")
pr("     the form rules 5 and 7 already have, and the LHC cannot tell it from the standard picture, because the link then still equals mass / 246 GeV.")

head("5. The Higgs value, dark energy and gravity as 'constants of space'")
v = 246.22; Mpl = 1.22089e19; Mred = Mpl/math.sqrt(8*math.pi); de = 2.25e-12; meet = 10**16.46
pr("  (a) Are they constant? Measured limits on change:")
pr("        electron/proton mass ratio (Higgs value against the strong balance): (-8 +/- 36)e-18 per year (Lange et al. 2021)")
pr("        gravity's strength: (7.1 +/- 7.6)e-14 per year (Hofmann & Mueller 2018)")
pr("        dark energy: pressure/density = -1.03 +/- 0.03 (Planck 2018 with others); DESI 2025 hints at change at 2.8 to 4.2 sigma (live test, 44.5)")
pr("      HOLDS as far as measured: all three are properties of space that do not depend on what matter is present.")
pr("  (b) Are they the same KIND of constant? The Higgs value melts: above a temperature of 159.5 +/- 1.5 GeV it is zero (D'Onofrio & Rummukainen 2016).")
pr("      Gravity's strength and the kept amount of the bag do not melt. The Higgs particle has spin 0 (measured; spin 2 excluded above 99.9%, ATLAS and")
pr("      CMS 2013-2015); gravity's waves have spin 2. And the Higgs does not count all energy (45.7, a factor 130). So the Higgs value is a STATE")
pr("      that space has settled into; gravity's strength and the bag's amount are fixed properties. Same shelf, different kind. FAIL as 'the same thing'.")
pr("  (c) Are they one number? Three scales: Higgs value 246 GeV; gravity's grain 1.22e19 GeV; dark energy 2.25 meV.")
pr(f"        Higgs^2 / grain            = {v*v/Mpl*1e12:.2e} meV  (dark energy scale is {de/(v*v/Mpl):.0f} times this)   FAIL  (Arkani-Hamed, Hall, Kolda & Murayama 2000)")
pr(f"        Higgs^2 / reduced grain    = {v*v/Mred*1e12:.2e} meV  (dark energy scale is {de/(v*v/Mred):.0f} times this)    FAIL")
pr(f"        Higgs^2 / meeting energy   = {v*v/meet*1e12:.2f} meV      (dark energy scale 2.25 meV: {100*(v*v/meet/de - 1):+.0f}%)")
pr(f"        the same with 174 GeV      = {vh*vh/meet*1e12:.2f} meV      ({100*(vh*vh/meet/de - 1):+.0f}%)")
pr("      The third line has the same form as the neutrino's mass in test 3 (Higgs^2 over a heavy energy), with the energy where the strong and weak")
pr("      forces meet (10^16.46 GeV, an earlier test) in place of the hidden partner's mass. ATTACK: I tried 3 heavy energies x 2 ways of writing the Higgs")
pr("      value, and one of the six landed; the meeting energy itself is known only to within a factor of about 2, so '7%' is not meaningful;")
pr("      a hit that close among 6 tries spread over a few powers of ten happens by chance roughly 1 time in 6; and there is no mechanism.")
pr("      COINCIDENCE, recorded so it is not rediscovered. Not used. Not a result.")
pr(f"      (Already recorded: the dark-energy scale lies inside TRPG's neutrino mass range, 0.36 to 50 meV. Same status.)")
pr("  (d) Does the third energy explain where the Higgs's hold fades (10^9.8 to 10^11.5 GeV, an earlier test)? The heaviest hidden partner sits in that window")
pr(f"      (10^{math.log10(Ms[2]):.1f}). ATTACK: the window is 1.7 of about 17 powers of ten, and there are three partners: chance {1 - (1 - 1.7/17)**3:.0%}. And the partner's")
pr(f"      link ({mtau/vh:.3f}) is far too weak to cause the fade (it would need about 0.4 or more; Elias-Miro et al. 2012). FAIL as a mechanism; weak coincidence.")
pr("  -> Grouping the Higgs value with gravity and dark energy HOLDS as a classification: they are properties of space, not turns of the matrix")
pr("     (turns give the forces; the far side settling gives the Higgs value; space itself gives gravity and the kept amount).")
pr("     It FAILS as 'one thing' and as 'one number'. It removes no measured number.")

head("6. The unplaced pieces: the six coloured far-side parts, and dark matter")
par16 = sorted(set(np.round(np.cos(3*math.pi*np.linalg.eigvalsh(herm(-1j*BLg))), 6).tolist()))
bl10 = np.linalg.eigvals(-1j*rho13[-1]).real; par10 = sorted(set(np.round(np.cos(3*math.pi*bl10), 6).tolist()))
pr(f"  Matter number of the far side's 10 parts: {sorted(set((np.round(bl10, 3) + 0.0).tolist()))}. The coloured ones carry colour, electric charge 1/3 and matter number 2/3:")
pr("     not dark. As particles they would let protons decay, so they must weigh more than about 1e11 GeV (Nath & Fileviez Perez 2007) or be read")
pr("     as entries, not particles (44.2). They are not the third energy and not dark matter. Unchanged.")
pr(f"  A sign that survives BOTH settlings: (-1)^(3 x matter number). On all 16 states: {par16}. On all 10 far-side parts: {par10}. The spot's pairing carries 2,")
pr("     so the sign is kept exactly (Krauss & Wilczek 1989; Martin 1992). Consequence: the lightest thing with sign -1 can never decay.")
sm = sum(mnu)*1e9; om = sm/93.12
pr(f"  In the main line every -1 thing is ordinary matter, and the lightest is the lightest neutrino. Neutrinos left from the Big Bang: {sm*1e3:.0f} meV in total")
pr(f"     gives {om:.1e} of the critical density (x h^2); dark matter is 0.120. That is {100*om/0.120:.1f}% of it, and it moves too fast to form galaxies")
pr("     (Cowsik & McClelland 1972; Planck 2018). FAIL as dark matter.")
GamN = (me/vh)**2*Ms[0]/(8*math.pi)
pr(f"  The lightest hidden partner: decays in {6.582e-25/GamN:.0e} s. FAIL.")
pr("  -> The matrix supplies the PROTECTION a dark particle would need (an exact sign) but no particle to protect: that would take a state with")
pr("     sign +1 and spin 1/2, or sign -1 and spin 0 (Kadastik, Kannike & Raidal 2010), and rule 6 ('nothing extra') has none. Dark matter still has NO candidate. The open item is")
pr("     now exact: either rule 6 gives way by one such state, or dark matter is not a particle of the matrix at all.")
Ue = (0.68, 0.30, 0.022); tms = [u*m*1e12 for u, m in zip(Ue, mnu)]
pr("  One prediction falls out of test 1: the third energy exists only if neutrinos are their own antiparticles (Majorana 1937). Then a rare decay")
pr(f"     (double beta decay with no neutrinos; Furry 1939) happens with an effective mass between {abs(tms[1] - tms[2] - tms[0]):.1f} and {sum(tms):.1f} meV. Today's limit: 28 to 122 meV")
pr("     (KamLAND-Zen 2024). Planned experiments reach about 10 to 20 meV. So it is a real prediction that cannot be checked for a long time.")

head("7. Summary")
pr("  HOLDS  : there is exactly one energy in the matrix that the Higgs does not give. It belongs to the observer's own spot and is allowed because")
pr("           rule 4 keeps the spot exactly. A third push/pull force (charge: matter number) goes with it and is given up in the same act.")
pr("  HOLDS  : (classification) the Higgs value, gravity's strength and dark energy are properties of space, not turns of the matrix.")
pr("  FAILS  : that energy as the source of the three families (cannot touch charged masses; three unpatterned values for the neutrinos; adds no states).")
pr("  FAILS  : a family-blind Higgs (LHC, about 70 sigma for muon against tau); the Higgs as the same kind of constant as gravity; any one-number relation.")
pr("  FAILS  : dark matter from any of the unplaced pieces.")
pr("  NOT FIXED: the size of the third energy (input); why three families; everything on the open list of 44.7 / 45.9.")
