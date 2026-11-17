"""(A) drop 'left-handed / right-handed' as a concept and keep only the observer's perspective: does everything
the main line has survive? (B) Then start with gravity's strength: gravity is tied to entanglement, so it may tie to the Higgs field, and to TRPG's dark
energy and dark matter; and quark mixing may come out. Each step is tested here and attacked. Uses three_way_relationship_2026-11-14.py (same folder).
Credit: Weyl 1929; Lee & Yang 1956; Wu 1957; Landau 1957; Salam 1957; Lueders 1954; Pauli 1955; Stueckelberg 1960; Renou et al. 2021; Georgi 1975;
Fritzsch & Minkowski 1975; Freudenthal 1954; Minkowski 1977; Gell-Mann, Ramond & Slansky 1979; Yanagida 1979; Bekenstein 1973; Hawking 1975;
Bombelli, Koul, Lee & Sorkin 1986; Srednicki 1993; Jacobson 1995; Sakharov 1967; Dvali 2007; Gross & Wilczek 1973; Politzer 1973; Wilczek 2001;
Georgi, Quinn & Weinberg 1974; Kaplunovsky 1988; Zee 1979; Atkins & Calmet 2013; Shaposhnikov & Wetterich 2010; Buttazzo et al. 2013; Einstein 1919;
Weinberg 1989; Padmanabhan; Verlinde 2017; McGaugh, Lelli & Schombert 2016; Planck 2018; DESI 2025; Dodelson & Widrow 1994; Asaka & Shaposhnikov 2005;
Gatto, Sartori & Tonin 1968; Fritzsch 1977; Kobayashi & Maskawa 1973; PDG 2024; FLAG 2024."""
import os, io, contextlib, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src3 = open(os.path.join(HERE, "three_way_relationship_2026-11-14.py")).read()
with contextlib.redirect_stdout(io.StringIO()):                       # build the three-way product and its 12 turns exactly as an earlier test did
    exec(src3[src3.index("import math, itertools"):src3.index('head("6. Attacks")')])

def full_commutant(mats):                                            # every real matrix that commutes with all of mats
    n_ = mats[0].shape[0]; I = np.eye(n_); Gm = np.zeros((n_*n_, n_*n_))
    for A in mats:
        M = np.kron(A, I) - np.kron(I, A.T); Gm += M.T@M
    w, Vv = np.linalg.eigh(Gm); k = int((w < 1e-9*w[-1]).sum()); return [Vv[:, i].reshape(n_, n_) for i in range(k)]
def invariant_bilinear(cmats):                                       # complex B with g^T B + B g = 0: a pairing of two states that needs nothing else
    Gm = np.zeros((256, 256), complex); I = np.eye(16)
    for g in cmats:
        M = np.kron(g.T, I) + np.kron(I, g.T); Gm += M.conj().T@M
    w, Vv = np.linalg.eigh(Gm); k = int((w < 1e-9*w[-1]).sum()); return [Vv[:, i].reshape(16, 16) for i in range(k)]

head("A1. How many 'i' does the three-way matrix have? (An 'i' = a way to turn every state a quarter turn that all the turns respect)")
c46 = full_commutant(K); c45 = full_commutant(Kd); c12 = full_commutant(G)
J = [m for m in c46 if abs(np.trace(m)) < 1e-6*np.linalg.norm(m)] or [c46[0] - np.trace(c46[0])/32*np.eye(32)]
Jm = c46[0] - np.trace(c46[0])/32*np.eye(32)
if np.linalg.norm(Jm) < 1e-8: Jm = c46[1] - np.trace(c46[1])/32*np.eye(32)
Jm = Jm/math.sqrt(-np.trace(Jm@Jm)/32)
pr(f"  Matrices that commute with all 46 turns of the three-way product: {len(c46)} (the number 1 and one 'i'). With the 45: {len(c45)}.")
pr(f"  That one 'i' squared = -1 (error {np.abs(Jm@Jm + np.eye(32)).max():.0e}). So the whole three-way structure has exactly ONE i, up to its sign.")
pr(f"  Matrices that commute with only the 12 force carriers: {len(c12)}. With the 12 alone, each of the six kinds of state could be given its own")
pr("     i, with its own sign: each kind could be made left or right separately. The whole structure does not allow that: one sign for all.")
pr("  -> There are two choices and only two: i or -i. That is the observer's facing direction, forward or backward.")

head("A2. Is the world its own mirror image? And which pairs of states can link with nothing else present?")
ysr = np.round(ys, 4) + 0.0; pos_ = sorted(ysr.tolist()); neg_ = sorted((-ysr + 0.0).tolist())
pr(f"  Charges of the 16 states: {sorted(set(pos_))}. Their mirror would need: {sorted(set(neg_))}. Same list: {pos_ == neg_}.")
G12 = [C(g) for g in G]; K45 = [C(g) for g in Kd]
b12 = invariant_bilinear(G12); b45 = invariant_bilinear(K45)
pr(f"  Pairings of two states allowed with nothing else present: under all 45 turns: {len(b45)}; under the 12 carriers: {len(b12)}.")
if b12:
    Bm = b12[0]; wgt = np.array([[abs(U[:, a_]@Bm@U[:, b_]) for b_ in range(16)] for a_ in range(16)]); a0, b0 = np.unravel_index(wgt.argmax(), wgt.shape)
    pr(f"     That one pairing joins a state of charge {rows[a0][2]+0.0:+.3f} (colours {'3' if rows[a0][0] > 1 else '1'}, {'pair' if rows[a0][1] > .5 else 'single'}) with itself: weight on it {wgt[a0, b0]**2/ (wgt**2).sum():.3f}.")
pr("  -> The list of states is not its own mirror image, and no state can pair off with another by itself, except the hidden neutrino partner")
pr("     with itself. That one exception is the standard reason the hidden partners are very heavy (Minkowski 1977; Yanagida 1979).")

head("A3. Mass without the words 'left' and 'right': the link through the far side of the matrix")
one = np.zeros(dim, complex); one[:3] = 1; tr = lambda z: z[0] + z[1] + z[2]
chk = max(np.abs(jp(one, emb(E16[i])) - emb(E16[i])).max() for i in range(16)); e3 = np.zeros(dim, complex); e3[2] = 1
half = max(np.abs(jp(e3, emb(E16[i])) - 0.5*emb(E16[i])).max() for i in range(16))
far = [0, 1] + [i for i in range(3, dim) if i not in lk]
pr(f"  Checks: (1,1,1) is the unit (error {chk:.0e}); the link states sit half-way from the spot (error {half:.0e}); the far side has {len(far)} parts.")
def cross(x, y):
    p = jp(x, y); return p - 0.5*tr(p)*one                           # Freudenthal's product for two traceless elements
Mm = np.zeros((len(far), 16, 16), complex); spotpart = 0.0
for a_ in range(16):
    for b_ in range(16):
        r_ = cross(emb(E16[a_]), emb(E16[b_])); Mm[:, a_, b_] = r_[far]; spotpart = max(spotpart, abs(r_[2]), np.abs(r_[lk]).max())
Mf = Mm.reshape(len(far), 256); rk = np.linalg.matrix_rank(Mf, tol=1e-9); Mp = np.linalg.pinv(Mf)
def act10(g):
    rhs = (np.einsum('kab,ac->kcb', Mm, g) + np.einsum('kab,bd->kad', Mm, g)).reshape(len(far), 256); rho = rhs@Mp
    return rho, np.abs(rho@Mf - rhs).max()
res10 = max(act10(g)[1] for g in K45)
pr(f"  Two link states combine (two-way) into the far side only (leak {spotpart:.0e}); the far side is reached in all {rk} of its parts.")
pr(f"  The 45 turns act consistently on the far side (error {res10:.0e}).")
r_c = [act10(C(c))[0] for c in cn]; r_w = [act10(C(c))[0] for c in wn]; r_y = act10(Yn)[0]
Cc10 = -sum(m@m for m in r_c); Cw10 = -sum(m@m for m in r_w); Y10 = -1j*r_y
wv, Uv = np.linalg.eig(Cw10 + math.pi*Cc10 + math.e*Y10); Ui = np.linalg.inv(Uv); tab10 = {}
for k_ in range(len(far)):
    key = (round(((Ui[k_]@Cc10@Uv[:, k_])).real, 2) + 0.0, round((Ui[k_]@Cw10@Uv[:, k_]).real, 2) + 0.0, round((Ui[k_]@Y10@Uv[:, k_]).real, 4) + 0.0); tab10[key] = tab10.get(key, 0) + 1
pr("  What the far side is, under the 12 carriers:")
for key, m in sorted(tab10.items()): pr(f"     colours {'3' if key[0] > 1 else '1'}, {'pair' if key[1] > .5 else 'single'}, charge {key[2]:+.4f}: {m} parts")
PH = np.eye(len(far)) - Cc10/(4/3)                                   # projector on the colourless part of the far side
kinds = {}
for a_ in range(16):
    key = (round(rows[a_][0], 2), round(rows[a_][1], 2), round(rows[a_][2], 4) + 0.0); kinds.setdefault(nm.get((abs(key[0]), abs(key[1]), key[2]), '?').lower(), []).append(a_)
names_ = ["quark pair", "lepton pair", "up partner", "down partner", "electron partner", "hidden neutrino partner"]
pr("  Strength of the link between two kinds of state through the COLOURLESS weak pair of the far side (0 = no link):")
pr("     " + " "*26 + "".join(f"{n_[:14]:>16}" for n_ in names_))
link = {}
for A in names_:
    line = f"     {A:<26}"
    for B_ in names_:
        s_ = 0.0
        for a_ in kinds[A]:
            for b_ in kinds[B_]:
                v_ = PH@np.einsum('kab,a,b->k', Mm, U[:, a_], U[:, b_]); s_ += np.vdot(v_, v_).real
        link[(A, B_)] = s_; line += f"{s_:>16.3f}"
    pr(line)
per = {B_: link[(A, B_)]/len(kinds[B_]) for A, B_ in [("quark pair", "up partner"), ("quark pair", "down partner"), ("lepton pair", "electron partner"), ("lepton pair", "hidden neutrino partner")]}
pr("  Per single state: " + "; ".join(f"{k_} {v_:.3f}" for k_, v_ in per.items()))
pr("  -> Every state has exactly one mass partner: a member of a weak pair links to a single state, through the colourless weak pair of charge 1/2")
pr("     that sits on the far side of the matrix (the Higgs field's numbers: 4 parts). No 'left' or 'right' was used. Mass = the link between a")
pr("     paired state and a single state, made through the part of the matrix away from the observer.")
pr("     ATTACK 1: the four links come out EQUAL in strength. Taken alone that says top = bottom = tau. Measured: top is 60 to 80 times bottom")
pr("     at every energy from the Z upward. So this gives WHICH links exist, not their sizes. Sizes stay where earlier tests left them.")
pr("     ATTACK 2: the far side also has 6 coloured parts (charge 1/3). If they existed as light particles the proton would decay. Known problem")
pr("     (Dimopoulos & Georgi 1981). In TRPG they are entries of the matrix, not particles, but that has to be said, not assumed.")

head("A4. So can 'handedness' be removed and replaced by the perspective? The reasoning, step by step")
pr("  1. Each of the 16 states is a spinning wave. Relativity's simplest spinning wave has ONE screw sense and is written with an i (Weyl 1929).")
pr("  2. Quantum mechanics uses one i for everything (Stueckelberg 1960; tested by Renou et al. 2021). An inherited rule, like relativity.")
pr("  3. The three-way matrix has exactly one i for all its states together, the observer's facing direction (A1).")
pr("  4. Same i for the wave and for the states: then all 16 states get the same screw sense relative to their charges. No state can be")
pr("     flipped separately (A1: one sign for all), and nothing is doubled.")
pr("  5. Facing the other way (i -> -i) then changes three names at once: matter <-> antimatter, left <-> right, i <-> -i. These are not three")
pr("     facts. And which sense is called 'left' is only a name: a world with every sense reversed is ours seen in a mirror.")
pr("  6. What is physical, and measured (Lee & Yang 1956; Wu 1957): the world is not its own mirror image. That follows from A2.")
pr("  7. ATTACK: give the wave its own separate i. Then every state comes in both screw senses: mirror states. A state and its mirror can")
pr("     always pair off with no Higgs, so the measured masses would have no reason; and rule 6 (nothing extra) forbids the doubling.")
pr("     ATTACK: use the 12 carriers alone, without the whole three-way structure. Then each kind of state could be flipped separately (A1: 14")
pr("     commuting matrices); most such choices fail the consistency sums. The whole structure allows no such choice.")
pr("  -> HOLDS. Nothing the main line has uses 'left' or 'right' as a separate idea: the 12 carriers, the charges, the consistency sums, the")
pr("     mass links, pion decay (37.6), neutrino masses and the matter excess all go through unchanged. The earlier added reading 'spin uses")
pr("     the same i as the states' is not removed: it becomes part of rule 4, 'one perspective, one i'. One reading fewer, same content.")
pr("     What stays inherited, not derived: that matter is made of spinning waves at all (step 1) and that there is one i (step 2).")

head("B1. Gravity's strength as a number")
hbarc = 0.1973269804; GN = 6.70883e-39                               # G in GeV^-2 (PDG 2024)
MPl = 1/math.sqrt(GN); mp = 0.93827208816; MPr = MPl/math.sqrt(8*math.pi)
pr(f"  Gravity has no strength 'by itself': G carries units. As a pure number it is a ratio. With the proton as the unit (the 'proton = 1'):")
pr(f"     grain energy (Planck) = {MPl:.4e} GeV = {MPl/mp:.3e} protons;  gravity between two protons, as a pure number: {GN*mp**2:.3e} (the electric force: 1/137).")
pr(f"     The honest size of the question is the logarithm: ln(grain / proton) = {math.log(MPl/mp):.2f}.")
pr("  In TRPG : entanglement across a surface grows with its area (Bombelli et al. 1986; Srednicki 1993), and that gives")
pr("     Einstein's equations with G = 1/(4 x entanglement per area) (Jacobson 1995). So gravity's strength IS the fineness of the links of space.")

head("B2. Is gravity's strength a separate number from the force strengths? (Two-loop running, as in an earlier test)")
MZ = 91.1876; a_em_inv = 127.951; s2 = 0.23122; a_s = 0.1180
a_inv0 = np.array([(3/5)*(1 - s2)*a_em_inv, s2*a_em_inv, 1/a_s]); b1 = np.array([41/10, -19/6, -7.0])
b2 = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12.0], [11/10, 9/2, -26.0]])
def run(t1, n_=4000):
    y = a_inv0.copy(); h = (t1 - 0)/n_
    f = lambda y: -b1/(2*math.pi) - b2@(1/y)/(8*math.pi**2)
    for _ in range(n_):
        k1 = f(y); k2 = f(y + h/2*k1); k3 = f(y + h/2*k2); k4 = f(y + h*k3); y = y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return y
lo, hi = math.log(1e10/MZ), math.log(1e19/MZ)
for _ in range(60):
    mid = (lo + hi)/2; y = run(mid, 800)
    if y[1] > y[2]: lo = mid
    else: hi = mid
tU = (lo + hi)/2; yU = run(tU); MU = MZ*math.exp(tU); yP = run(math.log(MPl/MZ))
pr(f"  Strong = weak at 10^{math.log10(MU):.2f} GeV, where 1/strength = {yU[2]:.1f} (hypercharge there: {yU[0]:.1f}).")
pr(f"  At the grain energy: 1/strength = {yP[0]:.1f} (hypercharge), {yP[1]:.1f} (weak), {yP[2]:.1f} (strong).")
L1 = 2*math.pi/(7*(1/yP[2])); pr(f"  The strong force grows slowly from there down to the proton: one-loop count of the logarithm = 2 pi/(7 x strength) = {L1:.1f},")
pr(f"     against ln(grain/proton) = {math.log(MPl/mp):.1f} (the difference is the heavy-quark steps and the proton's size in units of the strong scale).")
pr("  -> Known physics (Gross & Wilczek 1973; Politzer 1973; Wilczek 2001): the proton is light against the grain BECAUSE the strong force is")
pr(f"     weak at the grain (about 1/{yP[2]:.0f}) and strengthens only logarithmically. So gravity's strength and the force strengths are TIED: gravity's")
pr("     weakness is the same fact as the strong force's weakness at the grain of space. This is the proposed tie. It trades one number")
pr("     for another; the count of measured numbers does not drop, and nothing is derived yet.")
aU = 1/yU[2]; MPl_pred = MU/math.sqrt(aU)
pr(f"  TEST: 'gravity is as strong as the other forces where they meet each other'. Gravity's number at energy E is (E/grain)^2. At the meeting")
pr(f"     energy it is {(MU/MPl)**2:.1e}, against {aU:.4f} for the others: {aU/(MU/MPl)**2:,.0f} times too weak. Put the other way: this predicts a grain")
pr(f"     energy of 10^{math.log10(MPl_pred):.2f} GeV against the measured 10^{math.log10(MPl):.2f}. FAIL as a number ({MPl/MPl_pred:.0f}x in energy), though on the")
pr(f"     logarithmic scale it lands at {math.log(MPl_pred/mp)/math.log(MPl/mp)*100:.0f}% of the way from the proton to the grain. (Same known gap as Kaplunovsky 1988.)")
Estar = MPl*math.sqrt(1/yP[1]); pr(f"     Gravity actually reaches the others' strength at about 10^{math.log10(Estar):.2f} GeV, {Estar/MU:.0f} times above the meeting energy.")
Ndof = 3*16*2 + 12*2 + 4
pr(f"  TEST: 'gravity's strength from counting TRPG's own states' (Sakharov 1967; Dvali 2007). States: 3 x 16 matter x 2 + 12 carriers x 2 + 4 Higgs = {Ndof}.")
pr(f"     Grain = Planck/sqrt(N) = 10^{math.log10(MPl/math.sqrt(Ndof)):.2f} GeV to Planck x sqrt(12 pi/N) = 10^{math.log10(MPl*math.sqrt(12*math.pi/Ndof)):.2f} GeV, depending on how the cut is made.")
pr("     The count turns a grain SIZE into a strength; it does not give the size. No TRPG scale sits there. FAIL as a derivation.")

head("B3. Gravity and the Higgs field, in both directions")
v = 246.21965; xi = MPr**2/v**2
pr(f"  Direction 1, 'the Higgs field sets gravity's strength' (Zee 1979): needs a Higgs-to-curvature link of {xi:.1e}.")
pr(f"     Measured limit from Higgs rates at the LHC: below 2.6e15 (Atkins & Calmet 2013). FAIL by a factor {xi/2.6e15:.0e}.")
pr(f"     With an ordinary link (1/6), the Higgs field changes gravity's strength by {v**2/(6*MPr**2):.0e} of itself. Nothing.")
mt, mh = 172.57, 125.20; crit = 129.6 + 2.0*(mt - 173.34) - 0.5*(a_s - 0.1184)/0.0007; dcrit = 1.2     # error as in an earlier test (includes the top-mass theory error)
pr("  Direction 2, 'gravity's grain sets the Higgs' (Shaposhnikov & Wetterich 2010, who predicted 126 GeV two years before the discovery):")
pr(f"     the Higgs's hold must fade to zero exactly at the grain. With today's top mass that needs a Higgs of {crit:.1f} +/- {dcrit:.1f} GeV (Buttazzo et al. 2013);")
pr(f"     measured {mh}. Off by {(crit - mh)/mh*100:.1f}%, {(crit - mh)/dcrit:.1f} sigma (the same test as an earlier test, B4). A near miss, strict FAIL. The hold fades near")
pr("     10^10.5 GeV, the region of TRPG's hidden neutrino partners (about 10^9 to 10^11 GeV), not the grain. That stays a lead.")
rhoL4 = 2.25e-12; stored = 104.0
pr(f"  The Higgs and dark energy: the Higgs balance stores (104 GeV)^4. Dark energy is ({rhoL4*1e12:.2f} meV)^4. Ratio {(stored/rhoL4)**4:.1e}.")
pr("     If stored balance energy acted as dark energy the universe would not exist as it is (Weinberg 1989). In gravity-from-entanglement it does")
pr("     not: there the dark energy is a constant fixed once, not the energy of any field (Jacobson 1995; Einstein 1919). In TRPG words: it is a")
pr("     kept amount of the bag (rule 2), set at the pop. HOLDS as the form (and it is required, or the Higgs would wreck gravity). No value.")
pr("     Prediction of this form: dark energy never changes (w = -1 exactly). DESI 2025 hints at a change (2.8 to 4.2 sigma). Live test.")

head("B4. TRPG's dark matter: 'the gravity-only class'. Does the main line have a member?")
pr("  The only state with no force charge at all is the hidden neutrino partner (A2), the very state the observer keeps exactly.")
hb = 6.582119569e-25; age = 4.35e17; vv = 174.10
for nme, M, mnu in [("lightest (matter-excess fit)", 3.6e4, 0.36e-12), ("lightest", 7e5, 0.36e-12), ("middle", 1.3e9, 8.66e-12), ("heaviest", 6.3e10, 50.1e-12)]:
    Gam = M**2*mnu/(8*math.pi*vv**2); pr(f"     {nme:<29} {M:.1e} GeV: it decays through the same link that gives the neutrino its mass; lifetime about {hb/Gam:.0e} s or less")
pr(f"     Age of the universe: {age:.1e} s. FAIL by more than 30 powers of ten: they are gone in the first instant. A partner made stable would")
pr("     leave its neutrino massless, against TRPG's own lightest neutrino (0.36 meV).")
cH = 2.99792458e8*67.4e3/3.0857e22
pr(f"  The other reading, 'dark matter is entanglement's response, not particles' (Verlinde 2017): predicts an acceleration scale c H0/6 = {cH/6:.2e} m/s^2;")
pr("     galaxies show 1.20e-10 +/- 0.24e-10 (McGaugh, Lelli & Schombert 2016). Right size. But published tests of the detailed galaxy pattern, the")
pr("     Solar System and galaxy clusters go against it (Lelli, McGaugh & Schombert 2017; Hees, Famaey & Bertone 2017; Halenka & Miller 2020),")
pr("     and it has no account of the pattern in the cosmic background. Not carried forward.")
pr("  -> The main line has no dark matter candidate among its states. That is an open item the list did not show. So is dark energy's amount:")
pr("     one measured number that the input table (40.5) left out. Honest count.")

head("B5. Quark mixing: what the main line already gives, strictly")
ms_md = 27.42*(1 + 0.465)/2; dms = 27.42*0.024/2; pV = math.sqrt(1/ms_md); dpV = pV*0.5*dms/ms_md   # FLAG 2024
tests = [("|V_us| = sqrt(m_d/m_s)  (Gatto, Sartori & Tonin 1968)", pV, dpV, 0.2243, 0.0008),
         ("|V_td/V_ts| = sqrt(m_d/m_s)  (Fritzsch 1977)", pV, dpV, 0.207, 0.003),
         ("|V_ub/V_cb| = sqrt(m_u/m_c)  (Fritzsch 1977)", math.sqrt(2.16e-3/1.10), 0.001, 3.82e-3/41.1e-3, 0.0056)]
for nme, p_, dp, m_, dm in tests: pr(f"     {nme:<52} predicted {p_:.4f}, measured {m_:.4f} +/- {dm:.4f}: {(p_ - m_)/math.hypot(dp, dm):+.1f} sigma")
pr("  The Cabibbo number and the down mass are already tied in TRPG (one follows from the other, an earlier test). The same rule on the other two fails.")
pr("  A3 says the simplest mass link gives NO mixing at all; measured mixing is small (13.0, 2.4 and 0.2 degrees), so that is the right starting")
pr("  point, and what is open is the small tilt between the up and down placements. It did not come out of gravity or the Higgs here.")
