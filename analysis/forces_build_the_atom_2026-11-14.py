"""Forces first, the atom as the result. The frame: two kinds of thing, a push/pull that sets the energy
and a change of perspective; the strong and electric forces bind, the weak force changes what a particle is; a third force with a third facing may sit
on the Higgs, where the 'charge' is mass and the families differ; this might answer the three generations. Each statement is tested here and attacked.
Uses perspective_and_gravity_2026-11-14.py and three_way_relationship_2026-11-14.py (same folder).
Credit: Lee & Yang 1956; Wu 1957; Georgi 1979 (survival hypothesis); Vafa & Witten 1984; Michel & Radicati 1971; Higgs 1964; Englert & Brout 1964;
Weinberg 1967; Salam 1968; Glashow 1961; Georgi & Glashow 1974; Georgi 1975; Fritzsch & Minkowski 1975; T. D. Lee 1973; Gross & Wilczek 1973;
Politzer 1973; Yukawa 1935; Shifman, Vainshtein & Zakharov 1978; Bressi et al. 2011; MICROSCOPE 2022; Cartan 1935; Freudenthal 1953;
Jordan, von Neumann & Wigner 1934; Albert 1934; Loos 1977; Koide 1982; Brannen 2006; Kobayashi & Maskawa 1973; PDG 2024."""
import os, io, contextlib, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src4 = open(os.path.join(HERE, "perspective_and_gravity_2026-11-14.py")).read()
with contextlib.redirect_stdout(io.StringIO()):                       # the three-way product, its 12 turns, the far side and the mass map
    exec(src4[src4.index("import os, io"):src4.index('head("A4.')])

col8 = [C(c) for c in cn]; wk3 = [C(c) for c in wn]; Yg = Yn; gens12 = col8 + wk3 + [Yg]
def facing_blind(gens):
    """True if the 16 states, under these turns, are their own mirror image (an invertible S with S conj(g) = g S for every turn g)."""
    I = np.eye(16); Gm = np.zeros((256, 256), complex)
    for g in gens:
        M = np.kron(I, g.conj().T) - np.kron(g, I); Gm += M.conj().T@M
    w, Vv = np.linalg.eigh(Gm); k = int((w < 1e-9*max(w[-1], 1e-30)).sum())
    if k == 0: return False
    best = 0.0
    for _ in range(6):
        S = (Vv[:, :k]@(rng.normal(size=k) + 1j*rng.normal(size=k))).reshape(16, 16); sv = np.linalg.svd(S, compute_uv=False); best = max(best, sv[-1]/sv[0])
    return best > 1e-6

head("1. 'Blind to the facing': which of the 12 turns treat the world and its mirror image alike?")
t3 = wk3[0]; chk3 = np.round(np.linalg.eigvalsh(-1j*t3), 3); pr(f"  (check: one weak turn has values {sorted(set((chk3 + 0.0).tolist()))} on the states)")
cases = [("all 12 (strong + weak + the single one)", gens12), ("strong + weak (11)", col8 + wk3), ("strong + the single one as it stands (9)", col8 + [Yg]),
         ("strong + one weak turn + the single one (10)", col8 + [t3, Yg]), ("strong + one weak turn alone (9)", col8 + [t3]), ("weak + the single one (4)", wk3 + [Yg]),
         ("strong alone (8)", col8), ("weak alone (3)", wk3)]
for nme, gs in cases: pr(f"     {nme:<48} blind to the facing: {facing_blind(gs)}")
good = [a for a in np.round(np.arange(-3, 3.0001, 0.125), 3) if facing_blind(col8 + [Yg + a*t3])]
pr(f"  Strong + ONE mixture (the single one + a x a weak turn), a scanned from -3 to 3 in steps of 1/8: blind to the facing only for a = {[float(a_) for a_ in good]}.")
idx = kinds["down partner"]; P3 = U[:, idx]; A3 = [P3.conj().T@g@P3 for g in col8]
Ms = np.array([np.concatenate([a_.imag.ravel()]) for a_ in A3]).T; un, sn, vn = np.linalg.svd(Ms); nul = vn[(sn > 1e-9).sum():]
so3 = [sum(c*g for c, g in zip(cv, col8)) for cv in nul]
pr(f"  ATTACK: smaller sets. A 3-turn part of the strong force (turning the colours like an ordinary rotation: {len(so3)} turns) + all of weak: {facing_blind(so3 + wk3)} (6 turns).")
pr("  -> Subsets of the 12 with 9 or more turns: all 12; strong + weak (11); strong + two singles (10); strong + one single (9, a one-number family).")
pr("     Nothing else reaches 9, because the biggest proper part of the strong force has 4 turns. Of these, exactly ONE is blind to the facing:")
pr("     the strong force + the mixture with a = +1 or -1 (the same mixture seen from the two members of the weak pair).")
pr("     So 'the LARGEST part of the forces that is blind to the facing' is unique. Smaller blind sets exist (attack above), so 'largest' is needed.")

head("2. Let the far side of the matrix settle. What is kept, and is it the same set?")
keys10 = [(round((Ui[k_]@Cc10@Uv[:, k_]).real, 2) + 0.0, round((Ui[k_]@Cw10@Uv[:, k_]).real, 2) + 0.0, round((Ui[k_]@Y10@Uv[:, k_]).real, 4) + 0.0) for k_ in range(len(far))]
rho12 = [act10(g)[0] for g in gens12]
def kept_by(h):
    A = np.array([np.concatenate([(r_@h).real, (r_@h).imag]) for r_ in rho12]).T; u_, s_, v_ = np.linalg.svd(A); return v_[(s_ > 1e-9*s_[0]).sum():]
h0 = Uv[:, [k_ for k_ in range(len(far)) if keys10[k_] == (0.0, 0.75, 0.5)][0]]; hc = Uv[:, [k_ for k_ in range(len(far)) if keys10[k_][0] > 1][0]]
N0 = kept_by(h0); Nc = kept_by(hc); kept0 = [sum(c*g for c, g in zip(cv, gens12)) for cv in N0]; keptc = [sum(c*g for c, g in zip(cv, gens12)) for cv in Nc]
colour_kept = max(np.linalg.norm(r_@h0) for r_ in rho12[:8])
u_, s_, v_ = np.linalg.svd(np.array([np.concatenate([(r_@h0).real, (r_@h0).imag]) for r_ in rho12[8:]]).T); qc = v_[-1]; qc = qc/qc[3]
Qg = sum(c*g for c, g in zip(qc, wk3 + [Yg])); a_len = np.linalg.norm(qc[:3])
pr(f"  Settle the colourless weak pair (the Higgs numbers, 44.2): turns kept = {len(N0)}: all 8 strong ones (change {colour_kept:.0e}) and one single one,")
pr(f"     which is the earlier single one + {a_len:.3f} x a weak turn. Blind to the facing: {facing_blind(kept0)}. The 3 turns lost become the short-reach")
pr("     weak carriers (W+, W-, Z).")
qs = {}
for nme in names_:
    Pk = U[:, kinds[nme]]; qs[nme] = sorted(np.linalg.eigvalsh(Pk.conj().T@(-1j*Qg)@Pk).tolist())
for nme in names_: pr(f"     {nme:<26} value of the kept single turn: {sorted(set(round(q_, 4) + 0.0 for q_ in qs[nme]))}")
pr("     These are the electric charges: +2/3 and -1/3 for the quarks, 0 and -1 for the neutrino and electron, and the opposite for each partner.")
pr(f"  Settle a COLOURED part of the far side instead: turns kept = {len(Nc)}; blind to the facing: {facing_blind(keptc)}. (Colour would be broken.)")
pr("  -> The colourless direction keeps the larger set (9 against 7), and that set is exactly the one found in test 1:")
pr("     what survives the Higgs = the largest facing-blind part = the strong force + electric charge. These are the two forces that bind.")
pr("     HONEST NOTE: this is not two independent confirmations. If every state gets a mass partner, the turns that are kept MUST be blind to")
pr("     the facing, because the mass link itself maps the world onto its mirror image (Georgi 1979). So 'blind to the facing' and 'lets")
pr("     everything have mass' are one statement. The independent content is the converse: the largest blind set is one a settled far side keeps.")
pr("     Not derived: THAT the far side settles, and that it takes the direction keeping the most (a rule of thumb: Michel & Radicati 1971).")
S9 = np.concatenate([r_.T for r_ in [act10(g)[0] for g in kept0]], axis=0); u_, s_, v_ = np.linalg.svd(S9); ell = v_[(s_ > 1e-9*s_[0]).sum():].conj()
Rm = ell@Y10@np.linalg.pinv(ell); wv_, vv_ = np.linalg.eig(Rm.T); ells = [(vv_[:, i]@ell, wv_[i].real) for i in range(len(wv_))]
pr(f"  Settled directions that the 9 kept turns leave alone: {len(ell)} (charges {sorted(float(round(e_[1], 3)) + 0.0 for e_ in ells)}).")
for el, yv in ells:
    Mh = np.einsum('k,kab->ab', el, Mm); inv = max(np.abs(g.T@Mh + Mh@g).max() for g in kept0)/np.abs(Mh).max(); sv = np.linalg.svd(Mh, compute_uv=False); sv = sv/sv[0]
    blk = {(A, B_): np.linalg.norm(U[:, kinds[A]].T@Mh@U[:, kinds[B_]]) for A in names_ for B_ in names_}; tot = math.sqrt(sum(v_**2 for v_ in blk.values()))
    on = sorted({tuple(sorted(k_)) for k_, v_ in blk.items() if v_ > 1e-6*tot})
    pr(f"     direction of charge {yv:+.1f}: links {' and '.join(a_ + ' with ' + b_ for a_, b_ in on)}; {int((sv > 1e-6).sum())} states linked, sizes all equal: {bool(np.allclose(sv[sv > 1e-6], 1, atol=1e-6))}; unchanged by the 9 (error {inv:.0e})")
pr("  -> One settled direction links the up quark and the neutrino to their partners; the other links the down quark and the electron to theirs.")
pr("     If the far side is REAL (one Higgs pair, as 37.6 and 43.2 say), the two are the same size: top = bottom = tau (FAILS, 44.2).")
pr("     If the far side is COMPLEX, as the states of the three-way reading are, there are two Higgs pairs and top : bottom is free (T. D. Lee 1973),")
pr("     but bottom = tau stays (measured 0.6 to 0.73 at high energy, earlier tests: FAILS), and a second Higgs pair has not been seen.")
pr("     Either way the SIZES of the links are not given. Which states are linked is.")

head("3. The atom from the charges")
qu, qd = max(qs["quark pair"]), min(qs["quark pair"]); qe = min(qs["lepton pair"]); qn = max(qs["lepton pair"])
pr(f"  Proton = up + up + down: charge {2*qu + qd:+.12f}.  Neutron = up + down + down: {qu + 2*qd:+.12f}.  Electron: {qe:+.12f}.  Neutrino: {qn:+.12f}.")
pr(f"  Hydrogen atom = proton + electron: charge {2*qu + qd + qe:+.1e}. Measured: zero to 1 part in 10^21 (Bressi et al. 2011).")
pr("  Three quarks, one of each colour, are colour-neutral because there are exactly 3 colours. So the proton is 'white' and the atom is")
pr("  electrically neutral, both exactly, and both follow from the matrix. the 'neutral bundle' at the two scales is derived, not assumed.")

head("4. Locked in, hidden, or free: why the three charges are treated three ways")
nH = 1; b3 = -11 + (2/3)*kC*3; b2_ = -22/3 + (2/3)*kT*3 + (1/3)*0.5*nH; b1_ = (3/5)*((2/3)*kY*3 + (1/3)*2*0.25*nH)
pr(f"  From the main line's own count (3 families of 16 states, 12 carriers, one Higgs pair): how each strength changes with distance:")
pr(f"     strong {b3:+.3f} (grows at large distance); weak {b2_:+.3f} (would grow too); the single one {b1_:+.3f} (fades at large distance).")
L3 = 91.1876*math.exp(-2*math.pi*(1/0.1180)/(23/3)); L2 = 91.1876*math.exp(-2*math.pi*29.585/abs(b2_))
pr(f"  Strong: becomes unbreakable near {L3:.2f} GeV by the simple count (0.2 to 0.3 GeV in the full calculation): the proton's own scale. Colour is LOCKED IN.")
pr(f"  Weak: would lock in only near {L2:.0e} GeV. The Higgs side settles at 246 GeV, {math.log10(246/L2):.0f} powers of ten earlier. Its label is HIDDEN first.")
pr("  Electric: its strength fades with distance and nothing hides it (test 2). It runs FREE.")
pr(f"  Limits from the same count: colour stays locked only with {int((11*3/2)/kC)} or fewer families.")
pr("  -> The direction (grow or fade) follows from the count. Where each one happens uses the measured strengths. HOLDS, with that input.")

head("5. The Higgs link as the third force: its charge is mass")
v = 246.21965; mH = 125.20; hbarc_fm = 0.1973269804; reach = hbarc_fm/mH
pr(f"  Reach = {reach:.2e} fm = 1/{0.8409/reach:.0f} of the proton's radius. At the size of an atom it is cut by a factor e^-{52917.7/reach:.1e}. No leftover outside.")
crit = 0.8399
for nme, m in [("electron", 0.000511), ("top quark", 172.57)]:
    pr(f"     pull between two of the {nme}: strength {m**2/(4*math.pi*v**2):.1e} (electric: {1/137.036:.1e})")
mb = (crit*8*math.pi*v**2*mH)**(1/3); mw = 2*crit*80.377/(1/29.585)
pr(f"  Can it bind? A short-reach pull binds two equal particles only if strength x (half their mass) / carrier mass > {crit} (Yukawa 1935).")
pr(f"     The Higgs link would bind particles heavier than about {mb:.0f} GeV; the heaviest, the top, has {172.57/mb*100:.0f}% of that mass.")
pr(f"     The weak force would need about {mw:.0f} GeV. Neither binds anything that exists.")
pr("     So the pattern holds: the two facing-blind forces bind; the weak force and the Higgs link only change things (identity; mass).")
pr("     ATTACK: the reason they do not bind is their short reach, and the short reach is the Higgs side settling. So 'sees the facing' and")
pr("     'cannot bind' are tied through test 2, not independently.")
fl = 0.09; fN = fl + (2/9)*(1 - fl)
pr(f"  Is mass-as-charge the same as gravity's charge? Gravity: every kind of matter falls alike to 1.5 parts in 10^15 (MICROSCOPE 2022).")
pr(f"     Higgs charge per unit of mass: electron 1.00; proton about {fN:.2f} (Shifman, Vainshtein & Zakharov 1978), because 99% of the proton is balance")
pr("     energy that the Higgs link does not see. FAIL: the Higgs link is not gravity. They share one thing only: a single sign, so nothing cancels.")
pr("     Coincidence, NOT a result: the heavy quarks' share of the proton's Higgs charge is exactly 2/9 at first order; TRPG's twist is 2/9.")
pr("     One comes from counting colours and light quarks (2 x 3 / 27), the other is an angle. No link is known. Recorded so it is not rediscovered.")

head("6. The third facing and the three families")
x = rng.normal(size=dim) + 1j*rng.normal(size=dim); Lxx = np.array([tri27(x, x, np.eye(dim, dtype=complex)[k_]) for k_ in range(dim)]).T
ev = np.sort(np.linalg.eigvals(Lxx).real); cl = [[ev[0]]]
for e_ in ev[1:]:
    (cl[-1].append(e_) if abs(e_ - cl[-1][-1]) < 1e-7*abs(ev[-1]) else cl.append([e_]))
mult = [len(c_) for c_ in cl]; val = [float(np.mean(c_)) for c_ in cl]; sing = [v_ for v_, m_ in zip(val, mult) if m_ == 1]; octs = [v_ for v_, m_ in zip(val, mult) if m_ == 8]
pr(f"  The three-way product on the WHOLE matrix (27 parts). For a random element x, the operation 'x, x, .' has values with multiplicities {mult}.")
if len(sing) == 3 and len(octs) == 3:
    pairs_ = sorted((sing[i] + sing[j])/2 for i in range(3) for j in range(i + 1, 3)); pr(f"     Three single values {np.round(sing, 4).tolist()} and three eightfold ones {np.round(octs, 4).tolist()};")
    pr(f"     each eightfold value is the average of two single ones (error {max(abs(a_ - b_) for a_, b_ in zip(pairs_, sorted(octs))):.0e}).")
xs = rc(); Ls = L(xs, xs); evs = np.sort(np.linalg.eigvals(Ls).real); cls = [[evs[0]]]
for e_ in evs[1:]:
    (cls[-1].append(e_) if abs(e_ - cls[-1][-1]) < 1e-7*abs(evs[-1]) else cls.append([e_]))
pr(f"     The same on the 16 link states alone: multiplicities {[len(c_) for c_ in cls]} (two sizes).")
pr("  -> A relationship matrix with octonion entries, read three-way, has exactly THREE independent sizes: never four, because no larger such matrix")
pr("     exists (Jordan, von Neumann & Wigner 1934; Albert 1934; Cartan 1935; Loos 1977). From one spot the observer sees two of them.")
th = 2/9; sz = sorted(1 + math.sqrt(2)*math.cos(th + 2*math.pi*k_/3) for k_ in range(3)); me, mm, mt_ = 0.51099895, 105.6583755, 1776.93
pr(f"  TRPG's own rule 5 already has this shape: the three lepton masses are the three sizes of ONE 3x3 matrix, at three angles 120 degrees apart")
pr(f"     (the 'third facing' as an angle): sizes {np.round(sz, 4).tolist()}; squared ratios {sz[0]**2/sz[1]**2:.6f} and {sz[1]**2/sz[2]**2:.6f};")
pr(f"     measured electron/muon {me/mm:.6f}, muon/tau {mm/mt_:.6f}. (Koide 1982; Brannen 2006; an earlier test. The twist 2/9 is measured, an earlier test.)")
pr("  TEST: 'the three families are the three spots of the one matrix' (the observer's third choice). Each spot sees 16 link states, so three")
pr(f"     families need 48. The matrix has 3 links of 8, each written with an i: {3*8}. Each link would be shared by two families. FAIL by a factor 2.")
pr("     (From any one spot the third link is not matter at all: it is on the far side, with the Higgs numbers.)")
pr("  TEST: 'the Higgs link tells the families apart' in the main line as it stands (three placements, an earlier test). The 12 turns and the three-way")
pr("     product act the same on every placement and put no condition on how the link differs between them. It is free: that is where the")
pr("     measured rules 5 and 7, the heavy quark masses, the mixing tilt and the twist all sit. NOT DERIVED.")
pr("  -> 'No more than three': the matrix gives 'at most three sizes' as a theorem. It becomes 'at most three families' only if a family IS a size")
pr("     of one matrix, and one matrix cannot hold both the forces and the families (24 against 48). So it is a strong hint that the 3 of the")
pr("     families and the 3 of the matrix are the same 3, and it is NOT a derivation. The item stays open, resting on rule 6 and data.")

head("7. The switch tested: 'the Higgs counts all energy because it is a constant; gravity supplies much more of the mass'")
pr("  The idea behind it, taken in two halves and tested against measurement.")
pr("  HALF 1: the Higgs field is a constant filling space (true: 246 GeV everywhere), so it should count ALL energy.")
nl = 5; b_light = 11 - 2*nl/3; amp = (b_light/8)/(1/12)
pr(f"     Test: most of a proton's energy is gluon energy. A field that counted all of it would link to gluon energy {amp:.1f} times more strongly")
pr(f"     than the Higgs is predicted to (which sees gluons only through the top quark). The LHC makes Higgs particles mostly from gluons, so")
pr(f"     the rate would be {amp**2:.0f} times higher. Measured rate / standard prediction: 1.0 +/- 0.1 (ATLAS and CMS 2022). FAIL by a factor of about {amp**2:.0f}.")
pr("     (Georgi, Glashow, Machacek & Nanopoulos 1978; Shifman, Vainshtein & Zakharov 1978; Goldberger, Grinstein & Skiba 2008.)")
pr("     Second test: switch off the Higgs link to the up and down quarks, and the proton still weighs about 880 of its 938 MeV (94%).")
pr("     The proton's mass computed from the strong force alone agrees with measurement to about 2% (Duerr et al. 2008).")
pr(f"     Correction to my own '1%': the Higgs's share of the proton is 1% (up and down masses), about 9% (with the strange quark), or {fN*100:.0f}%")
pr("     if the heavy quarks' indirect effect is counted. So 1% undersold it. It is still far from all.")
pr("  HALF 2: gravity supplies much more of the mass.")
aG = 6.70883e-39*0.93827208816**2; Eg = aG*hbarc_fm/0.8409
pr(f"     As a force, no: the proton's own gravity holds {Eg:.1e} GeV, which is {Eg/0.938272:.0e} of its mass.")
pr("     As a statement about where the mass comes from, it holds as a READING in TRPG's own terms: 99% of the proton is balance energy, which")
pr("     is relationship energy (an earlier test; an earlier test), and gravity is the byproduct of relationship (entanglement). The proton's mass is")
pr("     defined in physics as exactly what gravity counts (the trace of its energy and stress). So most of the mass is relationship energy and")
pr("     only a small part comes from the Higgs link. NOT SHOWN: that the entanglement which gives gravity and the colour balance which gives")
pr("     the proton's mass are the same relationship. That would be a real result if it could be shown, and it is not shown here.")
pr("  WHAT IS CONSTANT AND COUNTS EVERYTHING: not the Higgs field but the kept amount of the bag (dark energy, 44.5), and gravity itself.")
pr("  -> The switch as stated FAILS on its first half (measured, a factor of about 130). Its second half holds as a reading, not as a calculation.")
pr("     Restated: gravity counts all the energy; most of that energy is relationship energy; the Higgs link adds between 1% and 29% of a proton.")
