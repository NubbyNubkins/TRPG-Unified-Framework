"""Left-handedness in the relationship-matrix picture, and the spiral picture (14 Nov 2026).
Question: once the weak force acts on pairs, do TRPG's rules force every pair to have ONE handedness?
the picture (14 Nov): a particle/wave is one spiral geometry, and handedness is the side you observe it from.
Uses the matrix tools of forces_from_entanglement_2026-11-14.py (same folder).
Credit: Weyl 1929; Lee & Yang 1956; Wu et al. 1957; Garwin, Lederman & Weinrich 1957; Frauenfelder et al. 1957; Goldhaber, Grodzins & Sunyar 1958;
Feynman & Gell-Mann 1958; Sudarshan & Marshak 1958; Landau 1957 (CP); Wigner; Majorana 1937; Stueckelberg 1960 (the i of quantum mechanics);
Georgi 1979; Krasnov 2019 (the two complex structures on the octonion pair); Todorov & Dubois-Violette 2018; Furey 2018; Dixon; Baez 2020-21;
Pati & Salam 1974; Bouchiat, Iliopoulos & Meyer 1972; Gross & Jackiw 1972; Witten 1982; TWIST 2011; PIENU 2015; ALEPH/DELPHI/L3/OPAL (LEP); PDG 2024."""
import math, itertools, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "forces_from_entanglement_2026-11-14.py")).read()
exec(src[src.index("import math, itertools"):src.index('head("1. Entanglement')])          # pr, head, mul, conj, restrict, commutant, derived, span_basis, rank_of
exec(src[src.index("def jordan_algebra(n):"):src.index("names = {2:")])                      # jordan_algebra, derivations

n = 8; dim, T = jordan_algebra(n); der = derivations(dim, T)
W = [0, 1, 2] + [3 + b*n + c for b in range(3) for c in (0, 1)]; notW = [i for i in range(dim) if i not in W]
face = restrict(der, lambda g: g[np.ix_(notW, W)].ravel()); spot = restrict(der, lambda g: g[:, 2]); both = restrict(face, lambda g: g[:, 2])
colour = restrict(both, lambda g: g[:, W].ravel()); weak4 = commutant(both, colour); cen = commutant(both, both)
links = list(range(3, 3 + 2*n)); V = lambda g: g[np.ix_(links, links)]
Gb = [V(g) for g in both]; Gs = [V(g) for g in spot]
E = np.eye(8); u = E[1]                                              # u: the direction the observer faces
Lu = np.array([mul(u, E[j]) for j in range(8)]).T; Ru = np.array([mul(E[j], u) for j in range(8)]).T
def blk(a, b):
    Z = np.zeros((16, 16)); Z[:8, :8] = a; Z[8:, 8:] = b; return Z
JL = blk(Lu, -Ru)      # 'multiply the link column by u from the left'  (the column is (conj(x2), x1))
JR = blk(Ru, -Lu)      # 'multiply the link column by u from the right'
Y = V(cen[0])
Cc = -sum(V(c)@V(c) for c in span_basis(colour)); w_, vec = np.linalg.eigh((Cc + Cc.T)/2)
Pl = vec[:, np.abs(w_) < 1e-8]@vec[:, np.abs(w_) < 1e-8].T; Pq = np.eye(16) - Pl
Yn = Y/math.sqrt(-np.trace(Pq@Y@Y)/12)/6
def charges(J): return -np.trace(J@Pq@Yn)/12, -np.trace(J@Pl@Yn)/4

head("1. The link states have only one 'i' of their own: the direction the observer faces")
basis = []
for k in range(256):
    b_ = np.zeros((16, 16)); b_.flat[k] = 1; basis.append(b_)
cm = restrict(basis, lambda g: np.concatenate([(g@s - s@g).ravel() for s in Gb]))
pr(f"  Everything that commutes with all 12 force carriers on the 16 link components: a {len(cm)}-dimensional family.")
for nm, J in (("left-multiply by u ", JL), ("right-multiply by u", JR)):
    pr(f"    {nm}: squares to -1: {np.allclose(J@J, -np.eye(16))}; commutes with the 12: {max(np.abs(g@J - J@g).max() for g in Gb) < 1e-9}; inside that family: {len(span_basis(cm + [J])) == len(cm)}")
pr(f"  That family is spanned by: the identity on the quark part, the identity on the lepton part, and these two. So every 'i' the link states")
pr("  can use without disturbing the forces is multiplication by the facing direction u, from the left or from the right. There is no other.")

head("2. Spin needs an 'i'. Two ways to supply it")
pr("  A spinning particle is described by complex numbers (Weyl 1929). Its 'i' is either")
pr("   (a) brought in from OUTSIDE the octonions, or (b) the link states' own i, the facing direction u.")
I16 = np.eye(16)
for nm, J in (("left", JL), ("right", JR)):
    Pp, Pm = (I16 - 1j*J)/2, (I16 + 1j*J)/2
    pr(f"  (a) with an outside i and u acting from the {nm}: (1 - iu)/2 and (1 + iu)/2 are both non-zero (sizes {np.linalg.matrix_rank(Pp)}, {np.linalg.matrix_rank(Pm)}) and multiply to {np.abs(Pp@Pm).max():.0e}.")
pr("      Two non-zero things multiplying to nothing: a hole (rule 1), the same hole that ruled out the joined number systems in an earlier test.")
pr("      The two halves are the pairs and their MIRROR pairs (opposite handedness, same charges). Option (a) = both handednesses = the hole.")
pr("  (b) with the link states' own i there is no second half: 8 one-handed spinning states, not 16. One handedness only.")
pr("  -> If rule 1 applies to the numbers a spinning particle is written in, the mirror particles ARE the hole, and matter is one-handed.")

head("3. Which one-handed arrangement? There are exactly two, and one is inconsistent")
pr(f"  {'the i is':22s} {'quark pair charge':>18s} {'lepton pair charge':>19s} {'weak-weak-charge sum (must be 0)':>34s} {'turns of the spot that respect it':>36s}")
for nm, J in (("u from the left", JL), ("u from the right", JR)):
    q, l = charges(J); keep = restrict(Gs, lambda g: (g@J - J@g).ravel())
    pr(f"  {nm:22s} {q:+18.4f} {l:+19.4f} {3*q + l:+34.3f} {len(keep):>30d} (rank {rank_of(keep)})")
keepL = restrict(Gs, lambda g: (g@JL - JL@g).ravel()); dL = derived(keepL)
pr(f"  Left: quark pair +1/6 and lepton pair -1/2. The sum vanishes: consistent. This is the measured left-handed world.")
pr(f"  Right: both pairs the same sign. The sum is -1: such a world is mathematically inconsistent (Bouchiat, Iliopoulos & Meyer 1972). It is the")
pr("  arrangement in which the quark pair and the lepton pair have OPPOSITE handedness.")
pr(f"  -> Consistency leaves one arrangement: all weak pairs the same handedness, charges +1/6 and -1/2. Nothing was chosen to get this.")
pr(f"  Side note: the turns that respect 'u from the left' number {len(keepL)} (self-interacting part {len(dL)}): four colours x the weak force (Pati & Salam's lepton as")
pr("  a fourth colour, an earlier test) WITHOUT a right-handed force. Facing one way then cuts the four colours to three plus the lepton, giving the 12.")
pr("  The turns that respect 'u from the right' are exactly the 12 (Krasnov 2019). So the two i's see different things: one fixes the forces, the other the matter.")

head("4. the spiral picture: handedness as the side you observe from")
q1, l1 = charges(JL); q2, l2 = charges(-JL)
pr(f"  Face the other way (u -> -u). The same link states now read: quark pair {q2:+.4f}, lepton pair {l2:+.4f} (were {q1:+.4f}, {l1:+.4f}).")
pr("  Every charge flips sign. A left-handed particle with charge +q is the same object as a right-handed antiparticle with charge -q, so an observer")
pr("  facing the other way describes the SAME world as right-handed, with 'matter' and 'antimatter' swapped (Landau 1957).")
pr("  -> In this picture left-handedness is not a property the particles have on their own. It is which way the observer faces.")
pr("")
pr("  One correction to the spiral as drawn. A spiral moving along its axis:")
om, v = 1.0, 0.6
for nm, look, V_obs in (("observer behind, at rest ", +1, 0.0), ("observer in front, at rest", -1, 0.0), ("observer who overtakes it ", +1, 0.9)):
    vz = (v - V_obs)/(1 - v*V_obs)                 # speed seen by the observer (relativistic)
    sense = om*look                                 # +: clockwise as seen along the line of sight
    hel = math.copysign(1, om*vz)                   # turning sense about the direction of motion
    pr(f"    {nm}: sees it {'receding' if vz*look > 0 else 'approaching'}, turning {'clockwise' if sense > 0 else 'counter-clockwise'}; handedness of the spiral = {'right' if hel > 0 else 'left'}")
pr("  The two observers at rest disagree about clockwise, but they agree on the spiral's handedness: a screw is right-handed from")
pr("  either end. What does flip the handedness is OVERTAKING the wave. That is possible only for a wave slower than light, that is, one with mass.")

head("5. The spiral picture tested against measurement: mass is what shows the other side")
pr("  A left-handed spinning state moving at speed v is seen as a left-turning spiral with probability (1 + v/c)/2.")
for nm, beta in (("neutrino (mass below 1 eV)", 1 - 1e-15), ("beta-decay electron, 1 MeV", math.sqrt(1 - (0.511/1.511)**2)), ("beta-decay electron, 100 keV", math.sqrt(1 - (0.511/0.611)**2)), ("slow electron, v = 0.1c", 0.1)):
    pr(f"    {nm:30s}: left-turning {100*(1+beta)/2:6.2f}%, net polarisation {-beta:+.3f}")
pr("  Measured: beta-decay electrons have polarisation -v/c (Frauenfelder et al. 1957 and many since); the neutrino's is -1 (Goldhaber, Grodzins & Sunyar 1958).")
me, mmu, mpi = 0.51099895, 105.6583755, 139.57039
R0 = (me/mmu)**2*((mpi**2 - me**2)/(mpi**2 - mmu**2))**2
pr(f"  Sharper test. A pion's decay to an electron needs the electron's 'other side', which only its mass supplies. Predicted rate relative to the muon:")
pr(f"    (m_e/m_mu)^2 x phase space = {R0:.4e}; with the known radiative correction (-3.76%) {R0*(1-0.0376):.4e}. Measured (PIENU 2015): 1.2344e-04 +/- 0.0030e-04.")
pr("    Without the handedness rule the electron mode would be about 5 times MORE common than the muon mode, not 8,000 times rarer.")
pr("  -> The picture holds quantitatively: the weak force sees one side of the spiral; mass is the link to the other side.")

head("6. What links the two sides: the Higgs balance point of an earlier test")
yQ, yL, yu, yd, ye = 1/6, -1/2, 2/3, -1/3, -1
pr(f"  Mass must join a pair member to its unpaired partner. The joining field needs charge: down {yQ - yd:+.3f}, electron {yL - ye:+.3f}, up {-(yQ - yu):+.3f}.")
pr("  All three are the same: ONE field, a weak pair of charge 1/2, gives every particle its mass. That is the Higgs, with the same charges as the lepton pair.")
pr("  Nothing extra is needed (rule 6). And because every pair is one-handed, no mass is possible WITHOUT that field: no mass 'by magic'.")

head("7. Attacks")
rest = [0, 1] + list(range(3 + 2*n, 3 + 3*n))
Yr = V(cen[0])*0; Yfull = cen[0]/math.sqrt(-np.trace(Pq@Y@Y)/12)/6
yv = np.round(np.sqrt(np.abs(np.linalg.eigvalsh(Yfull[np.ix_(rest, rest)]@Yfull[np.ix_(rest, rest)]))), 4); vals, counts = np.unique(yv, return_counts=True)
pr(f"7a. Are the partners in the matrix? Charges available outside the pairs: {vals.tolist()} (multiplicities {counts.tolist()}).")
pr("    Needed: 2/3 (up partner), 1/3 (down partner), 1 (electron partner), 0 (hidden neutrino partner). The matrix has room for the 1/3 and the 0 only.")
pr("    -> FAILS as 'one matrix holds a whole generation'. The up-type and electron-type partners are forced by consistency but are not in the matrix.")
pr("7b. Is the step 'the i of spin is the facing direction' forced? Only if rule 1 applies to the numbers a spinning particle is written in.")
pr("    A fully real description (Majorana 1937) avoids any outside i, but it then contains both this i and spin's own turning as separate things,")
pr("    which is option (a) again. So the step is one reading of rule 1, the same reading an earlier test used against the joined number systems. Not a theorem.")
pr("7c. Left or right from the left? Fixed by consistency (step 3), not by a rule. Two options existed; one is impossible.")
pr("7d. Data. The weak force is one-handed as far as measured:")
pr("    muon decay shape parameter rho = 0.74977 +/- 0.00026 (TWIST 2011); purely one-handed predicts 0.75000. Difference "
   f"{(0.74977-0.75)/0.00026:+.1f} sigma.")
pr("    light neutrino kinds from the Z's width: 2.9963 +/- 0.0074 (LEP); mirror pairs lighter than 45 GeV are excluded.")
pr("    no right-handed carrier found below about 5 TeV (LHC). Any right-handed weak effect, or any mirror particle, would break this result.")
