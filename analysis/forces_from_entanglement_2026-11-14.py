"""Forces as byproducts of relationship / entanglement, the way gravity is.
Question: after an earlier test, can the weak force, hypercharge and left-handedness arrive the way gravity does in TRPG, as something the
relationships do, and not as a symmetry put in by hand? Every step below is checked by direct computation and then attacked.

 1. Entanglement lives on the spheres of the bag's number systems (two qubits <-> a pair of quaternions; Mosseri & Dandoloff 2001).
 2. Gravity's input (the area law) holds in the full space, inner directions included.
 3. 'Inner gravity' (Kaluza 1921; Klein 1926): motion in an inner direction IS charge. What works and where it fails.
 4. The 3x3 relationship matrix over each allowed number system: its turns, with 'stand at one spot, face one way' (rule 4).
    Octonions: exactly SU(3) x SU(2) x U(1), 12 carriers (Dubois-Violette & Todorov 2018; Krasnov 2019; Baez 2020; Szangolies 2025).
 5. What the weak force acts on, and the charges (hypercharge).
 6. Quantum consistency (anomalies): what the right-handed partners must be.
 7. Attacks: rule 1 at matrix level; the right-handed force of the 4x4 picture; what is still assumed.

Credit: Jacobson 1995; Bombelli, Koul, Lee & Sorkin 1986; Srednicki 1993; Van Raamsdonk 2010; Swingle 2012; Maldacena & Susskind 2013;
Kaluza 1921; Klein 1926; DeWitt 1963; Kerner 1968; Witten 1981, 1983; Nielsen & Ninomiya 1981; Levin & Wen 2005; Wen 2013;
Harlow & Ooguri 2019; Banks & Seiberg 2011; Kirklin 2022/2026; Mosseri & Dandoloff 2001; Bernevig & Chen 2003; Wootters 1998;
Jordan, von Neumann & Wigner 1934; Freudenthal; Chevalley & Schafer 1950 (F4); Gunaydin & Gursey 1973; Dubois-Violette & Todorov 2018;
Todorov 2019-22; Krasnov 2019; Baez 2020-21; Boyle 2020; Szangolies 2025; Furey 2018; Bouchiat, Iliopoulos & Meyer 1972; Gross & Jackiw 1972;
Geng & Marshak 1989; Minahan, Ramond & Warner 1990; Babu & Mohapatra 1989; Foot, Joshi, Lew & Volkas 1990; Georgi 1979; PDG 2024."""
import math, itertools, numpy as np
rng = np.random.default_rng(35)
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
def conj(x): return np.concatenate([x[:1], -x[1:]])
def mul(a, b):
    n = len(a)
    if n == 1: return a*b
    h = n//2; p, q, r, s = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])
def rank_of(alg):
    if not alg: return 0
    best = 10**9
    for _ in range(4):
        X = sum(rng.normal()*g for g in alg)
        best = min(best, len(alg) - int(np.linalg.matrix_rank(np.array([(X@g - g@X).ravel() for g in alg]).T, tol=1e-7)))
    return best
def span_basis(mats, tol=1e-7):
    if not mats: return []
    A = np.array([m.ravel() for m in mats]); U, S, Vt = np.linalg.svd(A, full_matrices=False); r = int((S > tol*max(1, S[0])).sum())
    return [Vt[i].reshape(mats[0].shape) for i in range(r)]
def restrict(alg, cond):
    """Sub-family of alg on which the linear condition cond(D) (a vector) vanishes."""
    A = np.array([cond(g) for g in alg]).T
    if A.size == 0: return list(alg)
    U, S, Vt = np.linalg.svd(A); k = len(alg) - int((S > 1e-7).sum())
    return [sum(c*g for c, g in zip(Vt[-i-1], alg)) for i in range(k)]
def commutant(big, small): return restrict(big, lambda g: np.concatenate([(g@s - s@g).ravel() for s in small]))
def derived(alg): return span_basis([a@b - b@a for a in alg for b in alg])

head("1. Entanglement already lives on the bag's number systems")
errs = []; seps = []
for t in range(2000):
    c = rng.normal(size=4) + 1j*rng.normal(size=4); c /= np.linalg.norm(c)
    a = np.array([c[0].real, c[0].imag, c[1].real, c[1].imag]); b = np.array([c[2].real, c[2].imag, c[3].real, c[3].imag])   # a = c0 + c1 j, b = c2 + c3 j
    hopf = 2*mul(a, conj(b)); conc = 2*abs(c[0]*c[3] - c[1]*c[2]); errs.append(abs(np.hypot(hopf[2], hopf[3]) - conc))
for t in range(500):
    p = rng.normal(size=2) + 1j*rng.normal(size=2); q = rng.normal(size=2) + 1j*rng.normal(size=2); c = np.kron(p, q); c /= np.linalg.norm(c)
    a = np.array([c[0].real, c[0].imag, c[1].real, c[1].imag]); b = np.array([c[2].real, c[2].imag, c[3].real, c[3].imag]); h = 2*mul(a, conj(b)); seps.append(np.hypot(h[2], h[3]))
pr("  Two two-state systems (qubits) = 4 complex amplitudes = a PAIR of quaternions; all states of size 1 form the 7-sphere.")
pr(f"  The quaternion map (a, b) -> 2 a b* has two 'extra' components. For 2,000 random states they equal the entanglement (Wootters' concurrence)")
pr(f"  to within {max(errs):.1e}; for 500 unentangled states they are zero to within {max(seps):.1e}.")
pr("  -> The quaternion structure of a pair IS its entanglement (Mosseri & Dandoloff 2001). Three qubits do the same with a pair of octonions")
pr("     (Bernevig & Chen 2003), and the ladder stops at the octonions for the bag's own reason: the next system has holes.")

head("2. Gravity's input holds in the full space: the area law with an inner ring added")
def ent_block(N, n, L, m=0.05):
    idx = lambda i, j, k: (i*N + j)*n + k; K = np.zeros((N*N*n, N*N*n))
    for i in range(N):
        for j in range(N):
            for k in range(n):
                a = idx(i, j, k); K[a, a] += m*m
                nb = [(i+1, j, k), (i, j+1, k)] + ([(i, j, (k+1) % n)] if n > 2 else [])
                for ii, jj, kk in nb:
                    if ii < N and jj < N:
                        b = idx(ii, jj, kk); K[a, a] += 1; K[b, b] += 1; K[a, b] -= 1; K[b, a] -= 1
    w, V = np.linalg.eigh(K); X = 0.5*(V/np.sqrt(w))@V.T; P = 0.5*(V*np.sqrt(w))@V.T
    o = (N - L)//2; A = [idx(i, j, k) for i in range(o, o+L) for j in range(o, o+L) for k in range(n)]
    nu = np.sqrt(np.clip(np.linalg.eigvals(X[np.ix_(A, A)]@P[np.ix_(A, A)]).real, 0.25, None)); nu = nu[nu > 0.5 + 1e-9]
    return float(np.sum((nu+0.5)*np.log(nu+0.5) - (nu-0.5)*np.log(nu-0.5)))
pr(f"  {'inner ring size':>16s} {'block side L':>13s} {'entropy S':>10s} {'S / (boundary 4L x ring)':>26s}")
for n in (1, 3, 5):
    for L in (4, 6, 8):
        S = ent_block(14, n, L); pr(f"  {n:16d} {L:13d} {S:10.3f} {S/(4*L*n):26.4f}")
pr("  -> With an inner ring at every point the entanglement of a region still grows with its BOUNDARY (last column steady in L for each ring),")
pr("     and each wave that runs round the ring adds its own boundary term (smaller for tighter rings). So the input Jacobson's argument needs")
pr("     survives when the space has inner directions; his argument then applies to the whole space, inner directions included.")

head("3. 'Inner gravity': a particle moving in an inner direction behaves exactly as a charged particle (Kaluza 1921; Klein 1926)")
worst = 0
for t in range(200):
    A = rng.normal(size=4); eta4 = np.diag([-1.0, 1, 1, 1]); g = np.zeros((5, 5)); g[:4, :4] = eta4 + np.outer(A, A); g[:4, 4] = A; g[4, :4] = A; g[4, 4] = 1
    gi = np.linalg.inv(g); p = rng.normal(size=4); q = rng.normal(); P5 = np.concatenate([p, [q]])
    worst = max(worst, abs(P5@gi@P5 - ((p - q*A)@eta4@(p - q*A) + q*q)))
pr(f"  Energy function of a free particle in the 5-dimensional shape = that of a charge q in a field A, where q = momentum round the inner ring.")
pr(f"  Checked for 200 random fields: largest difference {worst:.1e}. Charge conservation = the ring's symmetry; the gauge freedom = relabelling the ring.")
G = 6.70883e-39; e = math.sqrt(4*math.pi/137.036); mKK = e/math.sqrt(16*math.pi*G)
pr("  What this gives for free: the gauge principle (a force for each direction the inner space can be turned without changing sizes).")
pr(f"  Where the PURE version fails (two published results):")
pr(f"    a) a particle whose charge is inner motion must weigh about e/sqrt(16 pi G) = {mKK:.1e} GeV; the electron weighs 5.1e-4 GeV. Off by {mKK/5.11e-4:.0e}.")
pr("    b) handed (left-only) matter cannot come out of pure geometry on a smooth inner space (Witten 1983).")
pr("  -> 'Everything, matter included, is gravity of the inner space' FAILS. What survives is the force side: inner turns that keep sizes are forces.")
pr("     The same conclusion comes from emergent gravity directly: where gravity comes from entanglement, an exact symmetry cannot stay 'global';")
pr("     it must come with a force field (Harlow & Ooguri 2019; proven in the anti-de Sitter setting, expected in general).")

head("4. The turns of the 3x3 relationship matrix, for each number system the bag allows")
def jordan_algebra(n):
    """3x3 hermitian matrices over the n-dimensional number system. Coordinates: d1,d2,d3, x1 (entry 23), x2 (entry 31), x3 (entry 12)."""
    dim = 3 + 3*n; one = np.eye(n)[0]
    def to_mat(v):
        d = v[:3]; x1, x2, x3 = v[3:3+n], v[3+n:3+2*n], v[3+2*n:]
        return [[d[0]*one, x3, conj(x2)], [conj(x3), d[1]*one, x1], [x2, conj(x1), d[2]*one]]
    def from_mat(M): return np.concatenate([[M[0][0][0], M[1][1][0], M[2][2][0]], M[1][2], M[2][0], M[0][1]])
    def mm(X, Y): return [[sum(mul(X[i][j], Y[j][k]) for j in range(3)) for k in range(3)] for i in range(3)]
    def jp(u, v):
        X, Y = to_mat(u), to_mat(v); A, B = mm(X, Y), mm(Y, X); return from_mat([[(A[i][k] + B[i][k])/2 for k in range(3)] for i in range(3)])
    E = np.eye(dim); T = np.zeros((dim, dim, dim))
    for i in range(dim):
        for j in range(i, dim): T[i, j] = T[j, i] = jp(E[i], E[j])
    return dim, T
def derivations(dim, T):
    rows = []
    for i in range(dim):
        for j in range(i, dim):
            for c in range(dim):
                row = np.zeros((dim, dim)); row[c, :] += T[i, j, :]; row[:, i] -= T[:, j, c]; row[:, j] -= T[i, :, c]; rows.append(row.ravel())
    M = np.array(rows); G = M.T@M; w, V = np.linalg.eigh(G); k = int((w < 1e-8*w[-1]).sum())
    return [V[:, i].reshape(dim, dim) for i in range(k)]
names = {2: "complex", 4: "quaternions", 8: "octonions"}; store = {}
pr(f"  {'number system':13s} {'all turns':>14s} {'+ face one way':>16s} {'+ stand at one spot':>20s} {'both (rule 4)':>16s}")
for n in (2, 4, 8):
    dim, T = jordan_algebra(n); der = derivations(dim, T)
    W = [0, 1, 2] + [3 + b*n + c for b in range(3) for c in (0, 1)]        # entries built from 1 and one chosen direction u: an ordinary complex 3x3 matrix
    notW = [i for i in range(dim) if i not in W]
    face = restrict(der, lambda g: g[np.ix_(notW, W)].ravel())
    spot = restrict(der, lambda g: g[:, 2])                                  # the state 'slot 3' stays put
    both = restrict(face, lambda g: g[:, 2])
    f = lambda a: f"{len(a)} (rank {rank_of(a)})"
    pr(f"  {names[n]:13s} {f(der):>14s} {f(face):>16s} {f(spot):>20s} {f(both):>16s}")
    store[n] = (dim, T, der, face, spot, both, W)
pr("  Measured: 12 carriers, 4 of them neutral, in blocks 8 | 3 | 1.")
dim, T, der, face, spot, both, W = store[8]
dalg = derived(both); cen = commutant(both, both)
colour = restrict(both, lambda g: g[:, W].ravel())               # turns that leave every complex (colourless) entry alone
weak4 = commutant(both, colour); weak = derived(weak4)
pr(f"\n  Octonion case in detail: {len(both)} carriers, {rank_of(both)} neutral. Self-interacting part {len(dalg)}, plus {len(cen)} that commutes with everything.")
pr(f"    turns that leave the colourless entries alone : {len(colour)} carriers, {rank_of(colour)} neutral   -> SU(3), colour")
pr(f"    turns that commute with colour               : {len(weak4)} carriers, {rank_of(weak4)} neutral   -> {len(weak)} self-interacting (SU(2), the weak force) + {len(weak4)-len(weak)} single (U(1), hypercharge)")
ok = max(np.abs(a@b - b@a).max() for a in colour for b in weak4)
pr(f"    the two blocks do not mix (largest cross-commutator {ok:.1e}); together {len(span_basis(colour + weak4))} = all {len(both)}.")
pr("  -> ONE meaning of symmetry (turns that keep how relationship matrices combine: the same meaning an earlier test used for colour), on ONE object")
pr("     (the 3x3 octonion relationship matrix of an earlier test), with rule 4 read as 'stand at one spot AND face one way', gives exactly 8 | 3 | 1.")
pr("     This is the published result of Dubois-Violette & Todorov (2018) and Krasnov (2019), reproduced here from scratch.")
pr("  Note the complex row: the weak block (4 carriers, 2 neutral) is already there for an ordinary complex 3x3 relationship matrix, the kind")
pr("  rule 5 gives the leptons. The weak force comes from the MATRIX (the relationships); colour comes from the octonion's extra directions.")

head("5. What the weak force acts on, and with what charges")
n = 8; links = list(range(3, 3 + 2*n))                    # x1, x2: the links from the chosen spot to the other two slots
Yg = cen[0]; Yl = Yg[np.ix_(links, links)]; ev = np.linalg.eigvals(Yl); ys = np.round(np.abs(ev.imag)/np.abs(ev.imag).min(), 6)
vals, counts = np.unique(ys, return_counts=True)
pr(f"  The 16 link components (8 complex states) carry the single carrier's charge in the ratio {vals.tolist()} with multiplicities {counts.tolist()} (real components).")
# normalise su(2) and compute its Casimir on the links
B = span_basis([w_[np.ix_(links, links)] for w_ in weak]); B = [b/np.sqrt(-np.trace(b@b)) for b in B]
c12 = B[0]@B[1] - B[1]@B[0]; s = np.sqrt(-np.trace(c12@c12)); Tn = [b/s for b in B]
Cw = -sum(t@t for t in Tn); cw = np.unique(np.round(np.linalg.eigvalsh((Cw + Cw.T)/2), 6))
Cc_basis = span_basis([c_[np.ix_(links, links)] for c_ in colour]); Cc = -sum(c_@c_ for c_ in Cc_basis); cc = np.linalg.eigvalsh((Cc + Cc.T)/2)
nsing = int((np.abs(cc) < 1e-8).sum())
pr(f"  Weak 'spin' of every link state: j(j+1) = {cw.tolist()} -> all are doublets (j = 1/2). Colour: {16-nsing} components feel it (3 colours x 2 x 2), {nsing} do not.")
vecs = np.linalg.eigh((Cc + Cc.T)/2)[1]; P0 = vecs[:, np.abs(cc) < 1e-8]; P3 = vecs[:, np.abs(cc) > 1e-8]
ysz = lambda P: math.sqrt(abs(np.trace(P.T@Yl@Yl@P))/P.shape[1])
pr(f"  The colourless doublet's charge is {ysz(P0)/ysz(P3):.3f} times the coloured doublet's in size.")
pr("  -> One colourless weak doublet and one coloured weak doublet, with charges 3 : 1 in size. That is the LEFT-HANDED lepton pair (neutrino,")
pr("     electron; hypercharge 1/2 in size) and quark pair (up, down; 1/6). The relative sign is only a naming (a weak pair and its antiparticle pair")
pr("     are the same object here); with the usual naming it is -1/2 and +1/6, and 3 x (1/6) = 1/2 is exactly what step 6 needs.")
pr("     Only doublets appear: the weak force acts on every state in this space.")

head("6. Quantum consistency: what must the partners of those doublets be?")
import sympy as sp
yQ, yL = sp.Rational(1, 6), sp.Rational(-1, 2)
pr(f"  Doublets alone: weak-weak-charge = {3*yQ + yL} (fine); colour-colour-charge = {2*yQ}; charge cubed = {6*yQ**3 + 2*yL**3}; so the doublets alone are NOT consistent.")
yu, yd, ye = sp.symbols("y_u y_d y_e")
sol = sp.solve([2*yQ - yu - yd, 6*yQ + 2*yL - 3*yu - 3*yd - ye, 6*yQ**3 + 2*yL**3 - 3*yu**3 - 3*yd**3 - ye**3], [yu, yd, ye], dict=True)
pr(f"  Adding weak-singlet partners (two coloured, one colourless) with unknown charges, the three consistency conditions have the solutions: {sol}")
pr("  -> Unique (up to which coloured one is called 'up'): +2/3, -1/3 and -1. Those are the measured right-handed up, down and electron.")
pr("     A neutral partner for the neutrino (TRPG's hidden partner) is allowed and, with its own mass, fixes the charges completely (Babu & Mohapatra 1989).")
bl15 = 6*sp.Rational(1, 3) - 3*sp.Rational(1, 3) - 3*sp.Rational(1, 3) - 2 + 1
pr(f"  Emergent gravity forbids exactly conserved 'global' numbers. Baryon-minus-lepton number: sum over 15 states = {bl15} (cannot be a force); with the")
pr(f"  hidden partner (16 states) = {bl15 + 1} (can), and its own mass then breaks it. Either way the hidden partner is the consistent choice. TRPG has it.")

head("7. Attacks")
pr("7a. Does the 3x3 octonion matrix obey rule 1 (no holes)?")
worst = 1e9
for t in range(2000):
    v = rng.normal(size=dim); sq = np.einsum("i,j,ijk->k", v, v, T); tr = sq[:3].sum(); worst = min(worst, tr/(v@v))
ii = np.zeros(8, dtype=complex); ii[1] = 1; ii[2] = 1j
pr(f"    Nothing squares to nothing: trace of X o X is at least {worst:.2f} x |X|^2 over 2,000 random matrices (it is the sum of squares of all entries).")
pr(f"    By contrast complex x octonions has (e1 + i e2)^2 = {np.abs(mul(ii, ii)).max():.1f}: something that squares to nothing.")
pr("    Two DIFFERENT slots do multiply to zero (they do not overlap), as in every matrix TRPG already uses. Rule 1, as 'no relationship squares to")
pr("    nothing and squares never cancel', holds (Jordan, von Neumann & Wigner 1934). The entries are single octonions: no joined number systems.")
pr("7b. The right-handed force of the 4x4 picture")
pr(f"    All turns of the 3x3 octonion matrix: {len(der)} carriers, only {rank_of(der)} neutral. The 4x4 picture with a right-handed force needs 5 neutral.")
pr("    -> It cannot fit. This route and the right-handed force (W_R, the 'right side's scale') exclude each other. TRPG has to choose.")
pr("       Choosing this route removes the earlier break 2 entirely: no 4x4 symmetry means no forced 'down = electron' mass relations, no lopsided")
pr("       content, no proton-decay tension, no monopoles. It also removes the right side's carriers that made the hidden partners (see the matter-excess note).")
pr("7c. Still assumed, stated plainly")
pr("    i)   'Stand at one spot' is an added reading of rule 4. Without it the count is 16 (rank 4), not 12.")
pr("    ii)  The three slots here are 'the spot + the two members of a weak pair'. In an earlier test the three placements are the three generations.")
pr("         Whether they can be both is NOT shown (Boyle 2020 proposes that the three choices of spot are the three generations).")
pr("    iii) Why the right-handed partners are weak singlets, and not mirror doublets, is argued (paired states freeze out heavy; Georgi 1979), not derived.")
pr("    iv)  That a symmetry must be a force: a theorem of relativity + quantum mechanics, and of emergent gravity where that is under control.")
pr("    v)   The remaining 40 turns: no carriers if the spot and facing are facts of observation (as in an earlier test); 40 heavy carriers if they are dynamical.")

pr("7d. Does it matter which spot and which facing direction?")
from scipy.linalg import expm
E33 = np.zeros(dim); E33[2] = 1
cnt_al, cnt_mis = [], []
for t in range(4):
    g = expm(sum(rng.normal()*d_ for d_ in face)); sp_al = g@E33               # another spot inside the plane you face
    cnt_al.append(len(restrict(face, lambda h: h@sp_al)))
    g = expm(sum(rng.normal()*d_ for d_ in der)); sp_mis = g@E33               # a spot outside that plane
    cnt_mis.append(len(restrict(face, lambda h: h@sp_mis)))
pr(f"    4 other spots inside the plane being faced: carriers {cnt_al} (same forces every time, as with colour in an earlier test).")
pr(f"    4 spots outside that plane: carriers {cnt_mis}. The spot must lie in the view: 'you stand inside what you are looking at'.")
pr("7e. What else is in the matrix besides the weak pairs? (the entries not linked to the spot)")
Yfull = cen[0]; ycol = ysz(P3); Yn = Yfull/ycol/6                               # scale: coloured doublet = 1/6
ipl = lambda a, b: -np.trace(a[np.ix_(links, links)]@b[np.ix_(links, links)])     # orthonormalise the weak carriers on the link states
wn = []
for w_ in span_basis(weak):
    for v_ in wn: w_ = w_ - ipl(w_, v_)*v_
    wn.append(w_/math.sqrt(ipl(w_, w_)))
wn = [w_/s for w_ in wn]
rest = [0, 1] + list(range(3 + 2*n, 3 + 3*n))
Yr = Yn[np.ix_(rest, rest)]; yv = np.round(np.sqrt(np.abs(np.linalg.eigvalsh(Yr@Yr))), 4); vals, counts = np.unique(yv, return_counts=True)
Cr = -sum(w_[np.ix_(rest, rest)]@w_[np.ix_(rest, rest)] for w_ in wn); cr = np.abs(np.round(np.linalg.eigvalsh((Cr + Cr.T)/2), 3)); v2, c2 = np.unique(cr, return_counts=True)
pr(f"    10 components: charge sizes {vals.tolist()} with multiplicities {counts.tolist()}; weak spin j(j+1) = {v2.tolist()} with multiplicities {c2.tolist()}.")
pr("    -> a neutral weak triplet (3), a coloured weak singlet of charge 1/3 in size (6: the size of the right-handed DOWN quark's charge), and 1 neutral singlet.")
pr("       Suggestive only; not used for anything here.")
pr("7f. A sharper version that could fail: if the whole matrix had ONE strength, the three forces would be tied together.")
T3 = wn[0]; k = np.trace(Yn@Yn)/np.trace(T3@T3); s2 = 1/(1 + k)
tm = (29.59 - 1/0.1180)/(7 - 19/6); mu = 91.1876*math.exp(2*math.pi*tm)
ia2 = 29.59 + (19/6)*tm; iaY = (5/3)*(59.00 - 4.1*tm); s2_need = ia2/(ia2 + iaY)
pr(f"    One strength for all 12 carriers predicts the weak mixing number sin^2(theta_W) = {s2:.3f} at the energy where the strong and weak strengths meet.")
pr(f"    Running the measured strengths up (one loop), they meet at {mu:.1e} GeV, where the measured-and-run value is {s2_need:.3f}.")
pr(f"    -> {'FAILS' if abs(s2 - s2_need) > 0.03 else 'passes'}: the route gives the KINDS of force and their charges, not their strengths. The three strengths stay measured inputs.")
pr("7g. Can the three slots ALSO be the three generations, with each generation 'standing at' its own slot? (Boyle 2020's proposal, plain form)")
def closure(alg):
    cur = span_basis(alg)
    while True:
        nxt = span_basis(cur + [a@b - b@a for a in cur for b in cur])
        if len(nxt) == len(cur): return nxt
        cur = nxt
weaks = []
for k in range(3):
    bk = restrict(face, lambda g: g[:, k]); wk4 = commutant(bk, colour); weaks.append(derived(wk4))
    pr(f"    standing at slot {k+1}: {len(bk)} carriers, of which colour {len(restrict(bk, lambda g: g[:, W].ravel()))} (the same 8 every time) and weak block {len(weaks[-1])} + {len(wk4) - len(weaks[-1])}")
same = len(span_basis(weaks[0] + weaks[1])); allw = closure(weaks[0] + weaks[1] + weaks[2])
mix = max(np.abs(a@b - b@a).max() for a in weaks[0] for b in weaks[1])
pr(f"    The weak blocks of slot 1 and slot 2 together span {same} directions (3 would mean 'the same weak force'); they do not commute (size {mix:.2f}).")
pr(f"    All three together close into {len(allw)} carriers (rank {rank_of(allw)}): a force that turns the generations into each other.")
pr(f"    Counting: the matrix has {3*n} link components in all; three generations of weak pairs need 3 x 16 = 48.")
pr("    Measured: ONE weak force, the same for all three generations (tested to about 0.1-0.2% in tau, muon and pion decays), and 3 weak carriers, not 8.")
pr("    -> FAILS in the plain form: the same 3x3 matrix cannot give both the weak force (slots = spot + pair) and the three generations (slots = generations).")
pr("       The published way round this doubles the matrix with complex numbers (Boyle 2020), and complex x octonions has holes (7a): not allowed by rule 1.")
pr("       So if the matrix gives the forces, the earlier 'at most 3 generations because octonion matrices stop at 3x3' is no longer available; and the reverse.")
