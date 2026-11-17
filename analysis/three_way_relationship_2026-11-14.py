"""The three-way idea: two relationships and not two things; a third gives the exact octonion's
worth of states; two of the same merge and seem to be one. Test: put a THREE-way product on the link states of the relationship matrix
(seen from the observer's spot) and ask (1) are there holes, (2) which turns keep two chosen states, (3) what the 16 link states then are.
Uses the tools of forces_from_entanglement_2026-11-14.py (same folder).
Credit: Jordan 1933; Jordan, von Neumann & Wigner 1934; Loos 1977 (Jordan triples, the 16-dimensional exceptional one); Koecher; Meyberg; Cartan 1935;
Chevalley 1954 (pure spinors); Georgi & Glashow 1974; Georgi 1975; Fritzsch & Minkowski 1975; Baez & Huerta 2010; Baez 2026; Baez & Schwahn 2026;
Dubois-Violette & Todorov 2018-19; Krasnov 2019; Furey 2015; Harari 1979; Shupe 1979; Pauli 1925; Adler; Bell & Jackiw; PDG 2024."""
import math, itertools, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "forces_from_entanglement_2026-11-14.py")).read()
exec(src[src.index("import math, itertools"):src.index('head("1. Entanglement')])
exec(src[src.index("def jordan_algebra(n):"):src.index("names = {2:")])
def restrict(alg, cond):                                   # memory-light version of the same tool
    A = np.array([cond(g) for g in alg]).T
    if A.size == 0: return list(alg)
    w, V = np.linalg.eigh(A.T@A); k = int((w < 1e-10*max(1.0, w[-1])).sum())
    return [sum(c*g for c, g in zip(V[:, i], alg)) for i in range(k)]
def commutant(big, small): return restrict(big, lambda g: np.concatenate([(g@s_ - s_@g).ravel() for s_ in small]))
n = 8; dim, T = jordan_algebra(n); Tc = T.astype(complex); lk = list(range(3, 19))
def jp(a, b): return np.einsum('i,j,ijk->k', a, b, Tc)
def emb(v):
    z = np.zeros(dim, complex); z[lk] = v; return z
def tri27(x, y, z):
    yc = np.conj(y); return jp(jp(x, yc), z) + jp(jp(z, yc), x) - jp(jp(x, z), yc)
E16 = np.eye(16, dtype=complex); TT = np.zeros((16, 16, 16, 16), complex); leak = 0.0
for i in range(16):
    for j in range(16):
        for k in range(16):
            r = tri27(emb(E16[i]), emb(E16[j]), emb(E16[k])); TT[i, j, k] = r[lk]; leak = max(leak, np.abs(np.delete(r, lk)).max())
def t(x, y, z): return np.einsum('a,b,c,abck->k', x, np.conj(y), z, TT)
def L(x, y): return np.einsum('a,b,abck->kc', x, np.conj(y), TT)
def R(M): return np.block([[M.real, -M.imag], [M.imag, M.real]])
def C(M): return M[:16, :16] + 1j*M[16:, :16]
rv = lambda v: np.concatenate([v.real, v.imag]); herm = lambda M: (M + M.conj().T)/2
rc = lambda: rng.normal(size=16) + 1j*rng.normal(size=16)

head("1. The three-way product on the link states")
pr("  The link states seen from the observer's spot: two octonions (16 numbers). Written with an outside i they are 16 complex states.")
pr("  Two-way product: x times y. Three-way product: {x, y, z} = (x.y*).z + (z.y*).x - (x.z).y*, where y* is y's mirror image (Jordan 1933; Loos 1977).")
a, b, x, y, z = rc(), rc(), rc(), rc(), rc()
res = np.abs(t(a, b, t(x, y, z)) - (t(t(a, b, x), y, z) - t(x, t(b, a, y), z) + t(x, y, t(a, b, z)))).max()
pr(f"  It never leaves the link states (leak {leak:.0e}) and obeys the standard consistency law for three-way products (error {res:.0e}).")
O8 = np.array([[mul(np.eye(8)[i], np.eye(8)[j]) for j in range(8)] for i in range(8)]).astype(complex); v8 = np.array([0, 1, 1j, 0, 0, 0, 0, 0])
sq = np.abs(np.einsum('i,j,ijk->k', v8, v8, O8)).max(); a_ = np.zeros(16, complex); a_[1] = 1; a_[2] = 1j; c3 = t(a_, a_, a_)
pr(f"  Two-way: the entry a = e1 + i e2 squares to nothing (size {sq:.0e}). That is the hole that ruled out joined number systems in earlier tests")
pr(f"  Three-way: the same entry gives {{a, a, a}} = {(c3[1]/a_[1]).real:.0f} a. Not nothing.")
H = np.array([[np.trace(L(E16[i], E16[j])) for j in range(16)] for i in range(16)]); ev = np.linalg.eigvalsh(herm(H))
mn = min(np.linalg.norm(t(x, x, x))/np.linalg.norm(x)**3 for x in (rc() for _ in range(20000)))
pr(f"  For every state: the weight of 'x with its own mirror' is positive (all 16 weights = {ev.min():.0f}); smallest |{{x,x,x}}| / |x|^3 in 20,000 random states: {mn:.2f}.")
pr("  -> With a three-way product NOTHING combines with itself into nothing. The doubled link space has no three-way holes (Loos 1977: it is 'positive').")

head("2. The turns of the three-way product")
basis = [E16[i] for i in range(16)] + [1j*E16[i] for i in range(16)]
K = span_basis([R(L(basis[i], basis[j]) - L(basis[j], basis[i])) for i in range(32) for j in range(i + 1, 32)], tol=1e-9); Kd = derived(K)
der = derivations(dim, T); W = [0, 1, 2] + [3 + b*n + c for b in range(3) for c in (0, 1)]; notW = [i for i in range(dim) if i not in W]
face = restrict(der, lambda g: g[np.ix_(notW, W)].ravel()); spot = restrict(der, lambda g: g[:, 2]); both12 = restrict(face, lambda g: g[:, 2])
V = lambda g: g[np.ix_(lk, lk)]; old12 = [R(V(g).astype(complex)) for g in both12]; old36 = [R(V(g).astype(complex)) for g in spot]
pr(f"  Turns that keep the three-way product: {len(K)} ({len(Kd)} self-interacting + {len(commutant(K, K))} single: the overall phase of every state at once), {rank_of(K)} neutral.")
pr(f"  The matrix's own 36 turns that keep the observer's spot are among them: {len(span_basis(K + old36, tol=1e-9)) == len(K)}. So the three-way product has MORE symmetry (46) than the matrix seen from the spot (36).")

head("3. The observer's spot and facing direction pick two states")
def st(slot, sign):
    v = np.zeros(16, complex); v[8*slot] = 1; v[8*slot + 1] = sign*1j
    return v/math.sqrt((np.vdot(v, t(v, v, v))/np.vdot(v, v)).real)
cands = {(s, g): st(s, g) for s in (0, 1) for g in (+1, -1)}
pr("  The colourless link states are built from 1 and u (the facing direction): link 1 x (1 +/- i u) and link 2 x (1 +/- i u). Each is a simplest possible state:")
for k_, v in cands.items():
    w = np.round(np.linalg.eigvals(L(v, v)).real, 4); pr(f"    link {k_[0] + 1} (1 {'+' if k_[1] > 0 else '-'} i u): {{e,e,e}} = e (error {np.abs(t(v, v, v) - v).max():.0e}); how the other states relate to it: {dict(zip(*[x_.tolist() for x_ in np.unique(w, return_counts=True)]))}")
exact = lambda v: (lambda g: g@rv(v))
def ray(v):
    def f(g):
        gv = C(g)@v; return rv(gv - v*np.vdot(v, gv)/np.vdot(v, v))
    return f
two = lambda c1, c2: (lambda g: np.concatenate([c1(g), c2(g)]))
def info(A): return f"{len(A):>2} turns = {len(derived(A)) if A else 0:>2} self-interacting + {len(commutant(A, A)) if A else 0} single; {rank_of(A)} neutral"
pr("  Two such states are either 'unrelated' (0) or 'half-related' (1/2):")
for p, q in itertools.combinations(cands, 2):
    lam = (np.vdot(cands[q], t(cands[p], cands[p], cands[q]))/np.vdot(cands[q], cands[q])).real
    pr(f"    link {p[0]+1} ({'+' if p[1]>0 else '-'}) with link {q[0]+1} ({'+' if q[1]>0 else '-'}): {lam:.1f}   turns of the 45 keeping the first exactly and the second's direction: {info(restrict(Kd, two(exact(cands[p]), ray(cands[q]))))}")
e1, e2 = cands[(0, 1)], cands[(1, -1)]

head("4. Which turns keep the two chosen states? Every reading")
for nm, alg in (("all 46 turns", K), ("the 45 self-interacting turns", Kd)):
    for n1, c1 in (("exactly", exact(e1)), ("as a direction", ray(e1))):
        for n2, c2 in (("exactly", exact(e2)), ("as a direction", ray(e2))):
            pr(f"    {nm:<30} state 1 {n1:<15} state 2 {n2:<15}: {info(restrict(alg, two(c1, c2)))}")
pr("  TRPG's rule 4 says 'stand at one spot and face one way'. A spot is kept exactly; a facing is a direction. Leaving out the overall phase (which no")
pr("  measurement can see), that reading gives 12 turns with 4 neutral: 11 self-interacting and 1 single.")
G = restrict(Kd, two(exact(e1), ray(e2))); D = derived(G); Z = commutant(G, G)
B = []
for g in D:
    w_ = g.copy()
    for v_ in B: w_ = w_ - np.trace(w_@v_.T)/np.trace(v_@v_.T)*v_
    B.append(w_/math.sqrt(np.trace(w_@w_.T)))
ad = lambda X: np.array([[np.trace((X@B[j] - B[j]@X)@B[i].T) for j in range(len(B))] for i in range(len(B))])
Cas = -sum(ad(b_)@ad(b_) for b_ in B); w, v = np.linalg.eigh((Cas + Cas.T)/2); vals = np.unique(np.round(w, 4)); grp = {}
for val in vals: grp[int((np.abs(w - val) < 1e-3).sum())] = [sum(v[i, k]*B[i] for i in range(len(B))) for k in range(len(w)) if abs(w[k] - val) < 1e-3]
pr(f"  The self-interacting 11 split into blocks of {sorted(grp, reverse=True)}: so the 12 are 8 | 3 | 1.")

head("5. What the 16 link states are under these 12")
col, wk = grp[8], grp[3]; Cc = -sum(C(c)@C(c) for c in col); Cw = -sum(C(c)@C(c) for c in wk); Y = -1j*C(Z[0])
cmax = np.linalg.eigvalsh(herm(Cc)).max(); wmax = np.linalg.eigvalsh(herm(Cw)).max(); Cc = Cc*(4/3)/cmax; Cw = Cw*(3/4)/wmax
ev, U = np.linalg.eigh(herm(Cw) + math.pi*herm(Cc) + math.e*herm(Y)); raw = []
for k in range(16):
    u_ = U[:, k]; raw.append(((u_.conj()@Cc@u_).real, (u_.conj()@Cw@u_).real, (u_.conj()@Y@u_).real))
yq = [r[2] for r in raw if r[0] > 1 and r[1] > 0.5][0]; rows = [(a, b, c/yq/6) for a, b, c in raw]; tab = {}
for a, b, c in rows: key = (round(a, 2), round(b, 2), round(c, 4)); tab[key] = tab.get(key, 0) + 1
nm = {(1.33, 0.75, 0.1667): "quark pair", (0.0, 0.75, -0.5): "lepton pair", (1.33, 0.0, 0.3333): "down partner", (1.33, 0.0, -0.6667): "UP PARTNER",
      (0.0, 0.0, 1.0): "ELECTRON PARTNER", (0.0, 0.0, 0.0): "hidden neutrino partner"}
pr(f"    {'colours':>8} {'weak':>8} {'charge':>8} {'states':>7}   particle")
for key, m in sorted(tab.items(), key=lambda q: -q[1]): pr(f"    {'3' if key[0] > 1 else '1':>8} {'pair' if key[1] > 0.5 else 'single':>8} {key[2]:>+8.4f} {m:>7}   {nm.get((abs(key[0]), abs(key[1]), key[2] + 0.0), '?')}")
ys = np.array([r[2] for r in rows]); cs = np.array([r[0] for r in rows]); ws = np.array([r[1] for r in rows])
pr(f"  Consistency sums (must all vanish): charges {ys.sum():+.1e}; cubes {(ys**3).sum():+.1e}; on coloured states {ys[cs > 1].sum():+.1e}; on weak pairs {ys[ws > 0.5].sum():+.1e}.")
q1 = (e1.conj()@Y@e1).real/np.vdot(e1, e1).real/yq/6; q2 = (e2.conj()@Y@e2).real/np.vdot(e2, e2).real/yq/6
pr(f"  The two chosen states themselves: state 1 has charge {q1:+.3f} (the hidden neutrino partner); state 2 has charge {q2:+.3f} (the electron partner).")
pr("  -> A WHOLE generation: both weak pairs and all four partners, with the measured charges. The up partner and the electron partner are there.")
wn = [c*math.sqrt((3/4)/wmax) for c in wk]; cn = [c*math.sqrt((4/3)/cmax) for c in col]; Yn = C(Z[0])/yq/6
kT = -np.trace(C(wn[0])@C(wn[0])).real; kC = -np.trace(C(cn[0])@C(cn[0])).real; kY = -np.trace(Yn@Yn).real
pr(f"  Strengths from this geometry (one strength for all 12): sin^2 = {1/(1 + kY/kT):.3f}; strong/weak = {kT/kC:.3f}. That is the whole-generation count of an earlier test:")
pr("     the near miss (10.6% off at the Z). Still a strict FAIL; nothing changes there.")

head("6. Attacks")
sh12 = len(G) + 12 - len(span_basis(G + old12, tol=1e-9)); sh36 = len(G) + 36 - len(span_basis(G + old36, tol=1e-9))
pr(f"6a. Is this the weak force of an earlier test? No. These 12 share {sh12} turns with the earlier 12 (the 8 colour turns and 1 more), and only {sh36} of them keep the")
pr("    matrix's own two-way product. The weak force here is a symmetry of the three-way product, not of the matrix. the earlier weak force made every link")
pr("    state a pair; this one makes half of them pairs and half singles, as measured.")
GA = restrict(K, two(exact(e1), exact(e2))); ZA = commutant(GA, GA); evA = np.linalg.eigvalsh(herm(-1j*C(ZA[0]))); evA = evA/np.abs(evA[np.abs(evA) > 1e-8]).min()
if evA[np.abs(evA) > 0.5].mean() < 0: evA = -evA
cntA = dict(zip(*[x_.tolist() for x_ in np.unique(np.round(evA, 3) + 0.0, return_counts=True)]))
pr(f"6b. Other readings. Keeping BOTH states exactly, with the overall phase allowed, also gives 12 turns, but the single carrier's charges are then {cntA}:")
pr("    it counts how many colour switches are on (0, 1, 2 or 3). That is the 'sets of three' and it is a real pattern (Georgi & Glashow 1974; Furey 2015),")
SA = np.array(sorted(np.round(evA).astype(int)))
pr(f"    but as a force it is inconsistent: its charges do not sum to zero ({int(SA.sum())}). Keeping both only as directions gives 13 or 14 turns: an extra carrier that")
pr("    is not seen. Only 'one exactly, one as a direction' gives the measured 12. It is the natural reading of rule 4, and it is still a choice.")
e1b, e2b = cands[(0, -1)], cands[(1, 1)]; Gb = restrict(Kd, two(exact(e1b), ray(e2b))); Zb = commutant(Gb, Gb); Yb = -1j*C(Zb[0])
qv = np.array([(u_.conj()@Yb@u_).real for u_ in U.T]); ratio = qv/np.array([r[2] for r in raw]); 
same = np.allclose(np.sort(np.abs(np.linalg.eigvalsh(herm(Yb)))), np.sort(np.abs(np.linalg.eigvalsh(herm(Y)))), atol=1e-6)
inter = len(G) + len(Gb) - len(span_basis(G + Gb, tol=1e-9))
pr(f"6c. Facing the other way (u -> -u) picks the other half-related pair of states. Its 12 turns share {inter} with the first 12, and its charges have the same sizes: {same}.")
orth = restrict(Kd, two(exact(e1), ray(cands[(1, 1)])))
pr(f"6d. Two UNRELATED states instead of half-related ones: {info(orth)}. That is a four-colour force (lepton as fourth colour), not the measured forces.")
pr("    So the two states must be 'half-related'. Half-related means: each is a relationship partner of the other, neither is the other's mirror.")
pr("6e. The two-way hole is still there in the entries (e1 + i e2 squares to nothing). This route needs rule 1 to be read for three-way combinations.")
pr("    Whether the joined number systems rejected in an earlier test would also pass a three-way reading has NOT been tested.")
pr("6f. Nothing extra: exactly 16 states and 12 carriers. No mirror states: the 16 are not their own mirror image, so matter is one-handed here without")
pr("    the hole argument of an earlier test (provided spin uses this same i, as in an earlier test).")
pr("6g. Not new mathematics: the 16 are the 'spinor' of Georgi 1975 and Fritzsch & Minkowski 1975, without their extra force carriers; the three-way")
pr("    form is Loos 1977 and Baez 2026. What is TRPG's is the reading: spot exactly, facing as a direction, and no three-way holes.")

head("7. The same reading on the three number systems (one meaning throughout)")
pr("  Same construction every time: link states seen from the spot, three-way product, keep one colourless state exactly and the half-related one as a direction.")
for n_ in (2, 4, 8):
    dim_, T_ = jordan_algebra(n_); Tq = T_.astype(complex); m = 2*n_; lq = list(range(3, 3 + m))
    jq = lambda a, b: np.einsum('i,j,ijk->k', a, b, Tq)
    def eq(v):
        z = np.zeros(dim_, complex); z[lq] = v; return z
    def t27(x, y, z):
        yc = np.conj(y); return jq(jq(x, yc), z) + jq(jq(z, yc), x) - jq(jq(x, z), yc)
    Em = np.eye(m, dtype=complex); Tm = np.zeros((m, m, m, m), complex)
    for i in range(m):
        for j in range(m):
            for k in range(m): Tm[i, j, k] = t27(eq(Em[i]), eq(Em[j]), eq(Em[k]))[lq]
    tq = lambda x, y, z: np.einsum('a,b,c,abck->k', x, np.conj(y), z, Tm); Lq = lambda x, y: np.einsum('a,b,abck->kc', x, np.conj(y), Tm)
    Rq = lambda M: np.block([[M.real, -M.imag], [M.imag, M.real]]); Cq = lambda M: M[:m, :m] + 1j*M[m:, :m]
    bq = [Em[i] for i in range(m)] + [1j*Em[i] for i in range(m)]
    Kq = span_basis([Rq(Lq(bq[i], bq[j]) - Lq(bq[j], bq[i])) for i in range(2*m) for j in range(i + 1, 2*m)], tol=1e-9); Kqd = derived(Kq)
    def sq_(slot, sign):
        v = np.zeros(m, complex); v[n_*slot] = 1; v[n_*slot + 1] = sign*1j
        return v/math.sqrt((np.vdot(v, tq(v, v, v))/np.vdot(v, v)).real)
    f1, f2 = sq_(0, 1), sq_(1, -1)
    def cond(g):
        gv = Cq(g)@f2; return np.concatenate([g@rv(f1), rv(gv - f2*np.vdot(f2, gv)/np.vdot(f2, f2))])
    Gq = restrict(Kqd, cond); Dq = derived(Gq) if Gq else []; Zq = commutant(Gq, Gq) if Gq else []
    Hq = np.array([[np.trace(Lq(Em[i], Em[j])) for j in range(m)] for i in range(m)]); pos = np.linalg.eigvalsh(herm(Hq)).min() > 1e-9
    blocks = []
    if Dq:
        Bq = []
        for g in Dq:
            w_ = g.copy()
            for v_ in Bq: w_ = w_ - np.trace(w_@v_.T)/np.trace(v_@v_.T)*v_
            Bq.append(w_/math.sqrt(np.trace(w_@w_.T)))
        adq = lambda X: np.array([[np.trace((X@Bq[j] - Bq[j]@X)@Bq[i].T) for j in range(len(Bq))] for i in range(len(Bq))])
        wq = np.linalg.eigvalsh(herm(-sum(adq(b_)@adq(b_) for b_ in Bq))); blocks = sorted(np.unique(np.round(wq, 3), return_counts=True)[1].tolist(), reverse=True)
    nm_ = {2: "complex numbers", 4: "quaternions", 8: "octonions"}[n_]
    pr(f"    {nm_:<16}: {m:>2} link states; no three-way holes: {pos}; turns of the three-way product {len(Kq):>2}; kept by the perspective: {len(Gq):>2} = blocks {blocks} + {len(Zq)} single; {rank_of(Gq) if Gq else 0} neutral")
pr("  -> One meaning for all three. Only the octonions give 8 | 3 | 1, and the octonions are what the two-way bag rules force for the entries.")
pr("     So the three-way reading does not replace the two-way rule. The two-way rule fixes the numbers in the matrix; the three-way rule is for the states.")
