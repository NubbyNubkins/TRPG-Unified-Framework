"""Falsification pass (14 Nov 2026). Standard: a result that can be proven wrong in this system fails.
Each open item is attacked three ways: (a) does it contradict TRPG's own starting rules; (b) do the theorems it leans on actually apply;
(c) does it clash with measured data when every number is taken at its stated precision.
earlier tests (this file): force groups, quark masses (colour rule, charm, bottom), the 2/9 twist.
an earlier test (matter excess) is in falsification_matter_excess_2026-11-14.py because it needs the ULYSSES code.
Data: PDG 2024 (masses, |Vus|); FLAG 2024 (quark-mass ratios); Huang & Zhou 2021 (running masses at M_Z); Super-Kamiokande 2020 (proton).
Credit: Hurwitz 1898; Cartan 1914; Killing; Gunaydin & Gursey 1973; Dixon; Furey 2018; Pati & Salam 1974; Georgi & Jarlskog 1979;
Georgi, Quinn & Weinberg 1974; Chang, Mohapatra & Parida 1984; Senjanovic & Mohapatra 1975; Zeldovich & Khlopov 1978; Preskill 1979;
Koide 1982; Foot 1994; Brannen 2006; Rivero 2011, 2026; Shulga 2026; Gatto, Sartori & Tonin 1968; Weinberg 1964; Llewellyn Smith 1973;
Cornwall, Levin & Tiktopoulos 1974; Coleman & Mandula 1967."""
import math, itertools, numpy as np
rng = np.random.default_rng(34)
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)

# ------------------------------------------------------------------ algebra tools (Cayley-Dickson; works with complex coefficients too)
def conj(x): return np.concatenate([x[:1], -x[1:]])
def mul(a, b):
    n = len(a)
    if n == 1: return a*b
    h = n//2; p, q, r, s = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])
def basis(n): return [np.eye(n)[i] for i in range(n)]
def derivations(n):
    """All linear maps D with D(xy) = D(x)y + xD(y): the infinitesimal 'turns that keep the multiplication'."""
    E = basis(n); T = np.array([[mul(E[i], E[j]) for j in range(n)] for i in range(n)]); rows = []
    for i in range(n):
        for j in range(n):
            for c in range(n):
                row = np.zeros((n, n))
                for k in range(n): row[c, k] += T[i, j, k]
                for a in range(n): row[a, i] -= T[a, j, c]; row[a, j] -= T[i, a, c]
                rows.append(row.ravel())
    M = np.array(rows); _, S, Vt = np.linalg.svd(M); k = int((S < 1e-9).sum()) + (M.shape[1] - len(S))
    return [Vt[-i-1].reshape(n, n) for i in range(k)]
def span_dim(mats): return int(np.linalg.matrix_rank(np.array([m.ravel() for m in mats]), tol=1e-8)) if mats else 0
def lie_closure(gens):
    cur = list(gens); d = span_dim(cur)
    while True:
        new = cur + [a@b - b@a for a in cur for b in cur]
        # keep an independent basis
        A = np.array([m.ravel() for m in new]); U, S, Vt = np.linalg.svd(A, full_matrices=False); r = int((S > 1e-8).sum())
        cur = [Vt[i].reshape(gens[0].shape) for i in range(r)]
        if r == d: return cur
        d = r
def rank_of(alg):
    """Rank = number of carriers that can all be 'neutral' together = dimension of a maximal commuting set (generic element's centraliser)."""
    if not alg: return 0
    best = 10**9
    for _ in range(4):
        X = sum(rng.normal()*g for g in alg)
        best = min(best, len(alg) - int(np.linalg.matrix_rank(np.array([(X@g - g@X).ravel() for g in alg]).T, tol=1e-8)))
    return best
def so(n):
    out = []
    for i in range(n):
        for j in range(i+1, n):
            m = np.zeros((n, n)); m[i, j] = 1; m[j, i] = -1; out.append(m)
    return out
def commutant(alg_big, alg_small):
    """Elements of alg_big that commute with everything in alg_small."""
    A = np.array([np.concatenate([(b@s - s@b).ravel() for s in alg_small]) for b in alg_big]).T
    _, S, Vt = np.linalg.svd(A); k = len(alg_big) - int((S > 1e-8).sum())
    return [sum(c*b for c, b in zip(Vt[-i-1], alg_big)) for i in range(k)]
def stabiliser(alg, vec):
    A = np.array([g@vec for g in alg]).T; _, S, Vt = np.linalg.svd(A); k = len(alg) - int((S > 1e-9).sum())
    return [sum(c*g for c, g in zip(Vt[-i-1], alg)) for i in range(k)]
def Lmat(u, n): return np.array([mul(u, e) for e in basis(n)]).T
def Rmat(u, n): return np.array([mul(e, u) for e in basis(n)]).T

head("PART 1. FORCE GROUPS: can 'all three forces follow from the bag rules' be shown to fail?")
pr("Needed (measured): 12 carriers in three blocks that do not mix: 8 (2 neutral) | 3 (1 neutral) | 1 (1 neutral). Rank 4 in total.")
pr("The 4x4 (Pati-Salam) version needs rank 5 (15 + 3 + 3 carriers).\n")
pr("1a. Every reading of 'a symmetry of the inner space', applied the SAME way to each number system the bag rules allow")
pr(f"    {'reading':58s} {'complex':>10s} {'quaternions':>12s} {'octonions':>11s}")
table = {}
for name, n in (("C", 2), ("H", 4), ("O", 8)):
    E = basis(n); der = derivations(n)
    u = np.concatenate([[0], rng.normal(size=n-1)]); u /= np.linalg.norm(u)
    der_u = stabiliser(der, u) if der else []
    son = so(n); son1 = stabiliser(son, E[0])
    J = Lmat(u, n); un = commutant(son, [J]); un1 = stabiliser(un, E[0]) if un else []
    comm = lambda alg: span_dim([a@b - b@a for a in alg for b in alg])
    gens = [Lmat(e, n) for e in E[1:]] + [Rmat(e, n) for e in E[1:]]; mult = lie_closure(gens)
    gensL = [Lmat(e, n) for e in E[1:]]; multL = lie_closure(gensL)
    table[name] = dict(der=(len(der), rank_of(der)), der_u=(len(der_u), rank_of(der_u)), son=(len(son), rank_of(son)), son1=(len(son1), rank_of(son1)),
                       un=(len(un), rank_of(un)), un1=(len(un1), rank_of(un1)), unss=comm(un), mult=(len(mult), rank_of(mult)), multL=(len(multL), rank_of(multL)))
    if name == "O": derO, derOu, uO, soO, JO, unO = der, der_u, u, son, J, un
def row(label, key): pr(f"    {label:58s} " + " ".join(f"{str(table[k][key][0])+' (rank '+str(table[k][key][1])+')':>11s}" for k in "CHO"))
row("A. turns that keep the multiplication (automorphisms)", "der")
row("B. the same, with one direction picked (the perspective)", "der_u")
row("C. turns that keep sizes (the literal content of rule 2)", "son")
row("D. the same, keeping the unit fixed", "son1")
row("E. the same, respecting the perspective's complex structure", "un")
row("E'. the same as E, also keeping the unit fixed", "un1")
row("F. multiplying by unit elements, left side only (closure)", "multL")
row("G. multiplying by unit elements, left and right (closure)", "mult")
pr("    Wanted in ONE row: 1 (rank 1) for hypercharge, 3 (rank 1) for the weak force, 8 (rank 2) for colour.")
pr("    -> No row has them. Colour's 8 (rank 2) appears only in row B on the octonions. The weak force's 3 appears in rows A, D, F on the")
pr("       quaternions. an earlier test took colour from row B and 'left and right' from row G: two different readings of rule 3.")

pr(f"    Row E in detail: the self-interacting parts have {table['H']['unss']} and {table['O']['unss']} carriers: quaternions give SU(2) x U(1) (4 carriers, 2 neutral: the measured")
pr("    electroweak count) and octonions give SU(4) x U(1) (the 4x4 colour group of an earlier test plus one). Row E is the only single reading that gets the")
pr("    weak block's count right, but it gives 16 colour-side carriers, not 8, one U(1) too many, and still needs two separate inner spaces. A lead, not a result.")
pr("    Also: the 4x4 colour group (15 carriers, rank 3) cannot sit inside the octonion's symmetries of row A (14, rank 2). So 'octonions give colour'")
pr("    (an earlier test, row B) and 'the lepton is the 4th colour' (an earlier test, which is row E) are two different structures, and both have been in use.")
pr("\n1b. Can ONE octonion inner space hold colour and the weak force side by side?")
pr(f"    Row A: the octonion's symmetries have rank {table['O']['der'][1]}; the measured forces need rank 4. Impossible, whatever the details.")
com = commutant(soO, derOu)
pr(f"    Row C (most generous, 28 turns, rank 4): turns that commute with colour's 8: {len(com)}, and they commute with each other: "
   f"{max(np.abs(a@b-b@a).max() for a in com for b in com) < 1e-9}.")
pr(f"    -> next to colour there is room for {len(com)} single (U(1)-type) carriers and no room for the weak force's 3. And rank 4 < 5: the 4x4 group cannot fit at all.")

pr("\n1c. So the weak force needs a second number system (the quaternions) next to the octonions. Does putting them together keep rule 1?")
i_ = 1j
def cplx(n, k): v = np.zeros(n, dtype=complex); v[k] = 1; return v
for name, n in (("complex x quaternions", 4), ("complex x octonions", 8)):
    a = cplx(n, 0) + i_*cplx(n, 1); b = cplx(n, 0) - i_*cplx(n, 1); p = mul(a, b)
    pr(f"    {name:22s}: (1 + i e1) x (1 - i e1) = {np.abs(p).max():.1f}   <- two non-zero things multiply to nothing: a hole in the glass")
# quaternion (x) octonion: tensor product, element = 4x8 array; product (a(x)b)(c(x)d) = ac (x) bd
def tmul(X, Y):
    EH, EO = basis(4), basis(8); out = np.zeros((4, 8))
    for i in range(4):
        for j in range(8):
            if X[i, j] == 0: continue
            for k in range(4):
                for l in range(8):
                    if Y[k, l] == 0: continue
                    out += X[i, j]*Y[k, l]*np.outer(mul(EH[i], EH[k]), mul(EO[j], EO[l]))
    return out
X = np.zeros((4, 8)); X[0, 0] = 1; X[1, 1] = 1          # 1(x)1 + i(x)e1
Y = np.zeros((4, 8)); Y[0, 0] = 1; Y[1, 1] = -1         # 1(x)1 - i(x)e1 ... (i e1)^2 = (+1): check
pr(f"    quaternions x octonions: (1 + i(x)e1) x (1 - i(x)e1) = {np.abs(tmul(X, Y)).max():.1f}")
pr("    The published octonion routes to the weak force (Dixon; Furey) work inside exactly these products and build particles from")
pr("    projectors P with P x (1 - P) = 0. They need the holes that rule 1 forbids.")
pr("    The only way left: keep three SEPARATE inner spaces at every point (complex, quaternion, octonion), never multiplied together.")
pr("    That keeps rule 1, but it is an added assumption (why three spaces, each once?), and row by row above it still needs a different")
pr("    reading of 'symmetry' for each space: unit multiplication for the first two (row F: 1 and 3) and row B for the octonions (8).")
pr("    Noted, not used: with a perspective picked in the quaternions too, row B leaves 1 carrier there and 8 in the octonions. That is the")
pr("    pattern of the two forces that stay unbroken in nature (electromagnetism and colour). Suggestive; it is not a derivation of the weak force.")
pr("    -> RESULT: 'rules 1-4 give SU(3) x SU(2) x U(1)' is FALSE as stated. Rules 1-2 allow one number system at a time;")
pr("       the octonions give colour and nothing else that commutes with it except single carriers. The weak force and hypercharge are NOT derived.")

pr("\n1d. What survives the attack")
ok = []
for t in range(5):
    u = np.concatenate([[0], rng.normal(size=7)]); u /= np.linalg.norm(u); st = stabiliser(derO, u); ok.append((len(st), rank_of(st)))
pr(f"    Colour: any of 5 new random perspectives leaves (carriers, neutral) = {sorted(set(ok))} -> SU(3), three colours. Survives.")
pr("    Uniqueness from counting: 8 carriers with 2 neutral is only SU(3); 4 with 2 neutral is only SU(2) x U(1). That is the")
pr("    Killing-Cartan classification, a theorem. It survives, but its inputs (the counts) are MEASURED, not derived from the bag rules.")
# anomaly check: the content is consistent (a test it passes)
Y = {"Q": (1/6, 3, 2), "u": (-2/3, 3, 1), "d": (1/3, 3, 1), "L": (-1/2, 1, 2), "e": (1, 1, 1), "N": (0, 1, 1)}   # left-handed Weyl fields (conjugates for right-handed)
s1 = sum(y*c*w for y, c, w in Y.values()); s3 = sum(y**3*c*w for y, c, w in Y.values())
s22 = sum(y*c for y, c, w in Y.values() if w == 2); s33 = sum(y*w for y, c, w in Y.values() if c == 3)
pr(f"    Consistency (anomalies) of the 16 particles per generation: sum Y = {s1:+.1e}, sum Y^3 = {s3:+.1e}, weak-weak-Y = {s22:+.1e}, colour-colour-Y = {s33:+.1e}: all zero. Passes.")

pr("\n1e. Left-handedness: does 'nothing extra' (rule 6, the minimal 4x4 content) survive the measured masses?")
bt = 0.734   # m_b/m_tau at the right side's scale 5.5e11 GeV (an earlier test run: 0.756 at 1e11, 0.697 at 1e13)
mc_ms, mb_mc, ms_md = 11.766, 4.577, 20.36                      # FLAG 2024 ratios (scale independent)
mtau, mmu, me = 1776.93, 105.6583755, 0.51099895
sm = bt/(mb_mc*mc_ms)*mtau/mmu; de = bt/(mb_mc*mc_ms*ms_md)*mtau/me
pr("    The minimal 4x4 Higgs content (one bidoublet) forces, at the right side's scale: down-type quark mass = charged lepton mass, generation by generation,")
pr("    and up-type quark mass = neutrino Dirac link (Pati & Salam 1974).")
pr(f"      bottom / tau     : predicted 1, measured {bt:.2f}")
pr(f"      strange / muon   : predicted 1, measured {sm:.2f}   (off by {1/sm:.1f}x)")
pr(f"      down / electron  : predicted 1, measured {de:.2f}   (off by {de:.1f}x)")
pr("    And with rule 5 (neutrino link = electron family's sizes) the same content forces top = tau: measured top/tau at that scale is about 50.")
pr("    -> The minimal content is CONTRADICTED by data and by rule 5. At least one more Higgs piece is required (the (15,2,2); Georgi & Jarlskog 1979).")
a1, a2, a3 = 1/59.00, 1/29.59, 0.1180; MZ = 91.1876; bSM = np.array([41/10, -19/6, -7])
def solve(b4, b2L, b2R):
    inv = np.array([1/a1, 1/a2, 1/a3])
    c4 = (inv[2], -bSM[2] + b4, -b4); c2L = (inv[1], -bSM[1] + b2L, -b2L)
    r0 = (5/3)*(inv[0] - 0.4*inv[2]); r1 = (5/3)*(-bSM[0] + 0.4*bSM[2]); c2R = (r0, r1 + b2R, -b2R)
    A = np.array([[c4[1]-c2L[1], c4[2]-c2L[2]], [c2L[1]-c2R[1], c2L[2]-c2R[2]]]); rhs = -np.array([c4[0]-c2L[0], c2L[0]-c2R[0]])
    tI, tU = np.linalg.solve(A, rhs); MI = MZ*math.exp(2*math.pi*tI); MU = MZ*math.exp(2*math.pi*tU)
    aU = 1/(c4[0] + c4[1]*tI + c4[2]*tU); tau = MU**4/(aU**2*0.938**5)*6.582e-25/3.156e7
    a2R = 1/(r0 + r1*tI)
    return MI, MU, tau, a2R
b_min = (-44/3 + 4 + 3, -22/3 + 4 + 1/3, -22/3 + 4 + 1/3 + 20/3)
cases = [("minimal (bidoublet + right triplet)", b_min, "handedness put in: the content itself is lopsided"),
         ("left-right symmetric (+ left triplet)", (b_min[0]+3, b_min[1]+20/3, b_min[2]), "handedness explained (spontaneous)"),
         ("minimal + (15,2,2) light", (b_min[0]+16/3, b_min[1]+5, b_min[2]+5), "fixes the mass relations; handedness put in"),
         ("symmetric + (15,2,2) light", (b_min[0]+3+16/3, b_min[1]+20/3+5, b_min[2]+5), "fixes the mass relations; handedness explained"),
         ("minimal + (15,2,2) at the top scale", b_min, "fixes the mass relations; handedness put in"),
         ("symmetric start, left triplet lifted to the top scale", b_min, "handedness explained; needs one parity-odd field (Chang, Mohapatra & Parida 1984)")]
pr("\n    Force running redone for each consistent content (one loop). Needs: right side >= heaviest partner (1.0e11-1.3e11 GeV); proton > 2.4e34 yr (Super-K).")
pr(f"    {'content':52s} {'right side':>11s} {'unify':>9s} {'proton (yr)':>12s}   verdict")
res = {}
for lab, b, note in cases:
    try:
        MI, MU, tau, a2R = solve(*b); good = MI < MU and MI > 1e3
    except Exception: good = False
    if not good or not np.isfinite(MI): pr(f"    {lab:52s} {'no solution':>11s}"); continue
    v = []
    if MI < 1.0e11: v.append("right side BELOW the heaviest partner")
    if tau < 2.4e34: v.append(f"proton likely {2.4e34/tau:.0f}x too fast")
    res[lab] = (MI, MU, tau, a2R)
    pr(f"    {lab:52s} 10^{math.log10(MI):5.2f}  10^{math.log10(MU):5.2f}  10^{math.log10(tau):5.1f}     {'; '.join(v) if v else 'OK'} | {note}")
pr("    -> Every version that passes the proton and partner tests either puts handedness in by hand or adds a field beyond 'nothing extra'.")
pr("       The fully symmetric versions are in tension with the proton limit at this rough one-loop level (a factor of about 3 in the unification")
pr("       scale; heavy-particle thresholds could move it either way, so this is 'likely too fast', not a firm exclusion).")
pr("       'Nothing left to assume beyond rule 6' is FALSE: rule 6 itself, read literally, is contradicted.")

pr("\n1f. The gauge principle 'follows from rule 2': do the theorems' conditions hold inside the bag rules?")
pr("    Counter-example: the octonion algebra at a point satisfies rules 1 and 2 exactly and contains no time, no relativity and no probabilities.")
pr("    So rules 1-2 cannot by themselves force a LOCAL symmetry; a global one satisfies them equally well.")
pr("    The theorems (Weinberg 1964; Llewellyn Smith 1973; Cornwall, Levin & Tiktopoulos 1974) assume special relativity, quantum mechanics")
pr("    (probabilities add to 1 in TIME) and, for heavy carriers, weak coupling. None of these is among rules 1-7.")
pr("    -> RESULT: the gauge principle is a theorem of relativity + quantum mechanics, which TRPG inherits. It is NOT a consequence of rule 2.")
pr("       Honest repair: list relativity and quantum mechanics as inherited starting rules. Then rule 3 can be dropped; the count of rules does not shrink.")

pr("\n1g. A condition the 4x4 picture brings with it (found while attacking it)")
MI = res["minimal (bidoublet + right triplet)"][0]; a2R = res["minimal (bidoublet + right triplet)"][3]
mM = MI/a2R; nos = 1e2*(MI/1.22e19)**3; Om = mM*nos*2891.2/1.0537e-5
pr(f"    Breaking the 4x4 group makes magnetic monopoles of mass about {mM:.1e} GeV. If the universe was ever hotter than the break ({MI:.1e} GeV)")
pr(f"    after it happened, the standard estimate (Zeldovich & Khlopov 1978; Preskill 1979) gives a monopole density of about {Om:.0f}x the critical")
pr("    density (allowed: under 0.12). So TRPG needs: never hotter than ~5e11 GeV after the break, yet at least ~1e11 GeV for the matter excess")
pr("    (falsification_matter_excess, E2). A window of about 5x. Rough estimate; it is a condition, not a failure.")

head("PART 2. QUARK MASSES: the colour rule, the bottom mass and charm")
def params(roots):
    roots = np.array(roots, float); A = roots.mean(); M = A*A; x = roots/A - 1; r = math.sqrt(2*np.mean(x*x))
    k = int(np.argmax(roots)); d = None
    c3 = 4*np.prod(x)/r**3            # cos(3 delta) = 4 prod(x_k)/r^3
    d = math.acos(max(-1, min(1, c3)))/3
    return M, r, d
def family(M, delta, r=math.sqrt(2)): return sorted([M*(1 + r*math.cos(delta + 2*math.pi*k/3))**2 for k in range(3)])
def Qv(roots): roots = np.array(roots); return float((roots**2).sum()/roots.sum()**2)
lep = [me, mmu, mtau]; Ml, rl, dl = params(np.sqrt(lep))
pr(f"Lepton family (measured): scale {Ml:.3f} MeV, lean r = {rl:.6f} (sqrt2 = {math.sqrt(2):.6f}), twist = {dl:.7f} (2/9 = {2/9:.7f})")
pred = family(3*Ml, 3*dl); ps, pc, pb = pred
def koide_third(m1, m2):
    a, b = math.sqrt(m1), math.sqrt(m2); B = -4*(a+b); C = 3*(m1+m2) - 2*(a+b)**2
    return ((-B + math.sqrt(B*B-4*C))/2)**2, ((-B - math.sqrt(B*B-4*C))/2)**2
pt = koide_third(pc, pb)[0]/1000
data = {"strange": (93.5, 0.8, ps), "charm": (1273.0, 4.6, pc), "bottom": (4183.0, 7.0, pb), "top (direct, GeV)": (172.57, 0.29, pt)}
pr("\n2a. The colour rule taken as an exact rule, against each measurement's own error bar")
chi_all = chi_noc = 0
for k, (m, s, p) in data.items():
    z = (p - m)/s; chi_all += z*z; chi_noc += z*z if k != "charm" else 0
    pr(f"    {k:18s} predicted {p:9.2f}   measured {m:8.2f} +- {s:<5.2f}  off {100*(p-m)/m:+5.2f}%  = {z:+5.1f} sigma")
from math import erfc
def pval(chi2, n):
    from scipy.stats import chi2 as c2; return float(c2.sf(chi2, n))
pr(f"    All four together: chi2 = {chi_all:.0f} for 4 numbers: probability {pval(chi_all, 4):.0e}. Without charm: chi2 = {chi_noc:.1f} for 3: probability {pval(chi_noc, 3):.3f}.")
pr("    -> As an exact rule it FAILS: charm alone is 19 sigma out, and even without charm the other three sit 1.5-2.4 sigma on the wrong side.")
pr("       The bottom cannot be kept while charm is dropped: both come out of the same two numbers (3 x scale, 3 x twist).")

pr("\n2b. Is there ANY single, consistent way of defining quark masses in which the rule holds?")
mZ = dict(t=168.26e3, b=2839.0, c=620.0); mZ["s"] = mZ["c"]/mc_ms
defs = {"mixed, as used in an earlier test (s at 2 GeV, c and b at own scale, t direct)": dict(s=93.5, c=1273.0, b=4183.0, t=172.57e3),
        "own scale for all (t = 162.5 GeV running)": dict(s=93.5, c=1273.0, b=4183.0, t=162.5e3),
        "all at 2 GeV (ratios from lattice QCD)": dict(s=93.5, c=93.5*mc_ms, b=93.5*mc_ms*mb_mc, t=None),
        "all at the bottom's scale": dict(s=4183/(mb_mc*mc_ms), c=4183/mb_mc, b=4183.0, t=None),
        "all at M_Z (Huang & Zhou 2021)": mZ,
        "pole masses (c 1.67, b 4.78, t 172.57 GeV)": dict(s=None, c=1670.0, b=4780.0, t=172.57e3)}
pr(f"    {'definition':66s} {'Q(-s,c,b)':>9s} {'lean r':>7s} {'twist/(2/9)':>11s} {'scale/lepton':>12s} {'Q(c,b,t)':>9s}")
for k, d in defs.items():
    if d["s"]:
        M, r, dd = params([-math.sqrt(d["s"]), math.sqrt(d["c"]), math.sqrt(d["b"])]); q1 = Qv([-math.sqrt(d["s"]), math.sqrt(d["c"]), math.sqrt(d["b"])])
        a = f"{q1:9.4f} {r:7.4f} {dd/(2/9):11.3f} {M/Ml:12.3f}"
    else: a = f"{'-':>9s} {'-':>7s} {'-':>11s} {'-':>12s}"
    q2 = f"{Qv([math.sqrt(d['c']), math.sqrt(d['b']), math.sqrt(d['t'])]):9.4f}" if d["t"] else f"{'-':>9s}"
    pr(f"    {k:66s} {a} {q2}")
pr(f"    {'THE RULE SAYS':66s} {2/3:9.4f} {math.sqrt(2):7.4f} {3.0:11.3f} {3.0:12.3f} {2/3:9.4f}")
# Monte Carlo errors for the common-scale shape (scale independent) and for the mixed convention
N = 200000
r1 = rng.normal(mc_ms, 0.030, N); r2 = rng.normal(mb_mc, 0.008, N)
roots = np.array([-np.ones(N), np.sqrt(r1), np.sqrt(r1*r2)]); Qc = (roots**2).sum(0)/roots.sum(0)**2
pr(f"    At ANY one common energy the shape is fixed by two lattice ratios: Q(-s,c,b) = {Qc.mean():.4f} +- {Qc.std():.4f}, against 2/3 = 0.6667: {(Qc.mean()-2/3)/Qc.std():.0f} sigma out.")
s_ = rng.normal(93.5, 0.8, N); c_ = rng.normal(1273, 4.6, N); b_ = rng.normal(4183, 7, N)
roots = np.array([-np.sqrt(s_), np.sqrt(c_), np.sqrt(b_)]); A = roots.mean(0); x = roots/A - 1; rr = np.sqrt(2*(x*x).mean(0)); dd = np.arccos(np.clip(4*x.prod(0)/rr**3, -1, 1))/3
pr(f"    In the mixed convention: scale/lepton = {(A*A/Ml).mean():.4f} +- {(A*A/Ml).std():.4f} ({((A*A/Ml).mean()-3)/(A*A/Ml).std():+.0f} sigma from 3); "
   f"twist/(2/9) = {(dd/(2/9)).mean():.4f} +- {(dd/(2/9)).std():.4f} ({((dd/(2/9)).mean()-3)/(dd/(2/9)).std():+.0f} sigma from 3); lean = {rr.mean():.4f} +- {rr.std():.4f} ({(rr.mean()-math.sqrt(2))/rr.std():+.0f} sigma from sqrt2).")
pr("    -> No consistent definition makes the rule hold. It comes closest only when each quark is quoted in a different convention,")
pr("       and even there every one of its three numbers (3, 3, sqrt2) is many sigma away. RESULT: the colour rule FAILS as a law.")
pr("       It stays what Rivero called it in 2011: a numerical observation, good to a few per cent.")

pr("\n2c. Charm: can any standard definition rescue it?")
want = {"from (b, t), 45-degree rule": koide_third(172.57e3, 4183.0)[1], "from the colour rule": pc}
sb = None
# charm wanted by (s, b) with Q(-s, c, b) = 2/3
from scipy.optimize import brentq
want["from (s, b), 45-degree rule"] = brentq(lambda c: Qv([-math.sqrt(93.5), math.sqrt(c), math.sqrt(4183.0)]) - 2/3, 800, 2500)
for k, v in want.items(): pr(f"    charm wanted {k:30s}: {v:7.1f} MeV")
pr("    charm in standard definitions: 986 (running, at 3 GeV), 1100 (at 2 GeV), 1273 +- 5 (at its own scale), 1548 (half the J/psi), 1670 +- 70 (pole).")
pr("    The wanted 1353-1372 matches none of them. A running mass equals 1360 at an energy of about 1.07 GeV, but nothing in TRPG picks that energy;")
pr("    since some definition can be found for ANY number between 0.99 and 1.67 GeV, 'a different definition' is an escape that cannot fail. Rejected.")
pr(f"    -> RESULT: charm FAILS ({(want['from (b, t), 45-degree rule']-1273)/4.6:.0f} to {(pc-1273)/4.6:.0f} sigma). The heavy-quark 45-degree relation is approximate (6-8%), not a law.")
pr("       The 'one leftover up-side factor' lead of an earlier test is withdrawn: the top value in that comparison is computed FROM charm, so it is not a second witness.")

pr("\n2d. Every link of the quark chain, each against its own error bar")
vus, dvus = 0.2243, 0.0008
st = brentq(lambda s: Qv([-math.sqrt(s), math.sqrt(1273.0), math.sqrt(4183.0)]) - 2/3, 20, 200)
st_err = np.std([brentq(lambda s: Qv([-math.sqrt(s), math.sqrt(c), math.sqrt(b)]) - 2/3, 20, 200) for c, b in zip(rng.normal(1273, 4.6, 400), rng.normal(4183, 7, 400))])
tt = koide_third(1273.0, 4183.0)[0]/1000
tt_err = np.std([koide_third(c, b)[0]/1000 for c, b in zip(rng.normal(1273, 4.6, 400), rng.normal(4183, 7, 400))])
links = [("top = full Higgs link (y_t = 1)", 174.10, 0.0, 172.57, 0.29, "GeV"),
         ("top from measured (c, b), 45 degrees", tt, tt_err, 172.57, 0.29, "GeV"),
         ("charm from measured (b, t), 45 degrees", want["from (b, t), 45-degree rule"], 3.0, 1273.0, 4.6, "MeV"),
         ("strange from measured (c, b), 45 degrees", st, st_err, 93.5, 0.8, "MeV"),
         ("strange/down from the Cabibbo cell (1/Vus^2)", 1/vus**2, 2*dvus/vus**3, ms_md, 0.27, ""),
         ("bottom from the colour rule", pb, 0.0, 4183.0, 7.0, "MeV")]
for lab, p, dp, m, dm, u in links:
    z = (p - m)/math.hypot(dp, dm)
    pr(f"    {lab:46s} predicted {p:8.2f} +- {dp:5.2f}  measured {m:8.2f} +- {dm:4.2f} {u:3s}  off {100*(p-m)/m:+6.1f}%  {z:+5.1f} sigma   {'PASSES' if abs(z) < 2 else ('borderline' if abs(z) < 3 else 'FAILS')}")
pr("    (Up quark from the eta -> 3 pi cell: 2.03 MeV, range 1.72-2.28, measured 2.16 +- 0.07: passes, but the test is weak.)")
sc = brentq(lambda s: Qv([-math.sqrt(s), math.sqrt(want["from (b, t), 45-degree rule"]), math.sqrt(4183.0)]) - 2/3, 20, 200)
pr(f"    Why the chain looked better before: it fed its own too-high charm (1353) into the strange link, which then gives {sc:.1f} MeV. With the measured")
pr("    charm the same link gives 84.5. The good strange value rested on the wrong charm.")
pr("    -> Strictly: the Cabibbo link passes. Everything that runs through the 45-degree rule for quarks fails or is borderline. The bottom stays an INPUT.")

head("PART 3. THE 2/9 TWIST")
pr("3a. As a measured rule for the charged leptons (their masses are directly measurable, so there is no definition problem)")
e_ = rng.normal(me, 1.6e-10, N); m_ = rng.normal(mmu, 2.3e-6, N); t_ = rng.normal(mtau, 0.09, N)
roots = np.sqrt(np.array([e_, m_, t_])); A = roots.mean(0); x = roots/A - 1; rr = np.sqrt(2*(x*x).mean(0)); dd = np.arccos(np.clip(4*x.prod(0)/rr**3, -1, 1))/3
Qm = (roots**2).sum(0)/roots.sum(0)**2
pr(f"    all three masses measured: Q = {Qm.mean():.7f} +- {Qm.std():.7f} (2/3 = 0.6666667: {(Qm.mean()-2/3)/Qm.std():+.1f} sigma);  twist = {dd.mean():.7f} +- {dd.std():.7f} (2/9 = 0.2222222: {(dd.mean()-2/9)/dd.std():+.1f} sigma)")
pr("    -> each rule on its own agrees with data inside the tau mass's error bar.")
pr("\n3b. The sharper test: BOTH rules exact at once (45 degrees AND 2/9). Then two mass ratios are predicted with nothing free.")
f = family(1.0, 2/9); re_mu = f[0]/f[1]; rt_mu = f[2]/f[1]
pr(f"    electron/muon : predicted {re_mu:.10f}   measured {me/mmu:.10f} +- {me/mmu*2.2e-8:.1e}   off {1e6*(re_mu/(me/mmu)-1):+.2f} parts per million = {(re_mu-me/mmu)/(me/mmu*2.2e-8):+.0f} sigma")
pr(f"    tau/muon      : predicted {rt_mu:.7f}      measured {mtau/mmu:.7f} +- {0.09/mmu:.7f}       off {1e6*(rt_mu/(mtau/mmu)-1):+.1f} parts per million = {(rt_mu-mtau/mmu)/(0.09/mmu):+.1f} sigma")
d_em = brentq(lambda d: family(1.0, d)[0]/family(1.0, d)[1] - me/mmu, 0.20, 0.24)
pr(f"    twist needed by the electron/muon ratio if 45 degrees is exact: {d_em:.11f}; 2/9 = {2/9:.11f}; difference {d_em-2/9:+.2e} ({1e6*(d_em/(2/9)-1):+.2f} ppm)")
mt_k = family(1.0, d_em); pr(f"    (with that twist the tau comes out at {mmu*mt_k[2]/mt_k[1]:.3f} MeV; measured 1776.93 +- 0.09)")
r_need = brentq(lambda r: family(1.0, 2/9, r)[0]/family(1.0, 2/9, r)[1] - me/mmu, 1.40, 1.43); fr = family(1.0, 2/9, r_need)
pr(f"    lean needed instead if 2/9 is exact: r = {r_need:.9f} (sqrt2 = {math.sqrt(2):.9f}), i.e. Q = {(1+r_need**2/2)/3:.9f}; tau then {mmu*fr[2]/fr[1]:.3f} MeV.")
pr("    Both alternatives fit the tau mass, so the data cannot say WHICH of the two is inexact, only that they are not both exact.")
pr("    -> RESULT: 'exactly 45 degrees and exactly 2/9' is RULED OUT by the electron/muon mass ratio, by about 8-9 parts per million")
pr("       (hundreds of sigma, because that ratio is known to 0.02 ppm). At least one of the two is approximate at the 1e-6 level.")
pr("       So 2/9 is a very good approximation (better than 1 in 100,000), not an exact number of nature in this form.")
pr("\n3c. Does the rule hold at other energies (running masses)?")
lz = [0.48307, 101.766, 1728.56]   # MeV at M_Z (Huang & Zhou 2021)
Mz_, rz, dz = params(np.sqrt(lz))
pr(f"    at M_Z: Q = {Qv(np.sqrt(lz)):.6f}, twist = {dz:.6f} (2/9 = 0.222222; off {1e2*(dz/(2/9)-1):+.3f}%). The rule is a statement about the directly measured masses only.")
pr("\n3d. The ring reading (two paths round a ring of 3 seats: (1/3)(2/3) = 2/9; Shulga 2026). What does it predict that could fail?")
ssum = sum((1-math.cos(2*math.pi*n/3))/n**2 for n in range(1, 2000000))/math.pi**2
pr(f"    ring sum = {ssum:.7f} (2/9 = {2/9:.7f}): the arithmetic is right.")
pr(f"    (i) Unit. A ring's natural unit is the full turn. 2/9 of a turn = {2/9*2*math.pi:.4f} rad; the measured twist is {dl:.4f} rad = {dl/(2*math.pi):.5f} of a turn.")
pr("        The reading only works if 'one ring = one radian', which nothing gives. For leptons it therefore predicts nothing: it was matched to the number.")
M_, r_, dq = params([-math.sqrt(93.5), math.sqrt(1273.0), math.sqrt(4183.0)]); roots = np.array([-np.sqrt(s_), np.sqrt(c_), np.sqrt(b_)]); A = roots.mean(0); x = roots/A - 1; rr = np.sqrt(2*(x*x).mean(0)); dqs = np.arccos(np.clip(4*x.prod(0)/rr**3, -1, 1))/3
M2, r2_, dq2 = params([-1.0, math.sqrt(mc_ms), math.sqrt(mc_ms*mb_mc)])
pr(f"    (ii) Its one extension, quarks = 3 x 2/9 = 0.6667: measured {dq:.4f} +- {dqs.std():.4f} in the mixed convention ({(dq-2/3)/dqs.std():+.0f} sigma); {dq2:.4f} at any common energy. FAILS.")
nu = [0.00036, 0.0087, 0.0501]
for sgn in (1, -1):
    Mn, rn, dn = params([sgn*math.sqrt(nu[0]), math.sqrt(nu[1]), math.sqrt(nu[2])])
    pr(f"    (iii) Neutrinos (TRPG's masses, lightest root sign {sgn:+d}): lean {rn:.4f}, twist {dn:.4f} = 2/9 {dn-2/9:+.4f}. Not 2/9, and the ring gives no reason for the difference.")
pr("    -> RESULT: the ring reading is a way of writing 2/9, not an explanation. It has a free unit, predicts nothing new for leptons,")
pr("       and its only extension fails. It is WITHDRAWN as an explanation. The twist stays one measured number.")
