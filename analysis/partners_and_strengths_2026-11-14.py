"""Are the two missing partners and the force strengths one problem? (14 Nov 2026)
Main line: the relationship-matrix picture. Uses the tools of forces_from_entanglement_2026-11-14.py (same folder).
Credit: Georgi & Glashow 1974; Georgi, Quinn & Weinberg 1974; Pati & Salam 1974; Chamseddine & Connes 1996; Baez & Huerta 2010; Baez & Schwahn 2026; Baez 2026;
Dubois-Violette & Todorov 2018-19; Krasnov 2019; Jordan, von Neumann & Wigner 1934; Albert 1934; Machacek & Vaughn 1983; Zee 1980; Dynkin 1952; PDG 2024."""
import math, itertools, os, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
hs = open(os.path.join(HERE, "handedness_and_perspective_2026-11-14.py")).read()
exec(hs[hs.index("import math, itertools"):hs.index('head("1. The link states have only one')])

head("1. Every particle of a generation as a list of five switches (3 colours, 2 weak directions)")
pr("  Known since Georgi & Glashow 1974 (Baez & Huerta 2010 for this way of saying it): each left-handed particle is an even number of the five")
pr("  basic directions switched on. Charge = -1/3 per colour switched on, +1/2 per weak direction switched on.")
rows = {}
for nc in range(4):
    for nw in range(3):
        if (nc + nw) % 2 == 0:
            y = -nc/3 + nw/2; mult = math.comb(3, nc)*math.comb(2, nw)
            rows[(nc, nw)] = (y, mult)
name = {(0, 0): "hidden neutrino partner", (2, 0): "UP partner", (1, 1): "quark pair", (0, 2): "ELECTRON partner", (3, 1): "lepton pair", (2, 2): "down partner"}
inm = {(0, 0): "yes (one real neutral state)", (2, 0): "NO", (1, 1): "yes (links)", (0, 2): "NO", (3, 1): "yes (links)", (2, 2): "yes (third link)"}
pr(f"    {'colours on':>10} {'weak on':>8} {'charge':>8} {'states':>7}   {'particle':<26} in the matrix?")
for k_ in sorted(rows, key=lambda t: (t[1] != 1, t)):
    pr(f"    {k_[0]:>10} {k_[1]:>8} {rows[k_][0]:>+8.3f} {rows[k_][1]:>7}   {name[k_]:<26} {inm[k_]}")
pr(f"  Total {sum(m for _, m in rows.values())} states. The matrix holds every state with exactly ONE weak direction on (the pairs), and the down partner.")
pr("  The two it lacks are the plainest ones: two colours on and nothing else (charge 2/3 in size), and both weak directions on and nothing else (charge 1).")
pr("  Together they are 3 + 1 = 4 states: exactly one more octonion's worth (an octonion is 1 + 3 when seen facing one way).")

head("2. What the matrix's own geometry says about the three strengths")
weak = derived(weak4)
def onb(gens, target):                    # orthonormalise under -trace on the links, each scaled to 'target'
    out = []
    for g in gens:
        w_ = g.copy()
        for v_ in out: w_ = w_ - (-np.trace(V(w_)@V(v_)))/(-np.trace(V(v_)@V(v_)))*v_
        out.append(w_)
    return [w_*math.sqrt(target/(-np.trace(V(w_)@V(w_)))) for w_ in out]
wk = onb(span_basis(weak), 4.0); col = onb(span_basis(colour), 2.0)
fY = 1/(math.sqrt(-np.trace(Pq@Y@Y)/12)*6); Yf = cen[0]*fY
Cw_links = -sum(V(w_)@V(w_) for w_ in wk); Cc_links = -sum(V(c)@V(c) for c in col)
pr(f"  Checks on the links: weak j(j+1) = {np.unique(np.round(np.linalg.eigvalsh((Cw_links+Cw_links.T)/2), 4)).tolist()} (3/4 = pairs); colour = {np.unique(np.round(np.linalg.eigvalsh((Cc_links+Cc_links.T)/2), 4)).tolist()} (4/3 = three colours).")
kap = lambda g: -np.trace(g@g)             # size of a turn, measured on the whole 27-component matrix
kT, kC, kY = kap(wk[0]), kap(col[0]), kap(Yf)
s2_geo = 1/(1 + kY/kT); r32_geo = kT/kC
pr(f"  Size of one weak turn, one colour turn and the hypercharge turn, on the whole matrix: {kT:.4f}, {kC:.4f}, {kY:.4f}.")
pr(f"  If the whole matrix has ONE strength: weak mixing number sin^2 = {s2_geo:.3f}; strong strength / weak strength = {r32_geo:.3f}.")
pr("  (an earlier test quoted the 0.750. The second number is added here: the colour turns and the weak turns are not the same size in the matrix.)")

head("3. What counting the matter says")
def count(states):   # states: (colour dim, weak dim, charge); returns sums of T3^2, colour^2, Y^2
    t = sum(c*(w*w - 1)*w/12 for c, w, y in states); cc = sum(w*(0.5 if c == 3 else 0) for c, w, y in states); yy = sum(c*w*y*y for c, w, y in states)
    return t, cc, yy
pairs = [(1, 2, -1/2), (3, 2, 1/6)]; gen = pairs + [(3, 1, -2/3), (3, 1, 1/3), (1, 1, 1.0), (1, 1, 0.0)]
for nm, st in (("pairs only (what the matrix holds)", pairs), ("pairs + down partner", pairs + [(3, 1, 1/3)]), ("whole generation", gen)):
    t, cc, yy = count(st); pr(f"    {nm:<36}: weak {t:.3f}, colour {cc:.3f}, hypercharge {yy:.3f} -> sin^2 = {1/(1 + yy/t):.3f}; strong/weak = {t/cc:.3f}")
pr("  If a force is as strong as the sum of everything it acts on (the rule Chamseddine & Connes 1996 use; no larger force group is needed),")
pr("  the pairs alone give 0.750 and the whole generation gives 0.375 with strong = weak. The difference between the two IS the missing partners.")

head("4. Against the measured strengths (two-loop running of the three forces, Standard Model particles only)")
b1 = np.array([41/10, -19/6, -7.0]); b2 = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12.0], [11/10, 9/2, -26.0]])
def run_up(inv0, t1):
    f = lambda t, x: -b1/(2*math.pi) - b2@(1/x)/(8*math.pi**2)
    return solve_ivp(f, [0, t1], inv0, rtol=1e-9, dense_output=True).sol
aem, s2z, a3z = 1/127.95, 0.23122, 0.1180
inv0 = np.array([(3/5)*(1 - s2z)/aem, s2z/aem, 1/a3z]); sol = run_up(inv0, math.log(1e19/91.1876))
s2_at = lambda t: (1/sol(t)[1])**-1/((1/sol(t)[1])**-1 + (5/3)*sol(t)[0]) if False else sol(t)[1]/(sol(t)[1] + (5/3)*sol(t)[0])   # sin^2 = a_Y/(a_Y + a_2)
s2_at = lambda t: (3/5)/sol(t)[0]/((3/5)/sol(t)[0] + 1/sol(t)[1])
E_ = lambda t: 91.1876*math.exp(t); tmax = math.log(1e19/91.1876)
t23 = brentq(lambda t: sol(t)[1] - sol(t)[2], 1, tmax); t12 = brentq(lambda t: sol(t)[0] - sol(t)[1], 1, tmax)
pr(f"  Strong = weak at 10^{math.log10(E_(t23)):.2f} GeV. There the measured weak mixing number is {s2_at(t23):.3f}.")
pr(f"     pairs only predicts 0.750 ({(0.75/s2_at(t23)-1)*100:+.0f}%); whole generation predicts 0.375 ({(0.375/s2_at(t23)-1)*100:+.0f}%).")
pr(f"  The measured mixing number equals 0.375 at 10^{math.log10(E_(t12)):.2f} GeV. There strong/weak = {sol(t12)[1]/sol(t12)[2]:.3f} (whole generation predicts 1).")
# the classic strict test: assume exact meeting of all three, predict the low-energy mixing number from the other two
def predict():
    def miss(tX):
        inv = sol(tX)[1]                                           # common inverse strength taken where weak sits; run 1 down from it
        f = lambda t, x: -b1/(2*math.pi) - b2@(1/x)/(8*math.pi**2)
        return inv
    tX = t23; common = sol(tX)[1]
    f = lambda t, x: -b1/(2*math.pi) - b2@(1/x)/(8*math.pi**2)
    down = solve_ivp(f, [tX, 0], [common, common, common], rtol=1e-9).y[:, -1]
    return (3/5)/down[0]/((3/5)/down[0] + 1/down[1])
s2pred = predict()
pr(f"  Strict test (whole-generation count, all three equal where strong = weak): predicted sin^2 at the Z = {s2pred:.4f}; measured {s2z:.5f} +/- 0.00004. Off by {(s2pred/s2z-1)*100:+.1f}%.")
tg = None
g = lambda t: sol(t)[1]/sol(t)[2] - r32_geo
if g(0.01)*g(tmax) < 0:
    tg = brentq(g, 0.01, tmax); pr(f"  The matrix's own geometry (strong/weak = {r32_geo:.2f}, sin^2 = {s2_geo:.2f}): the strength ratio is met at 10^{math.log10(E_(tg)):.2f} GeV, where sin^2 is {s2_at(tg):.3f}, not {s2_geo:.2f}.")
else:
    pr(f"  The matrix's own geometry (strong/weak = {r32_geo:.2f}): that ratio is never met below 10^19 GeV (measured ratio runs from {sol(0)[1]/sol(0)[2]:.2f} down to {sol(tmax)[1]/sol(tmax)[2]:.2f}).")
pr(f"  The measured mixing number never reaches 0.75 below 10^19 GeV (it is {s2_at(tmax):.3f} there).")
pr("  -> pairs only: FAILS badly. Whole generation: the classic near miss of Georgi, Quinn & Weinberg, FAILS strictly. The partners close most of the gap, not all.")

head("5. Where could the two missing partners live? Four homes tested")
pr("5a. As a relationship between two links (the product of two link states).")
pr("    Charges do add up: lepton pair x lepton pair -> 1; and two down partners -> 2/3. But the product of two spinning particles does not spin by a half:")
pr("    it is a force-like object, not a matter particle (a charged singlet of this kind is Zee's 1980 scalar). And all pairs of links together give")
n_extra = math.comb(8, 2)
pr(f"    {n_extra} combinations where 1 (or 4) is wanted. FAILS: wrong spin, and far too many (rule 6).")
pr("5b. Among the turns of the matrix that are NOT force carriers (52 - 12 = 40).")
B = span_basis(der); Bm = np.array([g.ravel() for g in B]); Bq, _ = np.linalg.qr(Bm.T); Bq = Bq.T     # orthonormal basis of the 52 turns
def ad(X): return np.array([[Bq[i]@(X@Bq[j].reshape(dim, dim) - Bq[j].reshape(dim, dim)@X).ravel() for j in range(len(Bq))] for i in range(len(Bq))])
aw = [ad(w_) for w_ in wk]; ac = [ad(c) for c in col]; aY = ad(Yf)
CW = -sum(a@a for a in aw); CC = -sum(a@a for a in ac); Y2 = -aY@aY
M = (CW + CW.T)/2 + math.pi*(CC + CC.T)/2 + math.e**2*(Y2 + Y2.T)/2; ev, vec = np.linalg.eigh(M)
lab = {}
for i in range(len(ev)):
    v_ = vec[:, i]; key = (round(float(v_@CW@v_), 2), round(float(v_@CC@v_), 2), round(math.sqrt(abs(float(v_@Y2@v_))), 3)); lab[key] = lab.get(key, 0) + 1
pr("    weak j(j+1), colour, charge size : how many of the 52")
for k_, m_ in sorted(lab.items()): pr(f"      {k_[0]:>5}  {k_[1]:>5}  {k_[2]:>6} : {m_}")
has23 = any(abs(k_[2] - 2/3) < 1e-2 and abs(k_[0]) < 1e-2 for k_ in lab); has1 = any(abs(k_[2] - 1) < 1e-2 and abs(k_[1]) < 1e-2 and abs(k_[0]) < 1e-2 for k_ in lab)
pr(f"    A coloured weak-single turn of charge 2/3: {'present' if has23 else 'absent'}. A colourless weak-single turn of charge 1: {'present' if has1 else 'absent'}.")
pr("    And a turn is a force-like object in any case. FAILS as a home for both partners.")
pr("5c. In the doubled link space (the links written with an outside i).")
pr("    Under the 12 carriers every link state is a weak pair, so doubling gives pairs and mirror pairs, never singles (an earlier test: the hole).")
pr("    The published way to get a whole generation there (Baez 2026; Dubois-Violette & Todorov 2019) places the weak force differently and joins complex")
pr("    numbers to the octonions. Under rule 1 as TRPG has used it that has holes. FAILS under the current reading.")
pr("5d. In a larger matrix (4x4), which would have room for more links.")
rng = np.random.default_rng(41)
def rnd(N):
    X = np.zeros((N, N, 8))
    for i in range(N):
        X[i, i, 0] = rng.normal()
        for j in range(i + 1, N): X[i, j] = rng.normal(size=8); X[j, i] = conj(X[i, j])
    return X
def mm(A, Bx):
    N = A.shape[0]; C = np.zeros_like(A)
    for i in range(N):
        for j in range(N):
            for k in range(N): C[i, j] += mul(A[i, k], Bx[k, j])
    return C
jp = lambda A, Bx: (mm(A, Bx) + mm(Bx, A))/2
for N in (3, 4):
    X, Yx = rnd(N), rnd(N); X2 = jp(X, X); d = np.abs(jp(jp(X2, Yx), X) - jp(X2, jp(Yx, X))).max()
    pr(f"    {N}x{N} octonion matrices: (x.x) and x fail to work together by {d:.1e} {'(consistent)' if d < 1e-9 else '(NOT consistent: powers of a relationship become ambiguous)'}")
pr("    FAILS: the 3x3 is the largest consistent octonion matrix (Albert 1934; Jordan, von Neumann & Wigner 1934), as an earlier test already used.")

head("6. Summary")
pr("  The two open items ARE one problem: the matrix holds half of each generation's weak singles, and the strength count is wrong by exactly that half.")
pr(f"  Supplying the partners moves the strength test from {(0.75/s2_at(t23)-1)*100:+.0f}% to {(0.375/s2_at(t23)-1)*100:+.0f}% at the meeting energy ({(s2pred/s2z-1)*100:+.1f}% at the Z).")
pr("  No home for the partners that respects rule 1 was found in four attempts. What is missing is now exact: one more octonion's worth of states,")
pr("  weak singles, charges 1 (colourless) and 2/3 (three colours).")
