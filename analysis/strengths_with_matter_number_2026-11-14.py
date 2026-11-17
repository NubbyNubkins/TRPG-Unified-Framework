"""Does the matter-number force of an earlier test (the 13th turn, a real force above the hidden partner's scale in
reading (ii)) close the 10.6% strength gap of earlier tests? One strength for all forces at the meeting (sin^2 = 3/8, strong = weak), as the main line's
geometry gives. Inputs: electric strength and strong strength at the Z; predicted: the weak angle at the Z.
Credit: Georgi, Quinn & Weinberg 1974; Machacek & Vaughn 1983; Davidson 1979; Marshak & Mohapatra 1980; Mohapatra & Senjanovic 1980;
Holdom 1986 (mixing of two U(1) forces); del Aguila, Coughlan & Quiros 1988; PDG 2024."""
import math, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
mZ = 91.1876; ainv = 127.951; a3inv = 1/0.1180; s2_meas = 0.23122
B_SM = np.array([41/10, -19/6, -7.0])
def predict(extra=None):
    """One loop. Choose M_U where weak = strong; the electric strength is fixed. Returns predicted sin^2 at the Z and M_U.
    extra = (M_X, db1): above M_X the combined hypercharge running gains db1 (GUT normalisation)."""
    # unknowns: a2inv(mZ) and L = ln(M_U/mZ). Conditions: a1inv(M_U) = a2inv(M_U) = a3inv(M_U); a1inv(mZ) = (5/3)(ainv - a2inv(mZ)).
    def resid(a2):
        L = 2*math.pi*(a2 - a3inv)/(B_SM[1] - B_SM[2]); aU = a2 - B_SM[1]/(2*math.pi)*L
        a1 = (3/5)*(ainv - a2); run1 = B_SM[0]*L
        if extra and math.log(extra[0]/mZ) < L: run1 += extra[1]*(L - math.log(extra[0]/mZ))
        return a1 - run1/(2*math.pi) - aU, L
    lo, hi = 20.0, 40.0
    for _ in range(200):
        mid = (lo + hi)/2; (lo, hi) = (mid, hi) if resid(lo)[0]*resid(mid)[0] > 0 else (lo, mid)
    a2 = (lo + hi)/2; L = resid(a2)[1]; return a2/ainv, mZ*math.exp(L)
head("1. The gap as it stands (one loop; two loops give 0.2067, an earlier test)")
s, MU = predict(); pr(f"  Predicted weak angle at the Z: sin^2 = {s:.4f}; measured {s2_meas}; off by {100*(s/s2_meas - 1):+.1f}%.  Meeting energy {MU:.1e} GeV.")
head("2. Add the matter-number force above its settling energy")
pr("  What changes in the running. Above the settling energy the two single turns (hypercharge and matter number) both act. The forces' strengths")
pr("  run because of the particles they touch:")
pr("   - the matter particles: identical. Over a whole generation (hidden partner included) hypercharge and the matter-number direction orthogonal")
pr("     to it share no cross term, so their combined running equals hypercharge's alone. Checked over the 16 states:")
st = [(6, 1/6, 1/3, 0), (3, -2/3, -1/3, -0.5), (3, 1/3, -1/3, 0.5), (2, -0.5, -1, 0), (1, 1, 1, 0.5), (1, 0, 1, -0.5)]   # mult, Y, B-L, T3R (partners as antiparticles)
pr(f"       sum over states of (right-handed charge x matter number) = {sum(m*t*bl for m, y, bl, t in st):+.3f}  (0 means no cross term)")
pr("   - the matter-number carrier itself: a single (abelian) carrier does not change any running.")
pr("   - the piece that settles it (the hidden partner's pairing, matter number 2, right-handed charge 1): adds (3/5)(1/3)(1) + (2/5)(1/3)(3/8)(4) = 2/5")
pr("     to the combined hypercharge running, and nothing to the strong or weak running.")
for MX in (7.3e5, 6.3e10, 1e12, 1.2e15):
    s_, MU_ = predict((MX, 0.4)); pr(f"     settling at {MX:8.1e} GeV: sin^2 = {s_:.4f} (was {s:.4f}); off by {100*(s_/s2_meas - 1):+.1f}%")
pr("  The gap gets WIDER, not narrower, whatever the settling energy.")
head("3. Why no extra single force can close it (the attack, made general)")
pr("  The prediction misses because hypercharge, run up from the Z, reaches the weak force's strength too early. Closing the gap needs hypercharge")
pr("  to run MORE SLOWLY, or the weak force to run more slowly, or the strong force faster. A single (abelian) carrier adds nothing; every new")
pr("  piece of matter or settling field adds a POSITIVE amount to the running of every force it touches. A field touching only hypercharge and")
pr("  matter number can only speed up hypercharge. So no extension by single forces closes the gap. To close it, something must touch the weak")
pr("  or strong force: new weak pairs or coloured states between the Z and the meeting (the 4x4 alternative's right-handed force did this, an earlier test;")
pr("  so do heavier partners of every particle in supersymmetry). The main line has none (rule 6).")
head("4. Summary")
pr(f"  FAIL: the matter-number force widens the gap ({100*(s/s2_meas - 1):+.1f}% -> about {100*(predict((6.3e10, 0.4))[0]/s2_meas - 1):+.1f}%). Shown in general: no single force can close it.")
pr("  SHARPENED: the strengths can only meet if something between the Z and the meeting energy touches the weak or strong force. In the main line")
pr("  that would be the six coloured far-side parts (44.2) taken as real particles below the meeting energy, or a break in rule 6.")
