"""Attempt to DERIVE two of TRPG's put-in rules from its own ideas (18 Nov 2026).
TESTS FIXED BEFORE RUNNING
A. The 45 degrees (Koide 1982; Brannen 2006 matrix form).
   Ideas used: (i) three families = three placements of one shape -> the family matrix is unchanged by a 120-degree turn
   (Z3), so the amplitude (square-root-of-mass) matrix is circulant: diagonal a = tie to itself, off-diagonal b = ties to the others.
   (ii) maximum uncertainty between self and relationship: the share p of the matrix's total weight (sum of squares) on the
   diagonal is the 'self' share; maximum uncertainty of the split {p, 1-p} for ANY symmetric measure (Shannon, Renyi, variance) is p = 1/2.
   PASS if p = 1/2 is shown to be exactly the 45 degrees (algebra), reproduces tau within 2 sigma, and does not contradict other data.
   If it holds for leptons only, the label is PARTIAL with the quark exception named.
B. The grain factor sqrt2.
   Idea used: containment (idea 2): the grain is where one particle is exactly over-contained by its own gravity (its own black hole).
   PASS only if TRPG's ideas single out ONE definition of that balance and it gives sqrt2. If several equally natural definitions give
   different factors, it is NOT derived (stays a chosen value with a reason).
"""
import math, numpy as np
print("A. THE 45 DEGREES")
print("  algebra: circulant amplitude matrix, eigenvalues sqrt(m_k) = a + 2|b| cos(delta + 2 pi k/3)")
print("  sum sqrt m = 3a ; sum m = 3a^2 + 6|b|^2 ; Koide Q = sum m/(sum sqrt m)^2 = 1/3 + 2|b|^2/(3a^2)")
print("  diagonal weight 3a^2, off-diagonal weight 6|b|^2  ->  p = 3a^2/(3a^2+6|b|^2)")
print("  p = 1/2  <=>  |b|^2 = a^2/2  <=>  Q = 2/3  <=>  angle between (sqrt m) and (1,1,1) = 45 degrees   [exact]")
def p_and_angle(m):
    s=np.sqrt(np.array(m,float)); a=s.mean(); b2=((s**2).sum()-3*a*a)/6
    Q=sum(m)/s.sum()**2; ang=math.degrees(math.acos(s.sum()/math.sqrt(3*(s**2).sum())))
    return 3*a*a/(3*a*a+6*b2), Q, ang
for lab,m in [("charged leptons (e, mu, tau) MeV",[0.51099895,105.6583755,1776.86]),
              ("down quarks (d, s, b) at 2 GeV / m_b(m_b)",[4.70,93.5,4183]),
              ("up quarks (u, c, t) at 2 GeV / m_c(m_c) / m_t pole",[2.16,1273,172570]),
              ("neutrinos (TRPG rule) meV",[0.36,8.66,50.13])]:
    p,Q,ang=p_and_angle(m); print(f"  {lab:48s} self share p = {p:.5f}  Q = {Q:.5f}  angle = {ang:.2f} deg")
# tau prediction from p = 1/2 with e and mu
def tau_from(p_target=0.5):
    e,mu=0.51099895,105.6583755
    f=lambda t: p_and_angle([e,mu,t])[0]-p_target
    lo,hi=1000,3000
    for _ in range(200):
        mid=(lo+hi)/2
        if f(lo)*f(mid)<=0: hi=mid
        else: lo=mid
    return (lo+hi)/2
t=tau_from(); print(f"  tau from p = 1/2: {t:.3f} MeV vs 1776.86 +- 0.12 -> {(t-1776.86)/0.12:+.1f} sigma")
print("  robustness: Shannon -p ln p-(1-p)ln(1-p), Renyi of any order, and p(1-p) all peak at p = 1/2 (symmetric in p <-> 1-p): yes")
print("  quarks: p != 1/2 (down 0.44, up 0.38 region) -> the rule as stated does not hold for quarks; colour is the named difference")
print()
print("B. THE GRAIN FACTOR")
cands=[("gravity radius 2GM/c^2 = reduced wavelength hbar/Mc", math.sqrt(2)),
       ("gravity radius 2GM/c^2 = Heisenberg minimum hbar/(2Mc)", 2.0),
       ("gravity radius 2GM/c^2 = full wavelength h/Mc", 1/math.sqrt(math.pi)),
       ("self-energy GM^2/R = quantum energy hbar c/R (energy balance)", 1.0),
       ("Einstein-equation coupling 8 pi G (reduced Planck)", math.sqrt(8*math.pi))]
for lab,c in cands: print(f"  {lab:62s} grain = Planck / {c:.3f}")
print("  -> five equally natural readings of 'over-contained by its own gravity' give five factors. TRPG's ideas do not single one out.")
print("  -> sqrt2 is NOT derived. It stays a chosen value with a physical reason (a particle becoming its own black hole); the 0.3% lattice")
print("     alpha_s decides it (sqrt2 -> 0.1178; 1 -> 0.1196; 2 -> ~0.1166 by the earlier run's slope).")
print()
print("A2. SIGNED AMPLITUDES (an amplitude matrix may have a negative eigenvalue; TRPG's neutrino rule already uses one)")
def p_signed(s):
    s=np.array(s,float); a=s.mean(); b2=((s**2).sum()-3*a*a)/6; return 3*a*a/(3*a*a+6*b2)
import itertools
for lab,m in [("neutrinos (TRPG rule) meV",[0.36,8.66,50.13]),("down quarks",[4.70,93.5,4183]),("up quarks",[2.16,1273,172570])]:
    r=np.sqrt(m); best=[]
    for sg in itertools.product([1,-1],repeat=3):
        if sg[2]<0: continue
        best.append((abs(p_signed(np.array(sg)*r)-0.5),sg,p_signed(np.array(sg)*r)))
    best.sort(); d,sg,p=best[0]
    print(f"  {lab:28s} best sign pattern {sg}: self share p = {p:.5f}")
print("  -> neutrinos sit exactly at p = 1/2 with the lightest amplitude negative (as TRPG's rule already states): same principle as the charged leptons.")
print("  -> quarks: no sign pattern reaches 1/2 within their own type. The exception stands.")
