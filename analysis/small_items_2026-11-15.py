"""TRPG small items (audit follow-up), 15 Nov 2026.
1 Hubble tension from 'thicker and thinner space' (local void)
2 Like with like
3 Right-angled unitarity triangle (66.3)
4 Middle fractions shrinking by sqrt(3/2) (61.2 P3)
5 Cone slices (56.1): one common tilted 45-degree cone for both quark ladders
6 Strong-CP angle: Nelson-Barr form of a phased fold (demo)
"""
import numpy as np
from fractions import Fraction as F
rng=np.random.default_rng(20261115)
np.set_printoptions(precision=4,suppress=True)
print("="*78); print("1. HUBBLE TENSION FROM A LOCAL UNDER-DENSE REGION ('thinner space')"); print("="*78)
H_late, sH_late = 73.04, 1.04     # SH0ES (Riess et al. 2022)
H_cmb,  sH_cmb  = 67.4, 0.5       # Planck 2018
Om, h, ns, s8 = 0.315, 0.674, 0.965, 0.811
need = H_late/H_cmb-1
f = Om**0.55
delta_need = -3*need/f
print(f"needed local expansion excess dH/H = {need:.4f}; growth rate f = Om^0.55 = {f:.3f}")
print(f"linear theory dH/H = -f*delta/3  ->  needed density contrast delta = {delta_need:.3f}")
print("observed: KBC void (Keenan, Barger & Cowie 2013) delta = -0.46 +/- 0.06 out to ~300 Mpc")
# LCDM rms fluctuation in top-hat spheres (BBKS transfer, Bardeen et al. 1986)
def T(k):
    q=k/(Om*h)            # k in h/Mpc
    return np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
k=np.logspace(-5,2,40000)
def W(x): return 3*(np.sin(x)-x*np.cos(x))/x**3
P0=k**ns*T(k)**2
def sig(R,A=1.0): return np.sqrt(A*np.trapezoid(P0*W(k*R)**2*k**2,k)/(2*np.pi**2))
A=(s8/sig(8.0))**2
for Rmpc in (300,640):
    R=Rmpc*h; s=sig(R,A)
    from math import erfc,sqrt
    p=0.5*erfc(abs(delta_need)/s/sqrt(2))
    print(f"R = {Rmpc} Mpc: LCDM rms delta = {s:.4f}; a void of {delta_need:.2f} is {abs(delta_need)/s:.0f} sigma (Gaussian chance {p:.0e})")
print("with the survey's own error (-0.46 +/- 0.06) the gap to anything LCDM makes is (0.46-0.05)/0.06 ~ 7 sigma;\n"
      "Haslbauer, Banik & Kroupa (2020) find 6.0 sigma. Other surveys find a shallower void (delta ~ -0.1 to -0.2).")
print("SN Hubble-diagram check: a void this deep and this large leaves a step in the SN\n"
      "distances near its edge; none is seen (Kenworthy, Scolnic & Riess 2019; Camarena et al. 2022).")
print("Verdict: the void he describes EXISTS (KBC, delta ~ -0.46, matching the needed -0.47),\n"
      "but ordinary growth cannot make it (>~10 sigma), and the SN data show no edge. It works only\n"
      "if large-scale structure grows much faster than in LCDM (Haslbauer, Banik & Kroupa 2020, MOND),\n"
      "which an earlier test found fails for the CMB and clusters. Main line: FAIL as the cause of the tension.")

print(); print("="*78); print("2. LIKE WITH LIKE"); print("="*78)
# The five 'directions' a single object can carry: 3 colours (Y=-1/3) and 2 weak (Y=+1/2)
dirs=[('colour',F(-1,3))]*3+[('weak',F(1,2))]*2
from itertools import combinations
pairs={}
for i,j in combinations(range(5),2):
    a,b=dirs[i][0],dirs[j][0]; key=('like '+a+'-'+b) if a==b else 'unlike colour-weak'
    pairs.setdefault(key,[]).append(dirs[i][1]+dirs[j][1])
names={'like colour-colour':'up partner (anti-colour, weak single)','like weak-weak':'electron partner (colourless, weak single)','unlike colour-weak':'quark pair (colour x weak pair)'}
for kk,v in pairs.items():
    print(f"{kk:22s}: {len(v)} states, Y = {sorted(set(v))}, charge Q=T3+Y -> {names[kk]}")
print("singles (one direction, conjugated): down partner Y=+1/3 (3 states), lepton pair Y=-1/2 (2 states)")
print("measured (left-handed conjugate fields): u^c Y=-2/3 (3), e^c Y=+1 (1), Q Y=+1/6 (6), d^c Y=+1/3 (3), L Y=-1/2 (2)")
print("Result: like-with-like antisymmetric pairs give EXACTLY the two partners the two-way reading lacked;\n"
      "unlike pairs give the quark pair; singles give the down partner and lepton pair. Total 10+5+1 = 16.\n"
      "This is the Georgi-Glashow (1974) grouping 16 = 10 + 5bar + 1. Exact identity, no fitted number.\n"
      "Antisymmetric = Pauli's rule for identical matter states. Does NOT bring SU(5) carriers back\n"
      "(main line has no quark-lepton carriers); it is a second reading of the same 16 states.")
# anomaly check
fs=[(F(-2,3),3),(F(1),1),(F(1,6),6),(F(1,3),3),(F(-1,2),2),(F(0),1)]
print("anomaly sums over the 16: sum Y =",sum(y*n for y,n in fs)," sum Y^3 =",sum(y**3*n for y,n in fs))

print(); print("="*78); print("3. RIGHT-ANGLED UNITARITY TRIANGLE (alpha = 90 deg)"); print("="*78)
alpha_fit, s_af = 91.6, 1.5        # CKM global fit (PDG 2024 / CKMfitter)
alpha_dir, s_ad = 85.2, 4.6        # direct measurement (PDG 2024)
sin2b, s_s2b = 0.709, 0.011        # HFLAV
beta=np.degrees(np.arcsin(sin2b))/2; sbeta=np.degrees(s_s2b/np.sqrt(1-sin2b**2))/2
gamma_dir, s_g = 65.9, 3.0         # PDG 2024 direct average (LHCb 64.6 +/- 2.8)
print(f"alpha fit {alpha_fit} +/- {s_af}: {(alpha_fit-90)/s_af:+.1f} sigma from 90;  direct {alpha_dir} +/- {s_ad}: {(alpha_dir-90)/s_ad:+.1f} sigma")
gpred=90-beta
print(f"beta from sin2beta = {beta:.2f} +/- {sbeta:.2f} deg -> right angle predicts gamma = {gpred:.2f} +/- {sbeta:.2f}")
print(f"measured gamma {gamma_dir} +/- {s_g}: {(gamma_dir-gpred)/s_g:+.1f} sigma. LHCb+Belle II aim for ~1 deg by ~2030: decisive.")
# chance: window allowed (+-2 sigma of fit) / plausible range, times number of 'special' angles
win=4*s_af; specials=[30,45,60,90,120]
print(f"luck: a 2-sigma window ({win:.0f} deg) over a 0-180 range is {win/180:.1%}; with {len(specials)} 'special' angles to choose from ~{len(specials)*win/180:.0%}")
print("Link to TRPG: in the Fritzsch-Xing quark parametrisation the one phase phi ~ alpha; phi = 90 deg is\n"
      "the 'maximal' phase, a quarter turn: the same quarter turn (rule 4's one i) that an earlier test\n"
      "found gives the main tilt 0.228 vs 0.2250 (Fritzsch & Xing 1999). So 'right-angled triangle' and\n"
      "'light corners a quarter turn apart' are one statement. Reading only.\n"
      "Verdict: ALLOWED (inside 1.1 sigma), not derived; LEAD with a sharp test: gamma = 67.4 deg.")

print(); print("="*78); print("4. MIDDLE FRACTIONS SHRINKING BY sqrt(3/2) (61.2 P3)"); print("="*78)
def frac(m1,m2,m3): return np.log(m2/m1)/np.log(m3/m1)
def mc(fn,vals,n=200000):
    s=[np.exp(rng.normal(np.log(v),e/v,n)) for v,e in vals]; x=fn(*s); return x.mean(),x.std()
footA={'lep':[(0.51099895,1e-8),(105.6583755,1e-6),(1776.86,0.12)],
       'up':[(2.16,0.07),(1273.,4.6),(172570.,290.)],
       'down':[(4.70,0.07),(93.5,0.8),(4183.,7.)]}
# all at M_Z (Huang & Zhou 2021)
footB={'lep':[(0.48307,1e-5),(101.766,1e-3),(1728.56,0.1)],
       'up':[(1.23,0.21),(620.,17.),(168260.,750.)],
       'down':[(2.67,0.19),(53.16,4.61),(2839.,26.)]}
r=np.sqrt(1.5)
for nm,ft in (('A (usual masses)',footA),('B (all at the Z)',footB)):
    fl,_=mc(frac,ft['lep']); fu,su=mc(frac,ft['up']); fd,sd=mc(frac,ft['down'])
    pu,pd=fl/r,fl/r**2
    print(f"footing {nm}: lepton {fl:.4f} | up {fu:.4f}+/-{su:.4f} (pred {pu:.4f}, {(fu-pu)/su:+.1f} sigma) | down {fd:.4f}+/-{sd:.4f} (pred {pd:.4f}, {(fd-pd)/sd:+.1f} sigma)")
# look-elsewhere: how many simple factors fit footing B as well
cands={'sqrt(3/2)':np.sqrt(1.5),'6/5':1.2,'5/4':1.25,'2^(1/4)':2**0.25,'golden^(1/2)':np.sqrt((1+5**.5)/2),'sqrt(4/3)':np.sqrt(4/3),'cbrt(2)':2**(1/3)}
fl,_=mc(frac,footB['lep']); fu,su=mc(frac,footB['up']); fd,sd=mc(frac,footB['down'])
ok=[k for k,c in cands.items() if abs(fu-fl/c)<2*su and abs(fd-fl/c**2)<2*sd]
print("simple factors that pass footing B within 2 sigma:",ok)
print("Verdict: footing A FAILS (up quark 24 sigma, down 2.5 sigma). Footing B passes within ~1 sigma but the quark\n"
      "masses at the Z are loose, so several factors pass: not evidence. A rule must name its footing\n"
      "; with the common-energy footing it stays ALLOWED / LEAD at most.")

print(); print("="*78); print("5. CONE SLICES: both quark ladders on ONE 45-degree cone with a tilted axis"); print("="*78)
def unit(v): return v/np.linalg.norm(v)
def common_axes(v1,v2,ang=45.0):
    """unit n with n.v1 = n.v2 = cos(ang) (v1,v2 unit). Returns list of solutions."""
    c=np.cos(np.radians(ang)); M=np.vstack([v1,v2])
    # n = a + t*w, with M a = c*1, w perpendicular to v1,v2
    a=np.linalg.lstsq(M,np.array([c,c]),rcond=None)[0]; w=unit(np.cross(v1,v2))
    disc=1-a@a
    if disc<0: return []
    return [a+np.sqrt(disc)*w, a-np.sqrt(disc)*w]
cube=[unit(np.array(d,float)) for d in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,-1,0),(1,0,-1),(0,1,-1),(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]]
def best_cube(n): return min(np.degrees(np.arccos(min(1,abs(n@d)))) for d in cube)
for nm,ft in (('A',footA),('B',footB)):
    u=unit(np.sqrt([x[0] for x in ft['up']])); d=unit(np.sqrt([x[0] for x in ft['down']]))
    ax=common_axes(u,d)
    one=unit(np.ones(3))
    print(f"footing {nm}: up opening {np.degrees(np.arccos(u@one)):.2f}, down {np.degrees(np.arccos(d@one)):.2f} deg")
    for n in ax:
        print(f"   common axis {n}: tilt from (1,1,1) {np.degrees(np.arccos(abs(n@one))):.2f} deg, nearest cube direction {best_cube(n):.2f} deg")
# control: random hierarchical triplets (log-uniform), how often a common axis exists within X deg of a cube dir
X=None
u=unit(np.sqrt([x[0] for x in footB['up']])); d=unit(np.sqrt([x[0] for x in footB['down']]))
real=min(best_cube(n) for n in common_axes(u,d))
cnt=0;tot=0;exists=0
for _ in range(20000):
    a=np.sort(np.exp(rng.uniform(0,np.log(2e5),3))); b=np.sort(np.exp(rng.uniform(0,np.log(2e3),3)))
    ax=common_axes(unit(np.sqrt(a)),unit(np.sqrt(b)))
    tot+=1
    if ax: exists+=1; cnt+= min(best_cube(n) for n in ax)<=real
print(f"control (random ladders, footing B): a common 45-deg axis exists {exists/tot:.0%} of the time;\n"
      f"one as close to a cube direction as the real one ({real:.2f} deg) {cnt/tot:.0%} of the time")
# is the 45-degree tilt of the common axis from (1,1,1) special?
tilt_real=None
for n in common_axes(u,d):
    t=np.degrees(np.arccos(abs(n@unit(np.ones(3)))))
    if abs(t-45)<abs((tilt_real or 999)-45): tilt_real=t
hit=0;tot=0
for _ in range(20000):
    a=np.sort(np.exp(rng.uniform(0,np.log(2e5),3))); b=np.sort(np.exp(rng.uniform(0,np.log(2e3),3)))
    ax=common_axes(unit(np.sqrt(a)),unit(np.sqrt(b)))
    if not ax: continue
    tot+=1; ts=[np.degrees(np.arccos(abs(n@unit(np.ones(3))))) for n in ax]
    hit+= min(abs(t-45) for t in ts)<=abs(tilt_real-45)
print(f"footing B: the common axis sits {tilt_real:.2f} deg from the lepton axis (45 would be 'a cone on a cone');")
print(f"random ladders do as well {hit/tot:.1%} of the time; footing A gives 47.97 deg instead")
print("Verdict: a slice rule CAN be stated (one cone, one tilted axis), but it always has solutions\n"
      "and the axis is not special (cube direction 51% chance; 45-deg tilt a few % and footing-dependent).\n"
      "FAIL as a rule; the 45-deg tilt in footing B is recorded as a coincidence to watch.")

print(); print("="*78); print("6. STRONG-CP ANGLE: phased fold in Nelson-Barr form (demo)"); print("="*78)
md=np.diag([2.67,53.16,2839.])*1e-3     # GeV, real, at the Z (d_L x d_R)
Mh=np.diag([2.0e5,3.0e5])                # two heavy D partners (E6 has one per family), real
Fm=Mh[0,0]*np.array([[0.9,0.5,1.2*np.exp(1j*1.0)],[0.4,1.1*np.exp(-1j*0.7),0.8]])  # complex links D_L-d_R
full=np.zeros((5,5),complex); full[:3,:3]=md; full[3:,:3]=Fm; full[3:,3:]=Mh  # d_L-D_R block ZERO
print(f"arg det(full down mass matrix) = {np.angle(np.linalg.det(full)):.1e} rad  (tree-level theta-bar)")
def light_left(full):
    H=full@full.conj().T; w,V=np.linalg.eigh(H); Vd=V[:3,:3]
    U,_,Wh=np.linalg.svd(Vd); return np.sqrt(np.abs(w[:3])), U@Wh
m,V=light_left(full)
J=np.imag(V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1]))
print(f"light down masses {np.round(m*1e3,3)} MeV; |V12| {abs(V[0,1]):.3f} |V23| {abs(V[1,2]):.4f} |V13| {abs(V[0,2]):.4f}; J = {J:.1e}")
print("(illustration: link values not fitted; the point is J != 0 while arg det = 0)")
m1=np.zeros((4,4),complex); m1[:3,:3]=md; m1[3,:3]=Fm[0]; m1[3,3]=Mh[0,0]
_,V1=light_left(m1); J1=np.imag(V1[0,1]*V1[1,2]*np.conj(V1[0,2])*np.conj(V1[1,1]))
print(f"one heavy D on a base already diagonal in the up basis: J = {J1:.1e} (no phase). With a non-diagonal\n"
      "real base (Bento, Branco & Parada 1991; the earlier base) one D is enough.")
full2=full.copy(); full2[:3,3:]=0.05*Mh[0,0]*np.exp(1j*0.4)
print(f"with the other block switched on: arg det = {np.angle(np.linalg.det(full2)):.3f} rad -> theta-bar not protected")
print("Verdict: a phased fold through a heavy down partner (E6's D) gives a CKM phase with theta-bar = 0\n"
      "at tree level IF one of its two link blocks vanishes\n"
      "(Nelson 1984; Barr 1984; Bento, Branco & Parada 1991). That is the\n"
      "standard Nelson-Barr solution of the strong-CP problem, and it uses the same part as item 2.\n"
      "No TRPG rule yet sets the zero block: LEAD. Otherwise theta-bar stays a measured input (<~1e-10).")
