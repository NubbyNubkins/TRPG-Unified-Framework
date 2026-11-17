"""The 'resting place' rule (rule 7 generalised).
Natural extrema of the invariants of the mass matrices (Michel & Radicati 1971; Cabibbo & Maiani 1970;
Alonso, Gavela, Isidori & Maiani 2013). Tests:
 A. numerically minimise generic invariant potentials for two 3x3 ladders and classify the resting points;
 B. geometry: the opening = angle between (sqrt m1, sqrt m2, sqrt m3) and the cube's body diagonal, so the
    resting points are the cube's own directions (edge / face diagonal / body diagonal);
 C. how far each real family sits from its resting point, both footings;
 D. chance: random (anarchic) matrices as the comparison;
 E. does the first step away from rest give the openings?"""
import numpy as np, math
from scipy.optimize import minimize
from scipy.linalg import expm
rng=np.random.default_rng(5)
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
def opening(m):
    r=np.sqrt(np.abs(np.asarray(m,float))); c=r.sum()/math.sqrt(3)/np.linalg.norm(r); return math.degrees(math.acos(min(1,c)))
# ---------- A
P("A. Resting points of generic invariant potentials (two ladders u,d; V = -T + lam T^2 + lam' Tr X^2 + g Tr(Xu Xd))")
def unit(p):
    H=np.zeros((3,3),complex); k=0
    for i in range(3):
        H[i,i]=p[k]; k+=1
        for j in range(i+1,3): H[i,j]=p[k]+1j*p[k+1]; H[j,i]=p[k]-1j*p[k+1]; k+=2
    return expm(1j*H)
def V(x,lp_u,lp_d,g):
    a2=x[:3]**2; b2=x[3:6]**2; U=unit(x[6:15]); M=np.abs(U)**2
    return -(a2.sum()+b2.sum())+0.5*(a2.sum()**2+b2.sum()**2)+lp_u*(a2**2).sum()+lp_d*(b2**2).sum()+g*(a2@M@b2)
tally={}
for trial in range(60):
    lp_u,lp_d=rng.choice([-0.3,0.3]),rng.choice([-0.3,0.3]); g=rng.choice([-0.2,0.2])
    best=None
    for s in range(6):
        r=minimize(V,rng.normal(size=15),args=(lp_u,lp_d,g),method='BFGS')
        if best is None or r.fun<best.fun: best=r
    a2=best.x[:3]**2; b2=best.x[3:6]**2; M=np.abs(unit(best.x[6:15]))**2
    def cls(v):
        v=np.sort(v)/v.max(); return 'one heavy' if v[1]<1e-3 else ('all equal' if v[0]>0.99 else ('two heavy' if v[0]<1e-3 and v[1]>0.99 else 'other'))
    # only overlaps between occupied directions matter
    iu=np.where(a2/a2.max()>1e-3)[0]; idd=np.where(b2/b2.max()>1e-3)[0]
    if cls(a2)=='all equal' or cls(b2)=='all equal': mix='mixing has no meaning (a ladder is all-equal)'
    else:
        ov=M[np.ix_(iu,idd)].sum()/min(len(iu),len(idd))   # 1 = fully aligned, 0 = fully apart
        mix='aligned (no mixing)' if ov>0.999 else ('fully apart (no mixing)' if ov<1e-3 else 'IN-BETWEEN ANGLE')
    key=(f"lam'_u {'<0' if lp_u<0 else '>0'}",f"lam'_d {'<0' if lp_d<0 else '>0'}",f"g {'<0' if g<0 else '>0'}",cls(a2),cls(b2),mix); tally[key]=tally.get(key,0)+1
for k,v in sorted(tally.items()): P("  ",v,"x",k)
nin=sum(v for k,v in tally.items() if k[-1]=='IN-BETWEEN ANGLE'); noth=sum(v for k,v in tally.items() if 'other' in k)
P(f"  -> ladders rest only at 'one heavy', 'two heavy' or 'all equal' (other: {noth}); resting points at an in-between mixing angle: {nin} of {sum(tally.values())}.")
# ---------- B
P("\nB. Geometry (exact): opening = angle between (sqrt m) and the cube's body diagonal")
for nm,v in [("one heavy (0,0,1) = cube EDGE",[0,0,1]),("two heavy (0,1,1) = FACE diagonal",[0,1,1]),("all equal (1,1,1) = BODY diagonal",[1,1,1])]:
    P(f"   {nm}: {opening(v):.4f} deg")
P(f"   edge + face diagonal = {opening([0,0,1])+opening([0,1,1]):.4f} deg, so their exact midpoint is 45.0000: the lepton opening (rule 7) is exactly halfway between 'one heavy' and 'two heavy'.")
# ---------- C
P("\nC. Where the real families sit (edge 54.74 = one heavy; 45 = midpoint; face diagonal 35.26)")
fam={'up, own scale':[2.16,1273,172570],'up, at M_Z':[1.23,620,168260],'down, own scale':[4.70,93.5,4183],'down, at M_Z':[2.67,53.2,2839],
     'charged leptons':[0.51099895,105.658,1776.93],'neutrinos (TRPG m1 0.36 meV)':[3.6e-4,8.62e-3,5.02e-2],'neutrinos (m1 = 0)':[0,8.61e-3,5.01e-2]}
for nm,m in fam.items():
    th=opening(m); P(f"   {nm:<30} opening {th:6.2f}  heaviest share {max(m)/sum(m)*100:5.1f}%  from edge {54.7356-th:5.2f} deg  fraction of the way 45->edge {(th-45)/9.7356:+.2f}")
ck=[13.02,2.38,0.20]; pm=[33.4,49.0,8.6]
P(f"   quark mixing angles {ck} deg (resting point: 0, 0, 0);  lepton mixing {pm} deg (Majorana resting point: one angle 45, with two equal neutrino weights)")
# ---------- D
P("\nD. Chance: random complex matrices (no rule at all), 200,000 pairs")
n=200000; hits=0; sh=[]; mx=[]
for _ in range(n//1000):
    A=rng.normal(size=(1000,3,3))+1j*rng.normal(size=(1000,3,3)); B=rng.normal(size=(1000,3,3))+1j*rng.normal(size=(1000,3,3))
    ua,sa,_=np.linalg.svd(A); ub,sb,_=np.linalg.svd(B)
    fa=sa[:,0]/sa.sum(1); fb=sb[:,0]/sb.sum(1)
    Vm=np.abs(np.einsum('nij,nik->njk',ua.conj(),ub))  # |V_ij|
    s13=Vm[:,0,2]; ang13=np.degrees(np.arcsin(np.clip(s13,0,1))); 
    c13=np.sqrt(1-np.clip(s13,0,1)**2); a12=np.degrees(np.arcsin(np.clip(Vm[:,0,1]/c13,0,1))); a23=np.degrees(np.arcsin(np.clip(Vm[:,1,2]/c13,0,1)))
    # heaviest paired with heaviest ordering (svd sorts descending): index 0 = heaviest
    big=np.maximum.reduce([np.minimum(a12,90-a12),np.minimum(a23,90-a23),np.minimum(ang13,90-ang13)])
    hits+=np.sum((fa>=0.977)&(fb>=0.977)&(big<=13.02)); sh+=list(fa[:50]); mx+=list(big[:50])
P(f"   random pairs with both ladders >= 97.7% in the heaviest AND all mixing <= 13 deg: {hits} of {n}  ({hits/n*100:.4f}%)")
P(f"   typical random ladder: heaviest share {np.median(sh)*100:.0f}%, largest mixing angle {np.median(mx):.0f} deg")
# ---------- E
P("\nE. First step away from rest: opening deficit vs the second rung")
for nm in ['up, own scale','up, at M_Z','down, own scale','down, at M_Z','charged leptons']:
    m=fam[nm]; b=math.sqrt(m[1]/m[2]); a=math.sqrt(m[0]/m[2]); th=opening(m)
    # first-order: tan(phi) = sqrt(3Q-1), Q ~ (1+b^2+a^2)/(1+a+b)^2
    P(f"   {nm:<18} sqrt(m2/m3) = {b:.4f}, sqrt(m1/m3) = {a:.5f} -> deficit from edge {54.7356-th:.2f} deg (first-order estimate {math.degrees((a+b)*3/(2*math.sqrt(2))):.2f})")
P("   The deficit is just the size of the lower rungs (about 61 deg x their square-root share). The rule gives WHERE rest is, not how far each family stepped off it.")
open('resting_points_2026-11-16.txt','w').write('\n'.join(out)+'\n')
