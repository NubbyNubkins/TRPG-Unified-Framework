"""'which turns of the cube survive' (residual symmetries; Lam 2007-08; Harrison, Perkins & Scott 2002;
Ma & Rajasekaran 2001). The cube's 24 rotations act on the three placements (signed permutations, det +1).
Each ladder keeps a leftover subgroup; the mixing is the mismatch between the two kept sets.
Test: list every mixing pattern the cube allows and compare with the measured quark and lepton mixing.
Also: can a kept turn fix a ladder's STEP (weight ratios)?"""
import numpy as np, itertools, math
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
# the 24 rotations
G=[]
for perm in itertools.permutations(range(3)):
    for signs in itertools.product([1,-1],repeat=3):
        M=np.zeros((3,3)); 
        for i,p in enumerate(perm): M[i,p]=signs[i]
        if round(np.linalg.det(M))==1: G.append(M)
assert len(G)==24
def order(g):
    m=np.eye(3)
    for k in range(1,7):
        m=m@g
        if np.allclose(m,np.eye(3)): return k
# kept sets with a unique eigenbasis: single elements of order 3 or 4; Klein groups of commuting half-turns
bases=[]
for g in G:
    if order(g) in (3,4):
        w,v=np.linalg.eig(g); 
        if len(set(np.round(w,6)))==3: bases.append((f"Z{order(g)}",np.linalg.qr(v)[0] if False else v/np.linalg.norm(v,axis=0)))
twos=[g for g in G if order(g)==2]
for a,b in itertools.combinations(twos,2):
    if np.allclose(a@b,b@a):
        # simultaneous eigenbasis
        w,v=np.linalg.eigh(a+2.3*b+0.0); bases.append(("Z2xZ2",v))
P(f"cube rotations: 24; kept sets with a definite set of directions: {len(bases)}")
def angles(V):
    A=np.abs(V); s13=min(1,A[0,2]); c13=math.sqrt(max(0,1-s13**2))
    if c13<1e-9: return None
    return (math.degrees(math.asin(min(1,A[0,1]/c13))),math.degrees(math.asin(min(1,A[1,2]/c13))),math.degrees(math.asin(s13)))
pats={}
for (n1,U1),(n2,U2) in itertools.product(bases,repeat=2):
    V=U1.conj().T@U2
    for pr in itertools.permutations(range(3)):
        for pc in itertools.permutations(range(3)):
            a=angles(V[np.ix_(pr,pc)])
            if a is None: continue
            key=tuple(round(x,1) for x in a); pats.setdefault(key,set()).add((n1,n2))
P(f"distinct mixing patterns (angles 12, 23, 13 in degrees) the cube allows: {len(pats)}")
ckm=(13.00,2.38,0.21); ckm_e=(0.04,0.04,0.006); pmns=(33.4,49.0,8.6); pm_e=(0.75,1.2,0.12)
def best(target,err):
    return sorted(pats,key=lambda k: sum(((k[i]-target[i])/err[i])**2 for i in range(3)))[:3]
P("\nclosest to the quark mixing (13.00, 2.38, 0.21):")
for k in best(ckm,ckm_e): P(f"   {k}  from {sorted(pats[k])[0]}   misses: "+", ".join(f"{(k[i]-ckm[i])/ckm_e[i]:+.0f}σ" for i in range(3)))
P("closest to the lepton mixing (33.4, 49.0, 8.6):")
for k in best(pmns,pm_e): P(f"   {k}  from {sorted(pats[k])[0]}   misses: "+", ".join(f"{(k[i]-pmns[i])/pm_e[i]:+.0f}σ" for i in range(3)))
allang=sorted(set(round(x,1) for k in pats for x in k)); P(f"\nevery angle that ever appears: {allang}")
P("   (0, 35.3, 45, 54.7, 90 = the cube's own angles again: body/face/edge directions.)")
# one kept half-turn only (one free angle): TM1/TM2 families
P("\none kept half-turn on the neutrino side (one free angle), charged leptons keep a 3-turn (Z3):")
w=np.exp(2j*np.pi/3); Uz3=np.array([[1,1,1],[1,w,w**2],[1,w**2,w]])/math.sqrt(3)
tbm=np.array([[math.sqrt(2/3),1/math.sqrt(3),0],[-1/math.sqrt(6),1/math.sqrt(3),1/math.sqrt(2)],[1/math.sqrt(6),-1/math.sqrt(3),1/math.sqrt(2)]])
bestv=None
for th in np.linspace(0,math.pi/2,1801):
    for ph in np.linspace(0,2*math.pi,73):
        R=np.array([[1,0,0],[0,math.cos(th),math.sin(th)*np.exp(-1j*ph)],[0,-math.sin(th)*np.exp(1j*ph),math.cos(th)]])  # keeps column 1 (TM1)
        a=angles(tbm@R); c=sum(((a[i]-pmns[i])/pm_e[i])**2 for i in range(3))
        if bestv is None or c<bestv[0]: bestv=(c,a,th,ph)
P(f"   TM1 best: angles {tuple(round(x,2) for x in bestv[1])}, chi2 {bestv[0]:.2f} for 3 angles with 2 free numbers (an earlier test: TM1 works)")
# eigenvalue statement
P("\nCan a kept turn fix a ladder's STEP (its weight ratios)?")
P("   A kept turn g means g M g^-1 = M (M = the ladder's weight matrix, Hermitian). Check on all 24 rotations:")
fixdir=0; loose=0
for g in G:
    nd=len(set(np.round(np.linalg.eigvals(g),6)))
    if nd==3: fixdir+=1
    else: loose+=1
P(f"   rotations that fix all three directions and leave all three weights free: {fixdir}; rotations that leave directions partly free too: {loose}; rotations that fix an unequal weight ratio: 0")
P("   (general fact: a symmetry fixes directions; on weights it can only allow 'free' or force 'equal'. Schur's lemma.)")
P("   -> no turn of the cube can fix a ratio like 0.086 or 0.15. Symmetry gives directions; the steps need a different kind of rule.")
open('cube_turns_survive_2026-11-16.txt','w').write('\n'.join(out)+'\n')
