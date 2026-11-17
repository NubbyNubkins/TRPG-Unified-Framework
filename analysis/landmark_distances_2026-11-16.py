"""Every distance between every landmark on the geometry (height = log10 weight),
which midpoints land on another landmark, the 7 x 5 = 34.5 search, and the piano-key (equal-step) test.
All 'distance' = decades of weight (the height on the cube). Midpoint of two = sqrt(product)."""
import numpy as np, itertools, math
rng=np.random.default_rng(1)
G=1.0  # GeV
LM={ # weight in GeV
 'horizon':1.44e-42,'nu1':3.6e-13,'DE':2.24e-12,'nu2':8.6e-12,'nu3':5.0e-11,
 'electron':5.11e-4,'muon':0.10566,'proton':0.938,'tau':1.777,'W':80.37,'Z':91.19,'Higgs':125.1,'v':246.2,'top':172.6,
 'N1':7.3e5,'MMG':2.41e8,'N2':1.3e9,'N3':6.3e10,'D':4.56e13,'grain':8.63e18}
BUILT={ # relations that hold by construction in TRPG / standard physics (not evidence)
 frozenset(['nu1','N1','electron']):'seesaw nu1*N1=me^2', frozenset(['nu2','N2','muon']):'seesaw nu2*N2=mmu^2',
 frozenset(['nu3','N3','tau']):'seesaw nu3*N3=mtau^2', frozenset(['MMG','grain','D']):'grain mirror MMG*grain=D^2',
 frozenset(['horizon','grain','DE']):'surface DE^2 ~ grain*horizon (middle 3.5 meV)'}
L={k:math.log10(v) for k,v in LM.items()}
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)

P("== an earlier test: distances on the geometry (decades of weight) ==")
names=list(L)
P("\n1. The comparisons (ratio, decades):")
for a,b in [('D','N1'),('D','N2'),('D','N3'),('N1','MMG'),('N2','MMG'),('N3','MMG'),('MMG','nu1'),('MMG','nu2'),('MMG','nu3'),('N2','N1'),('N3','N2'),('N3','N1')]:
    P(f"  {a:>5}/{b:<5} = {LM[a]/LM[b]:.3g}   ({L[a]-L[b]:.2f} decades)")
P("  MMG vs MMG: the three MMGs share one weight in the model (only the mean is fixed, an earlier test), so no gap is measurable yet.")

# midpoints landing on landmarks
P("\n2. Which midpoints sit on another landmark (within 0.30 decade = factor 2):")
hits=[]
for a,b in itertools.combinations(names,2):
    m=(L[a]+L[b])/2
    for c in names:
        if c in (a,b): continue
        d=m-L[c]
        if abs(d)<=0.30:
            tag=BUILT.get(frozenset([a,b,c]),'')
            hits.append((abs(d),a,b,c,d,tag))
for h in sorted(hits): P(f"  mid({h[1]},{h[2]}) ~ {h[3]}  off {h[4]:+.2f} dec  {('BUILT IN: '+h[5]) if h[5] else ''}")
nb=sum(1 for h in hits if not h[5])
# Three readings
def mid(a,b): return 10**((L[a]+L[b])/2)
P("\n   Three readings:")
P(f"   neutrinos <-> partners: mid(nu_i,N_i) = me, mmu, mtau exactly (seesaw, built in); overall mid of the two groups = {10**((np.mean([L['nu1'],L['nu2'],L['nu3']])+np.mean([L['N1'],L['N2'],L['N3']]))/2):.3g} GeV (Higgs 125: no)")
P(f"   mid(nu3,D) = {mid('nu3','D'):.3g} GeV ; mid(nu2,D) = {mid('nu2','D'):.3g} ; mid(nu1,N3) = {mid('nu1','N3'):.3g}")
P(f"   Higgs <-> grain: mid = {mid('Higgs','grain'):.3g} GeV (v: {mid('v','grain'):.3g}) vs N3 6.3e10 -> off {L['N3']-(L['Higgs']+L['grain'])/2:+.2f} / {L['N3']-(L['v']+L['grain'])/2:+.2f} dec")
P(f"   Higgs <-> D: mid = {mid('Higgs','D'):.3g} GeV (v: {mid('v','D'):.3g}) vs MMG 2.41e8 -> off {L['MMG']-(L['Higgs']+L['D'])/2:+.2f} / {L['MMG']-(L['v']+L['D'])/2:+.2f} dec")
P(f"   (with MMG*grain=D^2 built in, 'MMG at mid(v,D)' is the same as D^3 = v*grain^2: (v grain^2)^(1/3) = {(LM['v']*LM['grain']**2)**(1/3):.3g} vs D 4.56e13)")
# chance level
def count_nonbuilt(Ld):
    n=0
    for a,b in itertools.combinations(Ld,2):
        m=(Ld[a]+Ld[b])/2
        for c in Ld:
            if c in (a,b): continue
            if abs(m-Ld[c])<=0.30 and frozenset([a,b,c]) not in BUILT: n+=1
    return n
fixed=['horizon','nu1','nu2','nu3','electron','muon','proton','tau','W','Z','Higgs','v','top','grain','DE']
sim=[]
for _ in range(2000):
    Ld=dict(L)
    # re-throw the far-side landmarks (the ones whose places are TRPG's) log-uniformly between 1e4 and 1e17 GeV
    for k in ['N1','MMG','N2','N3','D']: Ld[k]=rng.uniform(4,17)
    sim.append(count_nonbuilt(Ld))
sim=np.array(sim)
P(f"\n   chance check: non-built-in midpoint hits = {nb}; with the 5 far-side landmarks thrown at random: mean {sim.mean():.1f}, {np.mean(sim>=nb)*100:.0f}% of random throws give as many or more")

# 3. 7 x 5 search
P("\n3. 'almost 7' x 'almost 5' = 34.5 (band 32.9-35.8) using the neutrino partners")
cands=[]
for a,b in itertools.permutations(names,2):
    if not ({a,b}&{'N1','N2','N3'}): continue
    r=LM[a]/LM[b]
    if r<=1: continue
    for kind,val in (('ratio',r),('decades',math.log10(r)),('e-folds',math.log(r)),('sqrt ratio',math.sqrt(r)),('cube root',r**(1/3))):
        cands.append((kind,a,b,val))
sev=[c for c in cands if 6.5<=c[3]<=7.5]; fiv=[c for c in cands if 4.5<=c[3]<=5.5]
P(f"   numbers near 7: {len(sev)}; near 5: {len(fiv)}")
for c in sev: P(f"     ~7: {c[0]} {c[1]}/{c[2]} = {c[3]:.3f}")
for c in fiv: P(f"     ~5: {c[0]} {c[1]}/{c[2]} = {c[3]:.3f}")
prods=[(s[3]*f[3],s,f) for s in sev for f in fiv if 32.92<=s[3]*f[3]<=35.75]
P(f"   products inside the 34.5 band: {len(prods)}")
for p,s,f in sorted(prods,key=lambda x:abs(x[0]-34.5)): P(f"     {s[0]} {s[1]}/{s[2]} ({s[3]:.3f}) x {f[0]} {f[1]}/{f[2]} ({f[3]:.3f}) = {p:.2f}")
def nprods(Ld):
    LMr={k:10**v for k,v in Ld.items()}; c=[]
    for a,b in itertools.permutations(Ld,2):
        if not ({a,b}&{'N1','N2','N3'}): continue
        r=LMr[a]/LMr[b]
        if r<=1: continue
        c+= [r,math.log10(r),math.log(r),math.sqrt(r),r**(1/3)]
    s=[x for x in c if 6.5<=x<=7.5]; f=[x for x in c if 4.5<=x<=5.5]
    return sum(1 for x in s for y in f if 32.92<=x*y<=35.75)
sim2=[]
for _ in range(1000):
    Ld=dict(L)
    for k in ['N1','N2','N3']: Ld[k]=rng.uniform(4,12)
    sim2.append(nprods(Ld))
sim2=np.array(sim2)
P(f"   chance check: with the three partners thrown at random (1e4-1e12 GeV): mean {sim2.mean():.1f} band-products, {np.mean(sim2>=len(prods))*100:.0f}% give as many or more")

# 4. piano keys / equal steps
P("\n4. Piano keys: do the landmarks sit on equal steps (a keyboard) better than random?")
P("   (an octave = 12 equal steps; 7 white + 5 black keys; 7 x 5 = 35)")
def lattice_score(vals,u):
    x=np.array(vals)/u; return np.sqrt(np.mean((x-np.round(x))**2))*u  # rms miss in decades
core=['nu1','nu2','nu3','electron','muon','tau','proton','Higgs','top','N1','MMG','N2','N3','D','grain']
vals=np.array([L[k] for k in core]); vals=vals-L['grain']
us=np.linspace(0.3,6,2000)
sc=[lattice_score(vals,u)/u for u in us]  # rms miss as fraction of a step
best=np.argmin(sc)
P(f"   best equal step: {us[best]:.3f} decades, rms miss {sc[best]*100:.1f}% of a step (random positions give ~29%)")
rb=[]
for _ in range(1000):
    r=rng.uniform(vals.min(),0,len(vals)); r[-1]=0
    rb.append(min(lattice_score(r,u)/u for u in us))
rb=np.array(rb); P(f"   random sets of 15 points, best step searched the same way: {np.mean(rb<=sc[best])*100:.0f}% do as well or better")
P("   white/black pattern (2,2,1,2,2,2,1 half-steps) on consecutive far-side gaps:")
far=['Higgs','N1','MMG','N2','N3','D','grain']; gaps=np.diff([L[k] for k in far])
P("   gaps (decades): "+', '.join(f'{far[i]}->{far[i+1]} {g:.2f}' for i,g in enumerate(gaps)))
open('landmark_distances_2026-11-16.txt','w').write('\n'.join(out)+'\n')
