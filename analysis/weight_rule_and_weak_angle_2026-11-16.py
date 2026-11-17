"""The 'only colour weighs by itself', the hidden neutrino partner as an overlapping wave (with
possibly hidden colour), and solving for the weak angle directly. the earlier gravity-leak test is left exactly as it was.
Re-uses the two-loop running of force_strengths_final_2026-11-15.py unchanged.
Credit: Minkowski 1977; Yanagida 1979; Gell-Mann, Ramond & Slansky 1979; Mohapatra & Senjanovic 1980 (two-wave seesaw);
Chikashige, Mohapatra & Peccei 1981 (partner weight borrowed from the matter-number piece); Pati & Salam 1974 (lepton as fourth colour);
Georgi, Quinn & Weinberg 1974; Machacek & Vaughn 1983-85; Shaposhnikov & Wetterich 2010 (Higgs weight from the grain); NuFIT; DESI 2025.
Re-run: python3 weight_rule_and_weak_angle_2026-11-16.py > weight_rule_and_weak_angle_2026-11-16.txt"""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def pr(*a): print(*a, flush=True)
def head(s): pr("\n" + "="*110 + "\n" + s + "\n" + "="*110)
mZ=91.1876; AINV=127.951; A3=0.1180; S2=0.23122; MT=172.57; YT=0.9369
dA3=0.0009; dS2=0.00004; dAINV=0.009; dYT=0.005
MPL=1.220890e19; MRED=2.435e18
b_SM=np.array([41/10,-19/6,-7.0]); B_SM=np.array([[199/50,27/10,44/5],[9/10,35/6,12],[11/10,9/2,-26]]); cT=np.array([17/10,3/2,2.0])
b_top=np.array([17/30,1/2,2/3])                                                    # the top's share, removed below m_t
dbL=np.array([2/5,2/3,0]); dBL=np.array([[9/50,9/10,0],[3/10,49/6,0],[0,0,0]])       # one spin-1/2 weak pair (2,+-1/2), Dirac
dbD=np.array([4/15,0,2/3]); dBD=np.array([[4/75,0,16/15],[0,0,0],[2/15,0,38/3]])     # one D partner (3,-1/3), Dirac
K=16*math.pi**2
def run(s2=S2,a3=A3,ainv=AINV,thr=(),two=True,topdec=True,stop='23',tmax=60.0,yt=YT,mt=MT):
    """thr: list of (mass, how many, db, dB). stop = '23' (strong = weak), '12' (weak = hypercharge) or a mass in GeV."""
    a=1/ainv; g=[math.sqrt(4*math.pi*a/(1-s2)*5/3), math.sqrt(4*math.pi*a/s2), math.sqrt(4*math.pi*a3)]
    cuts=sorted(set([math.log(mt/mZ)]+[math.log(m/mZ) for m,_,_,_ in thr]))
    tend=tmax if isinstance(stop,str) else math.log(stop/mZ)
    edges=[0.0]+[c for c in cuts if 0<c<tend]+[tend]
    y=np.array(g+[0.0]); tmt=math.log(mt/mZ)
    for i in range(len(edges)-1):
        t0,t1=edges[i],edges[i+1]; tm=0.5*(t0+t1)
        top=(tm>tmt) or (not topdec)
        if top and y[3]==0.0: y[3]=yt
        bb=b_SM.copy(); BB=B_SM.copy()
        if not top: bb=bb-b_top
        for m,n,db,dB in thr:
            if tm>math.log(m/mZ): bb=bb+n*db; BB=BB+n*dB
        def rhs(t,y):
            g=y[:3]; Y=y[3]; g2=g*g
            dg=g**3*bb/K
            if two: dg=dg+g**3/K**2*(BB@g2-cT*Y*Y)
            dY=Y/K*(4.5*Y*Y-8*g2[2]-2.25*g2[1]-17/20*g2[0])
            return [dg[0],dg[1],dg[2],dY]
        ev=None
        if stop=='23': ev=lambda t,y:y[1]-y[2]
        if stop=='12': ev=lambda t,y:y[0]-y[1]
        if ev is not None: ev.terminal=True
        sol=solve_ivp(rhs,[t0,t1],y,events=ev,rtol=1e-10,atol=1e-13,max_step=0.5)
        if ev is not None and sol.t_events[0].size:
            return mZ*math.exp(sol.t_events[0][0]), sol.y_events[0][0]
        y=sol.y[:,-1]
    if isinstance(stop,str): return None
    return stop,y
ai=lambda g:4*math.pi/g**2
def mismatch(**kw):
    """At the energy where strong = weak: (1/alpha1 - 1/alpha2, that energy, 1/alpha there). Zero means all three meet."""
    r=run(stop='23',**kw)
    if r is None: return None
    M,g=r; return ai(g[0])-ai(g[1]), M, ai(g[1])
def predict_s2(**kw):
    """Method A: inputs alpha and alpha_s; demand one meeting point; predict sin^2 at the Z."""
    return brentq(lambda s:mismatch(s2=s,**kw)[0],0.19,0.26,xtol=1e-9)
def thrL(ML,n=3): return [(ML,n,dbL,dBL)]
def thrLD(ML,MD,n=3): return [(ML,n,dbL,dBL),(MD,n,dbD,dBD)]
def ML_close(MD=None,n=3,lo=0.5,hi=13.5,**kw):
    f=lambda x: mismatch(thr=(thrL(10**x,n) if MD is None else thrLD(10**x,MD,n)),**kw)[0]
    hi2=hi if MD is None else min(hi,math.log10(MD))
    if f(lo)*f(hi2)>0: return None
    return 10**brentq(f,lo,hi2,xtol=1e-7)
pct=lambda s:100*(s/S2-1)

def s2at(mu, s2=S2, thr=()):
    r=run(s2=s2,thr=thr,stop=mu); g=r[1]; a=(3/5)*g[0]**2; return a/(a+g[1]**2)
print("SM sin2 (MS-bar) at footings:")
for nm,mu in (("top",172.57),("1 TeV",1e3),("far side 2.1e8",2.1e8),("meeting-ish 4.6e13",4.6e13),("grain 1.22e19",MPL)):
    print(f"  {nm:22s} {s2at(mu):.5f}")
# where does SM sin2 equal 1/4, 2/9, 3/8?
for tgt in (0.25,2/9,3/8):
    f=lambda x: s2at(10**x)-tgt
    xs=np.linspace(2.3,19,60); vals=[f(x) for x in xs]
    roots=[brentq(f,xs[i],xs[i+1]) for i in range(len(xs)-1) if vals[i]*vals[i+1]<0]
    print(f"  SM sin2 = {tgt:.4f} at: {[f'{10**r:.2e} GeV' for r in roots]}")
# rule: sin2 = 3/8 at the grain -> predict sin2(MZ) from alpha_em only
for lab,thr in (("SM only",()),("SM + far-side pairs at 1.52e8 (A-ii), D above",thrL(1.515e8)),("SM + pairs 1.52e8 + D at 1.6e15",thrLD(1.515e8,1.6e15))):
    for MU in (MPL,MRED):
        try:
            x=brentq(lambda s: s2at(MU,s2=s,thr=thr)-3/8,0.15,0.30,xtol=1e-9); print(f"  3/8 at {MU:.2e} [{lab}]: predicted sin2(MZ) = {x:.5f} ({100*(x/S2-1):+.2f}%, {(x-S2)/0.00004:+.0f} sigma exp)")
        except Exception as e: print("  fail",lab,MU,e)
# on-shell
mW,dmW,mZ_=80.3692,0.0133,91.1876; so=1-(mW/mZ_)**2; dso=2*mW*dmW/mZ_**2
print(f"on-shell sin2 = {so:.5f} +- {dso:.5f}; 2/9 = {2/9:.5f}: {(2/9-so)/dso:+.1f} sigma; 1/4: {(0.25-so)/dso:+.0f} sigma")
print("Angle = geometric 3/8 at a TRPG scale fixed beforehand (alpha_em only, strong not used):")
for nm,X in (("Higgs hold fades 10^10.5",10**10.5),("sqrt(v x grain)",math.sqrt(246.22*MPL)),("heaviest hidden partner = majoron scale 6.3e10",6.3e10),("right side breaks 10^11.7",10**11.7),("F of an earlier test (4.30e13)",4.30e13)):
    for lab,thr in (("SM",()),("SM+pairs",thrL(1.515e8))):
        x=brentq(lambda s: s2at(X,s2=s,thr=thr)-3/8,0.12,0.30,xtol=1e-9)
        print(f"  {nm:46s} [{lab:8s}]: sin2(MZ) = {x:.5f} ({100*(x/S2-1):+.2f}%)")
# neutrino partners as leak
me,mmu,mtau=0.51099895e-3,0.1056584,1.77686; MPL_=1.220890e19
mnu=np.array([0.36e-3,8.66e-3,50.34e-3])*1e-9  # GeV (TRPG)
N=np.array([me,mmu,mtau])**2/mnu; print("partners", N)
F=4.30e13; L=F**2/MPL_; print(f"leak F^2/M_P = {L:.3e}; geometric mean of partners {np.prod(N)**(1/3):.3e}; ratio {np.prod(N)**(1/3)/L:.2f}")
print("family-blind leak would need m_nu ratios = (lepton mass ratios)^2:", (mtau/mmu)**2, "vs TRPG/measured m3/m2", mnu[2]/mnu[1])
d21,d31=7.42e-5,2.517e-3   # eV^2 (NuFIT; TRPG's own masses give the same splittings)
d21=(8.66e-3**2-0.36e-3**2); d31=(50.34e-3**2-0.36e-3**2)
def gm(m1):
    m=np.array([m1,math.sqrt(m1**2+d21),math.sqrt(m1**2+d31)])*1e-9; return np.prod(np.array([me,mmu,mtau])**2/m)**(1/3)
m1=brentq(lambda x: gm(x)-L,1e-6,0.2); m=np.array([m1,math.sqrt(m1**2+d21),math.sqrt(m1**2+d31)])
print(f"leak sets the overall scale -> lightest neutrino {1e3*m1:.2f} meV, masses {np.round(1e3*m,2)} meV, sum {m.sum():.4f} eV (TRPG rule 5: 0.36 meV, sum 0.0594; DESI 2025: sum < 0.064 at 95%)")
for F_ in (5.21e13,):
    L2=F_**2/MPL_; m1=brentq(lambda x: gm(x)-L2,1e-6,0.2); m=np.array([m1,math.sqrt(m1**2+d21),math.sqrt(m1**2+d31)]); print(f"  with exact-closure F {F_:.2e}: lightest {1e3*m1:.2f} meV, sum {m.sum():.4f} eV")
# active-sterile mixing (camouflage)
th=np.sqrt(mnu/N); print("mixing angle light-heavy per family:", th)
