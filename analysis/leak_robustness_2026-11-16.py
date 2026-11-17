"""Does the earlier leak result survive changes of scale and footing? Exploratory, not official.
Uses slib (two-loop running of force_strengths_final_2026-11-15.py) unchanged. Run: python3 leak_robustness_2026-11-16.py"""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
SIG=0.00034
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

def solve(c=1.0, MP=MPL, **kw):
    def meet(ML):
        s_=predict_s2(thr=thrL(ML),**kw); return mismatch(s2=s_,thr=thrL(ML),**kw)[1]
    f=lambda x: math.log10(c*meet(10**x)**2/MP)-x
    M=10**brentq(f,4,13,xtol=1e-6); s=predict_s2(thr=thrL(M),**kw); return M,meet(M),s
rows=[("baseline",{}),("grain lends half (c=0.5)",dict(c=0.5)),("grain lends double (c=2)",dict(c=2.0)),
("reduced Planck (grain / sqrt(8 pi))",dict(MP=MRED)),("one-loop running only",dict(two=False)),
("strong strength +1 sigma",dict(a3=A3+0.0009)),("strong strength -1 sigma",dict(a3=A3-0.0009)),
("top mass +0.3 GeV",dict(mt=MT+0.3)),("top mass -0.3 GeV",dict(mt=MT-0.3))]
for lab,kw in rows:
    try:
        M,F,s=solve(**kw); print(f"{lab:38s} pairs {M:.3e}  F {F:.3e}  sin2 {s:.5f}  {100*(s/S2-1):+.2f}%  {(s-S2)/SIG:+.1f} sigma")
    except Exception as e: print(lab,"ERR",e)
