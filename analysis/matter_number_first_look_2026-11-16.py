"""Does "the far side settles along matter number" (Dimopoulos & Wilczek 1981 form)
give the weak pairs their second-look weight F^2/grain with three placements? Three-loop running as in leak_cause_and_grain_2026-11-16.py."""
import math, numpy as np
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

run2=run
# --- three-loop gauge terms (Mihaila, Salomon & Steinhauser 2012, PRL 108, 151602 and PRD 86, 096008), n_G=3, n_t=1 ---
# beta_i = mu^2 d(alpha_i/pi)/dmu^2 ;  x_i = alpha_i/pi, xt = alpha_t/pi = y_t^2/(4 pi^2), xl = lambda/(4 pi^2)  (L = -lambda (Phi^+ Phi)^2)
nG=3; nt=1; LAM0=0.126
def C3(x1,x2,x3,xt,xl):
    c1=(x1*x1*(489/512000-29*nG/2400-209*nG**2/8640)+x1*x2*(783/51200-7*nG/6400)+x1*x3*(-137*nG/14400)
        +x2*x2*(3401/20480+83*nG/1920-11*nG**2/960)+x2*x3*(-nG/320)+x3*x3*(275*nG/576-121*nG**2/2160)
        +nt*xt*(x1*(-2827/51200)+x2*(-471/2048)+x3*(-29/320))+nt*xt*xt*(339/5120+303*nt/2560)
        +xl*(x1*27/3200+x2*9/640-xl*9/320))
    c2=(x1*x1*(163/102400-7*nG/960-11*nG**2/2880)+x1*x2*(561/10240+13*nG/1280)+x1*x3*(-nG/960)
        +x2*x2*(-667111/110592+1603*nG/432-415*nG**2/1728)+x2*x3*(13*nG/64)+x3*x3*(125*nG/192-11*nG**2/144)
        +nt*xt*(x1*(-593/10240)+x2*(-729/2048)+x3*(-7/64))+nt*xt*xt*(57/1024+45*nt/512)
        +xl*(x1*3/640+x2*3/128-xl*3/64))
    c3=(-5714*x3*x3-(101/10)*x1*xt-(93/2)*x2*xt-160*x3*xt+18*xt*xt+42*xt*xt
        +nG*(-(13/30)*x1*x1-(1/10)*x1*x2+(241/6)*x2*x2+(308/45)*x1*x3+28*x2*x3+(20132/9)*x3*x3)
        +nG**2*(-(242/135)*x1*x1-(22/3)*x2*x2-(2600/27)*x3*x3))/256
    return np.array([c1,c2,c3])
def run(s2=S2,a3=A3,ainv=AINV,thr=(),two=True,three=True,topdec=True,stop='23',tmax=60.0,yt=YT,mt=MT):
    a=1/ainv; g=[math.sqrt(4*math.pi*a/(1-s2)*5/3), math.sqrt(4*math.pi*a/s2), math.sqrt(4*math.pi*a3)]
    cuts=sorted(set([math.log(mt/mZ)]+[math.log(m/mZ) for m,_,_,_ in thr]))
    tend=tmax if isinstance(stop,str) else math.log(stop/mZ)
    edges=[0.0]+[c for c in cuts if 0<c<tend]+[tend]
    y=np.array(g+[0.0,0.0]); tmt=math.log(mt/mZ)
    for i in range(len(edges)-1):
        t0,t1=edges[i],edges[i+1]; tm=0.5*(t0+t1)
        top=(tm>tmt) or (not topdec)
        if top and y[3]==0.0: y[3]=yt; y[4]=LAM0
        bb=b_SM.copy(); BB=B_SM.copy()
        if not top: bb=bb-b_top
        for m,n,db,dB in thr:
            if tm>math.log(m/mZ): bb=bb+n*db; BB=BB+n*dB
        def rhs(t,y):
            g=y[:3]; Y=y[3]; L=y[4]; g2=g*g
            dg=g**3*bb/K
            if two: dg=dg+g**3/K**2*(BB@g2-cT*Y*Y)
            if three:
                x=g2/(4*math.pi**2); dg=dg+4*math.pi**2*x*x*C3(x[0],x[1],x[2],Y*Y/(4*math.pi**2),L/(4*math.pi**2))/g
            dY=Y/K*(4.5*Y*Y-8*g2[2]-2.25*g2[1]-17/20*g2[0])
            gp2=0.6*g2[0]
            dL=(24*L*L+12*L*Y*Y-6*Y**4-3*L*(3*g2[1]+gp2)+0.375*(2*g2[1]**2+(g2[1]+gp2)**2))/K if Y>0 else 0.0
            return [dg[0],dg[1],dg[2],dY,dL]
        ev=None
        if stop=='23': ev=lambda t,y:y[1]-y[2]
        if stop=='12': ev=lambda t,y:y[0]-y[1]
        if ev is not None: ev.terminal=True
        sol=solve_ivp(rhs,[t0,t1],y,events=ev,rtol=1e-10,atol=1e-13,max_step=0.5)
        if ev is not None and sol.t_events[0].size:
            return mZ*math.exp(sol.t_events[0][0]), sol.y_events[0][0][:4]
        y=sol.y[:,-1]
    if isinstance(stop,str): return None
    return stop,y[:4]
SIG=0.00034
GR=MPL/math.sqrt(2)   # grain where no-escape radius 2GM/c^2 = leakage length hbar/(Mc)
def meet(ML,n=3):
    s_=predict_s2(thr=thrL(ML,n)); return mismatch(s2=s_,thr=thrL(ML,n))[1]

# Matter-number rule (Dimopoulos & Wilczek 1981 form): the far side settles along matter number (B-L).
# Each placement's far-side 10 holds D (B-L != 0) and the weak pair L (B-L = 0).
# Fermionic 10s: coupling y_ij psi_i 45 psi_j must be ANTISYMMETRIC in families (45 antisymmetric, Majorana bilinear symmetric)
# -> 3x3 antisymmetric matrix has rank 2: two D's weigh F, ONE D combination weighs only at second order, F^2/grain, with the 3 pairs.
def thr_case(ML,nDlight):
    t=[(ML,3,dbL,dBL)]
    if nDlight: t.append((ML,nDlight,dbD,dBD))
    return t
def meet_c(ML,nDl):
    th=thr_case(ML,nDl); s_=predict_s2(thr=th); return mismatch(s2=s_,thr=th)[1], s_
for lab,Mg in (("full Planck",MPL),("grain 2 (M_Pl/sqrt2)",GR)):
    for nDl in (0,1):
        f=lambda x: math.log10(meet_c(10**x,nDl)[0]**2/Mg)-x
        xs=list(np.arange(3,13.6,0.25)); v=[f(x) for x in xs]
        roots=[brentq(f,xs[i],xs[i+1],xtol=1e-6) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
        for r in roots:
            F,s=meet_c(10**r,nDl)
            print(f"{lab:22s} light D's={nDl}: pairs {10**r:.3e} GeV, F {F:.3e}, sin2 {s:.5f}  {(s-S2)/SIG:+.1f} sigma ({pct(s):+.2f}%)")
        if not roots: print(f"{lab:22s} light D's={nDl}: no self-consistent solution")
# what grain would rescue the rank-2 case?
print("\nRank-2 case: sin2 vs grain factor c (grain = M_Pl/c)")
for c in (0.25,0.5,1,math.sqrt(2),2,4,8,16):
    Mg=MPL/c; f=lambda x: math.log10(meet_c(10**x,1)[0]**2/Mg)-x
    xs=list(np.arange(3,13.6,0.25)); v=[f(x) for x in xs]
    roots=[brentq(f,xs[i],xs[i+1],xtol=1e-6) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
    for r in roots:
        F,s=meet_c(10**r,1); print(f"  c={c:6.3f}: pairs {10**r:.3e} F {F:.3e} sin2 {s:.5f} {(s-S2)/SIG:+.1f} sigma")
    if not roots: print(f"  c={c:6.3f}: none")
# Rescue: the left-behind D is the press's heavy partner. Its own weight = F^2/grain (second look), the press link sets
# 1-kappa = own^2/(own^2+link^2) = 8.4e-4 (an earlier test, band 7.8-9.2e-4) -> heavy D mass = own/sqrt(1-kappa) ~ 34.5 x own.
print("\n--- rescue by the press link ---")
for omk in (7.8e-4,8.4e-4,9.2e-4):
  r_=1/math.sqrt(omk)
  for lab,Mg in (("full Planck",MPL),("grain 2",GR)):
    def meet_r(ML):
        th=[(ML,3,dbL,dBL),(ML*r_,1,dbD,dBD)]; s_=predict_s2(thr=th); return mismatch(s2=s_,thr=th)[1], s_
    f=lambda x: math.log10(meet_r(10**x)[0]**2/Mg)-x
    xs=list(np.arange(3,13.6,0.25)); v=[f(x) for x in xs]
    roots=[brentq(f,xs[i],xs[i+1],xtol=1e-6) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
    for r in roots:
        F,s=meet_r(10**r); print(f"1-k={omk:.1e} {lab:12s}: pairs {10**r:.3e}, press D {10**r*r_:.3e}, F {F:.3e}, sin2 {s:.5f} {(s-S2)/SIG:+.1f} sigma")
    if not roots: print(omk,lab,"none")
# how heavy must the left-behind D be for the strengths to stay within 2 sigma (grain 2, pairs fixed by self-consistency)?
print("\nleft-behind D mass needed (grain 2):")
for lg in (10,11,12,13,13.5,14,15,16):
    MD=10**lg
    def meet_m(ML):
        th=[(ML,3,dbL,dBL),(MD,1,dbD,dBD)]; s_=predict_s2(thr=th); return mismatch(s2=s_,thr=th)[1], s_
    f=lambda x: math.log10(meet_m(10**x)[0]**2/GR)-x
    xs=list(np.arange(3,13.6,0.25)); v=[f(x) for x in xs]
    roots=[brentq(f,xs[i],xs[i+1],xtol=1e-6) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
    for r in roots:
        F,s=meet_m(10**r); print(f"  D {MD:.1e}: pairs {10**r:.3e} F {F:.3e} sin2 {s:.5f} {(s-S2)/SIG:+.1f} sigma")
