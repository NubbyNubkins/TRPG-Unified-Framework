"""Hybrid = surface-set dark energy (Lambda) + cosmic-web averaging (voids/walls).
Two-region Buchert model (Buchert 2000; separate-region approach of Rasanen 2006) in volume-average time t:
walls = matter + Lambda, voids = empty + Lambda; both on the timescape tracker at early times (Wiltshire 2007).
Analytic regions: a_w ~ sinh^(2/3)(1.5 sqrt(L) t), a_v ~ sinh(sqrt(L) t). Lapse gamma = Hbar/H_wall; dressed a = abar/gamma;
distances by the Appendix-A recipe of Dam, Heinesen & Wiltshire 2017 (d_A = c a_w(t) Int dt'/(gamma a_w)).
Lambda -> 0 reproduces the tracker (checked); no voids gives flat LCDM. Phenomenological extension (not a published model).
Scaling symmetry: fix L = 1, scan the void weight kappa; fv0 picks 'today'."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq, minimize_scalar
pts=[(0.295,'DV',7.942,0.075)]
pairs=[(0.510,13.588,0.167,21.863,0.425,-0.459),(0.706,17.351,0.177,19.455,0.330,-0.404),
       (0.934,21.576,0.152,17.641,0.193,-0.416),(1.321,27.601,0.318,14.176,0.221,-0.434),
       (1.484,30.512,0.760,12.817,0.516,-0.500),(2.330,38.988,0.531,8.632,0.101,-0.431)]
def model(kappa,L,fv0):
    sL=np.sqrt(L)
    aw_=lambda t: (4/9/L)**(1/3)*np.sinh(1.5*sL*t)**(2/3)
    av_=lambda t: np.sinh(sL*t)/sL
    fvf=lambda t: kappa*av_(t)**3/(aw_(t)**3+kappa*av_(t)**3)
    hi=min(40/sL,1e12)
    if not (fvf(1e-9/max(sL,1e-6)*0+1e-12)<fv0<fvf(hi)): return None
    t0=brentq(lambda t: fvf(t)-fv0,1e-12,hi)
    t=np.geomspace(t0*1e-7,t0,6000); aw,av=aw_(t),av_(t)
    Hw=np.sqrt(4/9/aw**3+L); Hv=np.sqrt(1/av**2+L)
    V=aw**3+kappa*av**3; Hb=(aw**3*Hw+kappa*av**3*Hv)/V; g=Hb/Hw
    a=V**(1/3)/g; z=a[-1]/a-1
    I=cumulative_trapezoid(1/(g*aw),t,initial=0); D=(1+z)*aw*(I[-1]-I)
    o=np.argsort(z); return (lambda q: np.interp(q,z[o],D[o])), Hb[-1]
def chi2_shape(Dfun):
    def obs(z):
        h=1e-3; DM=Dfun(z); DH=(Dfun(z+h)-Dfun(z-h))/(2*h); return DM,DH,(z*DM**2*DH)**(1/3)
    vecs=[]; z,k,v,s=pts[0]; vecs.append((np.array([obs(z)[2]]),np.array([v]),np.array([[s*s]])))
    for z,dm,sdm,dh,sdh,r in pairs:
        DM,DH,_=obs(z); C=np.array([[sdm**2,r*sdm*sdh],[r*sdm*sdh,sdh**2]])
        vecs.append((np.array([DM,DH]),np.array([dm,dh]),C))
    num=sum(m@np.linalg.solve(C,d) for m,d,C in vecs); den=sum(m@np.linalg.solve(C,m) for m,d,C in vecs)
    A=num/den; return sum((d-A*m)@np.linalg.solve(C,d-A*m) for m,d,C in vecs)
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
m=model(1.0,1e-10,0.685); P("check tiny Lambda, fv0=0.685: chi2 =",round(chi2_shape(m[0]),2),"(an earlier test analytic tracker: 15.67)")
P("\nscan: void fraction today fv0; Lambda share chosen by the fit (bare Omega_L = Lambda/(3 Hbar0^2))")
best=(1e9,None); res={}
for fv0 in (0.02,0.10,0.30,0.50,0.56,0.60,0.62,0.64,0.66,0.685,0.737):
    def f(lk):
        mm=model(10**lk,1.0,fv0); return 1e6 if mm is None else chi2_shape(mm[0])
    grid=np.linspace(-9,6,61); vals=[f(x) for x in grid]; j=int(np.argmin(vals))
    r=minimize_scalar(f,bounds=(grid[max(j-1,0)],grid[min(j+1,60)]),method='bounded',options={'xatol':0.005})
    mm=model(10**r.x,1.0,fv0); OL=1.0/mm[1]**2
    m0=model(1.0,1e-10,fv0); c0=chi2_shape(m0[0])
    res[fv0]=(r.fun,OL); P(f" fv0={fv0:.3f}: best chi2={r.fun:.2f} at bare Omega_L={OL:.3f}   (no Lambda: {c0:.2f})")
    if r.fun<best[0]: best=(r.fun,(fv0,OL))
P(f"\noverall best: chi2={best[0]:.2f} at fv0={best[1][0]}, Omega_L={best[1][1]:.2f}   (flat LCDM 10.54 with 1 shape number; hybrid has 2)")
P(f"Seifert SN value fv0=0.737 with Lambda allowed: chi2={res[0.737][0]:.2f} at Omega_L={res[0.737][1]:.2f}")
open('timescape_hybrid_2026-11-16.txt','w').write('\n'.join(out)+'\n')
