# Bottom surface (holographic dark energy, Li 2004, future event horizon) vs DESI DR2 (arXiv:2503.14738)
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import chi2 as C2, norm
OmL0=0.685
def hde(c):
    # dOm/dlna = Om(1-Om)(1+2 sqrt(Om)/c)  (matter + DE, radiation neglected); w = -1/3 - 2 sqrt(Om)/(3c)
    f=lambda x,O: O*(1-O)*(1+2*np.sqrt(np.clip(O,0,1))/c)
    s=solve_ivp(f,[0,np.log(0.4)],[OmL0],dense_output=True,rtol=1e-10)
    a=np.linspace(0.5,1,200); O=s.sol(np.log(a))[0]; w=-1/3-2*np.sqrt(O)/(3*c)
    A=np.vstack([np.ones_like(a),1-a]).T; w0,wa=np.linalg.lstsq(A,w,rcond=None)[0]
    return w0,wa,w[-1]
# DESI DR2 BAO+CMB+SN (w0, sig_w0, wa, sig_wa (symmetrised), quoted LCDM tension in sigma)
D={'PantheonPlus':(-0.838,0.055,-0.62,0.205,2.8),'Union3':(-0.667,0.088,-1.09,0.29,3.8),'DESY5':(-0.752,0.057,-0.86,0.215,4.2)}
def chi(p,d,r):
    w0,s0,wa,sa,_=d; x=np.array([p[0]-w0,p[1]-wa]); Cv=np.array([[s0*s0,r*s0*sa],[r*s0*sa,sa*sa]])
    return float(x@np.linalg.solve(Cv,x))
def sig(c2): return norm.isf(C2.sf(c2,2)/2)
print("HDE today w and CPL fit (a 0.5-1):")
for c in (0.8,0.9,1.0,1.1,1.2,1.5): w0,wa,wt=hde(c); print(f"  c={c}: w(today)={wt:.3f}  w0={w0:.3f} wa={wa:+.3f}")
print("\nDirection: in HDE w falls with time (wa>0); DESI prefers w rising with time (wa<0).")
for n,d in D.items():
    target=C2.isf(2*norm.sf(d[4]),2)
    try: r=brentq(lambda r: chi((-1,0),d,r)-target,-0.99,0.0)   # correlation reproducing DESI's quoted LCDM tension
    except ValueError: r=-0.65; print(f'  ({n}: quoted tension not reproducible with these rounded errors; using correlation -0.65)')
    best=min(((chi(hde(c)[:2],d,r),c) for c in np.linspace(0.5,2.5,81)))
    print(f"\n{n}: correlation (calibrated) {r:.2f}; LCDM {sig(chi((-1,0),d,r)):.1f} sigma")
    for c in (1.0,):
        print(f"  HDE c=1: {sig(chi(hde(c)[:2],d,r)):.1f} sigma")
    print(f"  HDE best c={best[1]:.2f}: {sig(best[0]):.1f} sigma")
