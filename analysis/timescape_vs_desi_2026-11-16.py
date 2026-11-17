"""Timescape (voids + walls, Buchert averaging) against DESI DR2 BAO.
Tracker solution: Dam, Heinesen & Wiltshire 2017, MNRAS 472, 835, Appendix A (eqs A3-A8),
from Wiltshire 2007b, 2009. Data: DESI DR2 BAO table (arXiv:2503.14738, Table IV),
typed from the published table (cloud access to arXiv blocked; re-check digits if used publicly).
Fits: amplitude A = c/(H0 r_d) is free in every model (as DESI does), so only the SHAPE of
the distance-redshift curve is tested. BAO caveat: DESI's fiducial template assumes FLRW;
timescape needs its own BAO reanalysis (not done here)."""
import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar

# --- DESI DR2 BAO (z, kind, value, sigma) and DM-DH correlations
pts = [(0.295,'DV',7.942,0.075)]
pairs = [(0.510,13.588,0.167,21.863,0.425,-0.459),
         (0.706,17.351,0.177,19.455,0.330,-0.404),
         (0.934,21.576,0.152,17.641,0.193,-0.416),
         (1.321,27.601,0.318,14.176,0.221,-0.434),
         (1.484,30.512,0.760,12.817,0.516,-0.500),
         (2.330,38.988,0.531,8.632,0.101,-0.431)]

# --- timescape tracker distances, units c/H0 (dressed)
def ts_D(z_arr, fv0):
    Hb = 1.0
    b = 2*(1-fv0)*(2+fv0)/(9*fv0*Hb)
    def fv(t): return 3*fv0*Hb*t/(3*fv0*Hb*t+(1-fv0)*(2+fv0))
    def zp1(t): return (2+fv(t))*fv(t)**(1/3)/(3*fv0**(1/3)*Hb*t)
    t0 = brentq(lambda t: zp1(t)-1, 1e-6, 10)
    def F(t):
        c=b**(1/3); s=t**(1/3)
        return 2*s + c/6*np.log((s+c)**2/(s*s-c*s+c*c)) + c/np.sqrt(3)*np.arctan((2*s-c)/(np.sqrt(3)*c))
    H0_over_Hb = (4*fv0**2+fv0+4)/(2*(2+fv0))
    out=[]
    for z in np.atleast_1d(z_arr):
        t = brentq(lambda t: zp1(t)-(1+z), 1e-12, t0) if z>0 else t0
        dA = t**(2/3)*(F(t0)-F(t))          # c/Hb
        out.append((1+z)*dA*H0_over_Hb)     # D_M in c/H0
    return np.array(out)

def lcdm_D(z_arr, Om, w0=-1.0, wa=0.0):
    zz=np.linspace(0,max(np.atleast_1d(z_arr))+0.01,4001)
    a=1/(1+zz)
    de=(1-Om)*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a))
    E=np.sqrt(Om*(1+zz)**3+de)
    Dc=np.concatenate([[0],np.cumsum(0.5*(1/E[1:]+1/E[:-1])*np.diff(zz))])
    return np.interp(z_arr,zz,Dc)

def obs(Dfun, z):  # returns DM, DH, DV (units c/H0) with DH = dDM/dz (flat effective)
    h=1e-3
    DM=Dfun(np.array([z]))[0]; DH=(Dfun(np.array([z+h]))[0]-Dfun(np.array([z-h]))[0])/(2*h)
    return DM, DH, (z*DM**2*DH)**(1/3)

def chi2_shape(Dfun):
    # analytic best amplitude A (linear), returns chi2, A
    vecs=[]; 
    z,k,v,s=pts[0]; DM,DH,DV=obs(Dfun,z); vecs.append((np.array([DV]),np.array([v]),np.array([[s*s]])))
    for z,dm,sdm,dh,sdh,r in pairs:
        DM,DH,DV=obs(Dfun,z)
        C=np.array([[sdm**2,r*sdm*sdh],[r*sdm*sdh,sdh**2]])
        vecs.append((np.array([DM,DH]),np.array([dm,dh]),C))
    num=sum(m@np.linalg.solve(C,d) for m,d,C in vecs); den=sum(m@np.linalg.solve(C,m) for m,d,C in vecs)
    A=num/den
    chi=sum((d-A*m)@np.linalg.solve(C,d-A*m) for m,d,C in vecs)
    return chi, A

L=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s); L.append(s)

P("== an earlier test: timescape vs DESI DR2 BAO (13 numbers, shape test, amplitude free) ==")
# LCDM
r=minimize_scalar(lambda Om: chi2_shape(lambda z: lcdm_D(z,Om))[0], bounds=(0.1,0.6), method='bounded')
c_l,A_l=chi2_shape(lambda z: lcdm_D(z,r.x)); P(f"flat LCDM: Om={r.x:.3f} chi2={c_l:.2f} (13-2 dof)  c/(H0 rd)={A_l:.2f}")
# CPL
rc=minimize(lambda p: chi2_shape(lambda z: lcdm_D(z,p[0],p[1],p[2]))[0] if p[1]+p[2]<0 else 1e6, [0.3,-0.8,-0.6], method='Nelder-Mead')
P(f"flat CPL (BAO only, poorly bounded): Om={rc.x[0]:.3f} w0={rc.x[1]:.2f} wa={rc.x[2]:.2f} chi2={rc.fun:.2f} (13-4 dof)")
# timescape
grid=np.linspace(0.50,0.90,81); cs=[chi2_shape(lambda z,f=f: ts_D(z,f))[0] for f in grid]
fbest=grid[int(np.argmin(cs))]
rt=minimize_scalar(lambda f: chi2_shape(lambda z: ts_D(z,f))[0], bounds=(max(0.5,fbest-0.02),min(0.9,fbest+0.02)), method='bounded')
c_t,A_t=chi2_shape(lambda z: ts_D(z,rt.x))
# 1-sigma range on fv0
ok=[f for f,c in zip(grid,cs) if c<=c_t+1.0]
P(f"timescape: fv0={rt.x:.3f} (1-sigma {min(ok):.3f}-{max(ok):.3f}) chi2={c_t:.2f} (13-2 dof)  c/(H0 rd)={A_t:.2f}")
c_s,A_s=chi2_shape(lambda z: ts_D(z,0.737)); P(f"timescape at Seifert+2025 SN value fv0=0.737: chi2={c_s:.2f} (13-1 dof)")
P(f"Delta chi2 (timescape - LCDM), same parameter count: {c_t-c_l:+.2f}")
# per point residual view at best fits
P("\nper-point pulls (data - model)/sigma, LCDM | timescape:")
for (z,dm,sdm,dh,sdh,rr) in pairs:
    a=obs(lambda q: lcdm_D(q,r.x),z); b_=obs(lambda q: ts_D(q,rt.x),z)
    P(f" z={z}: DM {(dm-A_l*a[0])/sdm:+.2f} | {(dm-A_t*b_[0])/sdm:+.2f}   DH {(dh-A_l*a[1])/sdh:+.2f} | {(dh-A_t*b_[1])/sdh:+.2f}")
z,k,v,s=pts[0]; a=obs(lambda q: lcdm_D(q,r.x),z); b_=obs(lambda q: ts_D(q,rt.x),z)
P(f" z={z}: DV {(v-A_l*a[2])/s:+.2f} | {(v-A_t*b_[2])/s:+.2f}")

# what DESI would read if the universe were timescape: fit CPL to noiseless timescape shape
zs=np.linspace(0.1,2.4,24)
def fitcpl(fv0,Om):
    T=ts_D(zs,fv0)
    def loss(p):
        if p[0]+p[1]>=0: return 1e6
        M=lcdm_D(zs,Om,p[0],p[1]); A=(M@T)/(M@M); return np.sum((T-A*M)**2/T**2)
    return minimize(loss,[-0.9,-0.3],method='Nelder-Mead').x
P("\nIf the universe were timescape, a CPL fit (Om fixed) to its distance shape reads:")
for Om in (0.30,0.31,0.32):
    w=fitcpl(rt.x,Om); w2=fitcpl(0.737,Om)
    P(f" Om={Om}: fv0 best -> w0={w[0]:.2f}, wa={w[1]:.2f};  fv0=0.737 -> w0={w2[0]:.2f}, wa={w2[1]:.2f}")
P(" DESI DR2 combined fits: w0 -0.84..-0.67, wa -0.62..-1.09 (wa<0).")

# dressed deceleration q0 and onset
def q(fv): return -(1-fv)*(8*fv**3+39*fv**2-12*fv-8)/(4+fv+4*fv**2)**2
P(f"\ndressed q0 at fv0={rt.x:.3f}: {q(rt.x):.3f}; at 0.737: {q(0.737):.3f} (apparent acceleration begins at fv=0.587)")
open('timescape_vs_desi_2026-11-16.txt','w').write('\n'.join(L)+'\n')
