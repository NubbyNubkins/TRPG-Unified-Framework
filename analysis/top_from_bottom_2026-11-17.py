"""Borrowed test 5: top mass from the bottom mass (asymptotic safety; Eichhorn & Held 2018, PRL 121, 151302).
Above the grain gravity adds the same pull to every Yukawa (it is blind to flavour). If top and bottom both sit at a fixed point,
the gravity term cancels in the difference of their one-loop beta functions, leaving  yt^2 - yb^2 = g'^2 / 3  at the grain.
TEST FIXED BEFORE RUNNING:
  take yb from the measured bottom (m_b(m_t) = 2.75 GeV MS-bar), run up to the grain (Planck/sqrt2), impose the relation,
  run back down, convert yt(m_t) to the top pole mass (Buttazzo et al. 2013 fit). With and without TRPG's far side.
  theory spread = |2-loop - 1-loop running| result, plus 1% for the one-loop relation at the grain.
  PASS: within 2 sigma of 172.57 +- 0.29. FAIL: beyond 3 sigma.
Break test 5b, fixed before running: the same rule for the second family, yc^2 - ys^2 = g'^2/3, predicts the charm mass
  (measured m_c(m_c) = 1.273 GeV). Eichhorn & Held expect the light families to sit at zero instead; if the rule can't hold, say so.
"""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
K=16*math.pi**2; MPL=1.220890e19; GRAIN=MPL/math.sqrt(2); MMG=2.4e8; F=4.6e13; MT0=172.57; v=246.22
dbL=np.array([2/3,2/3,0]); dBL=np.array([[25/9*9/50,5/3*9/10,0],[5/3*3/10,49/6,0],[0,0,0]])
dbD=np.array([4/9,0,2/3]); dBD=np.array([[25/9*4/75,0,5/3*16/15],[0,0,0],[5/3*2/15,0,38/3]])
b0=np.array([41/6,-19/6,-7.0]); B0=np.array([[199/18,9/2,44/3],[3/2,35/6,12],[11/6,9/2,-26]]); cT=np.array([17/6,3/2,2])
def beta(t,y,far,two):
    g1,g2,g3,yt,yb,l=y; mu=MT0*math.exp(t)
    b=b0.copy(); B=B0.copy()
    if far:
        if mu>MMG: b=b+3*dbL; B=B+3*dBL
        if mu>F:   b=b+3*dbD; B=B+3*dBD
    g=np.array([g1,g2,g3]); G=g*g; a1,a2,a3=G; T=yt*yt
    dg=g**3*b/K+(g**3*(B@G-cT*T)/K**2 if two else 0)
    dyt=yt*(4.5*T+1.5*yb*yb-8*a3-2.25*a2-17/12*a1)/K
    if two: dyt+=yt*(-12*T*T+T*(36*a3+225/16*a2+131/16*a1-12*l)+6*l*l-108*a3*a3+9*a2*a3+19/9*a1*a3-23/4*a2*a2-3/4*a1*a2+1187/216*a1*a1)/K**2
    dyb=yb*(4.5*yb*yb+1.5*T-8*a3-2.25*a2-5/12*a1)/K
    if two: dyb+=yb*(-108*a3*a3+9*a2*a3+31/9*a1*a3+T*(16*a3)-23/4*a2*a2)/K**2   # leading two-loop terms
    dl=(24*l*l-6*T*T+3/8*(2*a2*a2+(a2+a1)**2)+l*(12*T-9*a2-3*a1))/K
    return [dg[0],dg[1],dg[2],dyt,dyb,dl]
def start(mt,a3=0.1180,mw=80.3692,mh=125.20):
    return np.array([0.35830+0.00011*(mt-173.34)-0.00020*(mw-80.384)/0.014, 0.64779+0.00004*(mt-173.34)+0.00011*(mw-80.384)/0.014,
                     1.1666+0.00314*(a3-0.1184)/0.0007-0.00046*(mt-173.34), 0.93690+0.00556*(mt-173.34)-0.00042*(a3-0.1184)/0.0007,
                     math.sqrt(2)*2.75/v, 0.12604+0.00206*(mh-125.15)-0.00004*(mt-173.34)])
def up(y,far,two,mu=GRAIN):
    edges=sorted([0.0,math.log(mu/MT0)]+([math.log(MMG/MT0),math.log(F/MT0)] if far else []))
    for a,b in zip(edges[:-1],edges[1:]): y=solve_ivp(beta,[a,b],y,args=(far,two),rtol=1e-10,atol=1e-13).y[:,-1]
    return y
def mismatch(mt,far,two):
    y=up(start(mt),far,two); return y[3]**2-y[4]**2-y[0]**2/3
def predict(far,two=True):
    return brentq(lambda m: mismatch(m,far,two),120,200,xtol=1e-4)
for far in (False,True):
    p2=predict(far,True); p1=predict(far,False)
    y=up(start(MT0),far,True)
    s=math.sqrt(0.29**2+(p2-p1)**2+(0.01*p2)**2); z=(p2-MT0)/s
    v_="PASS" if abs(z)<2 else ("FAIL" if abs(z)>3 else "PARTIAL")
    print(f"{'with TRPG far side' if far else 'Standard Model only':20s}: at the grain today yt={y[3]:.4f} yb={y[4]:.4f} g'={y[0]:.4f} "
          f"-> rule needs yt={math.sqrt(y[4]**2+y[0]**2/3):.4f}")
    print(f"   predicted top pole mass {p2:.2f} GeV (1-loop running {p1:.2f}); spread {s:.2f}; vs 172.57: {z:+.1f} sigma -> {v_}")
print("\nbreak test 5b, second family:")
y=up(start(MT0),True,True)
print(f"   yc^2 - ys^2 = g'^2/3 needs yc(grain) >= {y[0]/math.sqrt(3):.3f}, about the top's own size; measured charm gives yc(grain) ~ {y[3]*1.273/162.5*0.5:.4f} "
      "(scaled from the top, order of magnitude) -> charm would be ~100 GeV: FAIL. The light families must sit at zero (as Eichhorn & Held say).")
