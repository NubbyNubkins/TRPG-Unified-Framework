"""Fewer-particle readings: numerical checks (19 Nov 2026). Exploration only.
Inputs: PDG 2026 central values as used in the v1.3.3 recalibration."""
import math
from fractions import Fraction as F
pi=math.pi
print("="*70); print("1. HARD COUNTS FROM DATA")
# Z invisible width -> number of light neutrinos
G_inv=499.0; G_inv_err=1.5   # MeV, LEP
G_nu_SM=167.21              # MeV per neutrino (SM, with rad. corr.)
N=G_inv/G_nu_SM; print(f"  Z invisible width: N_nu = {N:.3f} +/- {G_inv_err/G_nu_SM:.3f}  (LEP lineshape fit 2.9963+/-0.0074)")
# R ratio -> colours
for Nc in (1,2,3,4):
    R=Nc*(4/9+1/9+1/9+4/9+1/9)*(1+0.118/pi)
    print(f"  R(e+e- -> hadrons, 5 quarks) with {Nc} colour(s): {R:.2f}   measured ~3.85")
# pi0 -> 2 gamma (rate ~ Nc^2)
G3=7.78
for Nc in (1,3): print(f"  pi0 width with {Nc} colours: {G3*(Nc/3)**2:.2f} eV  measured 7.80+/-0.12")
# gluon self-coupling: colour factors at LEP (OPAL/ALEPH combined approx)
print("  LEP colour factors: C_A = 2.89+/-0.21 (SU(3): 3), C_F = 1.30+/-0.09 (SU(3): 4/3)")
print("  -> 3 light neutrinos, 3 colours, 8 self-interacting gluons are measured, not assumed.")

print("="*70); print("2. ANOMALY CANCELLATION: is one generation reducible?")
# Left-handed Weyl fields (hypercharge Y, Q = T3 + Y), with colour multiplicity
gen=[("Q",3,2,F(1,6)),("u^c",3,1,F(-2,3)),("d^c",3,1,F(1,3)),("L",1,2,F(-1,2)),("e^c",1,1,F(1,1))]
def anomalies(fields):
    a={}
    a["Y^3"]=sum(c*w*Y**3 for _,c,w,Y in fields)
    a["grav-Y"]=sum(c*w*Y for _,c,w,Y in fields)
    a["SU2^2 Y"]=sum(c*Y for _,c,w,Y in fields if w==2)
    a["SU3^2 Y"]=sum(w*Y for _,c,w,Y in fields if c==3)
    a["SU2 doublets (Witten, must be even)"]=sum(c for _,c,w,Y in fields if w==2)
    return a
def show(name,fields):
    a=anomalies(fields); ok=all(v==0 for k,v in a.items() if "Witten" not in k) and a["SU2 doublets (Witten, must be even)"]%2==0
    print(f"  {name:38s} "+", ".join(f"{k.split(' (')[0]}={v}" for k,v in a.items())+f"  -> {'CONSISTENT' if ok else 'BROKEN'}")
show("full generation (15 states)",gen)
show("quarks only",[f for f in gen if f[1]==3])
show("leptons only",[f for f in gen if f[1]==1])
show("drop right-handed electron",[f for f in gen if f[0]!="e^c"])
show("drop down-type quark (u only)",[f for f in gen if f[0]!="d^c"])
print("  -> A generation is ONE unit: quarks alone or leptons alone are inconsistent. Supports 'one object', blocks dropping parts.")

print("="*70); print("3. FORCE STRENGTHS: does REMOVING particles help them meet?")
aem=1/127.951; s2=0.23122; a3=0.1180
a1=(5/3)*aem/(1-s2); a2=aem/s2
need=(1/a2-1/a3)/(1/a1-1/a2)
print(f"  1/alpha at MZ: {1/a1:.2f}, {1/a2:.2f}, {1/a3:.2f}.  Meeting needs B = (b2-b3)/(b1-b2) = {need:.3f}")
def B(b): return (b[1]-b[2])/(b[0]-b[1])
SM=(41/10,-19/6,-7.)
cases={"Standard Model content":SM,
 "remove the Higgs doublet":(41/10-1/10,-19/6-1/6,-7),
 "gluons as quark links (no gluon self-loop)":(41/10,-19/6,-7+11),
 "remove a whole generation":(41/10-4/3,-19/6-4/3,-7-4/3),
 "quarks without colour (1 copy each)":None,
 "ADD a 2nd Higgs doublet":(41/10+1/10,-19/6+1/6,-7),
 "MSSM (doubles the particles)":(33/5,1.,-3.)}
# quarks without colour: per generation quark contributions divided by 3
qg=(F(4,3)*F(1,1),0,0)  # placeholder not used
for k,b in cases.items():
    if b is None:
        # per generation: quark part of b1=(4/3)*(17/20+... ) compute directly
        # b_i(fermions)=4/3 per generation split: b1 quarks 11/15*... use explicit sums
        def bf(colour):
            # Dirac-equivalent sums of (Y^2) in GUT norm: b1 += (2/3)*(3/5)*sum Y^2 over Weyl; b2 += (2/3)*(1/2)*#doublets; b3 += (2/3)*(1/2)*#triplets
            fl=[("Q",colour,2,1/6),("u",colour,1,2/3),("d",colour,1,1/3),("L",1,2,1/2),("e",1,1,1)]
            b1=sum((2/3)*(3/5)*c*w*Y*Y for _,c,w,Y in fl)
            b2=sum((2/3)*0.5*c for _,c,w,Y in fl if w==2)
            return b1,b2
        b1q,b2q=bf(1); b=(3*b1q+1/10,-22/3+3*b2q+1/6,-11.+0)
        k+= " [b3: no quarks feel colour]"
    d=abs(B(b)-need)-abs(B(SM)-need)
    tag='(reference)' if b==SM else ('same as SM' if abs(d)<1e-9 else ('closer than SM' if d<0 else 'further than SM'))
    print(f"  {k:52s} B = {B(b):.3f}   {tag}")
print("  -> Removals that data allow make it worse or leave it unchanged. The one removal that helps (colourless quarks)")
print("     is excluded by the R ratio and pi0 width above. ADDING states (2HDM, MSSM, TRPG's far-side MMG) is what moves it closer.")
# alpha_s running test of 'no gluon self-interaction'
print("  Strong-force running check (one loop, 5 flavours throughout, so only rough): alpha_s(m_tau)=0.315 measured, alpha_s(MZ)=0.118.")
def run(a0,b0,mu0,mu): return 1/(1/a0+b0/(2*pi)*math.log(mu/mu0))
print(f"    with gluon self-loops (b0=11-2*5/3=7.67): 0.315 at 1.777 GeV -> {run(0.315,11-10/3,1.777,91.19):.3f} at MZ")
print(f"    gluons as passive links (b0=-2*5/3):        0.315 at 1.777 GeV -> {run(0.315,-10/3,1.777,91.19):.3f} at MZ  (grows; FAILS)")

print("="*70); print("4. ONE THING IN THREE STATES: Koide test (three generations as one object)")
def koide(m): return sum(m)/sum(math.sqrt(x) for x in m)**2
lep=(0.51099895,105.6583755,1776.93)
print(f"  charged leptons Q = {koide(lep):.6f}  (2/3 = 0.666667)  -> holds to {abs(koide(lep)-2/3)/(2/3)*1e6:.0f} ppm")
up=(2.16e-3,1.2730,172.57); dn=(4.70e-3,0.0935,4.183)
print(f"  up quarks       Q = {koide(up):.4f}   down quarks Q = {koide(dn):.4f}  -> not 2/3 (quarks don't share the leptons' pattern)")

print("="*70); print("5. HIGGS AS A TOP-ANTITOP PAIR (top condensation, Bardeen-Hill-Lindner 1990)")
# one-loop SM RGEs, compositeness condition y_t -> large at Lambda
def gauge(t):  # t = ln(mu/MZ)
    b=(41/10,-19/6,-7.); a0=(a1,a2,a3)
    return [ (1/(1/a0[i]-b[i]/(2*pi)*t)) for i in range(3)]
for L in (1e6,1e10,1e15,1e19):
    tL=math.log(L/91.19); y=8.0; n=20000; dt=-tL/n; t=tL
    for _ in range(n):
        al=gauge(t); g1s,g2s,g3s=[4*pi*x for x in al]
        dy=y/(16*pi**2)*(4.5*y*y-8*g3s-2.25*g2s-(17/20)*g1s)
        y+=dy*dt; t+=dt
    mt=y*246.22/math.sqrt(2)
    print(f"  compositeness scale {L:.0e} GeV -> top mass {mt:.0f} GeV (measured 172.6); Higgs ~2 m_t in simplest form = {2*mt:.0f} (measured 125.1)")
print("  -> minimal 'Higgs = top pair' FAILS (top too heavy, Higgs far too heavy).")

print("="*70); print("6. PREON CHARGE PATTERN (Harari-Shupe rishons: T charge +1/3, V charge 0)")
for name,comb in [("positron","TTT"),("up","TTV"),("anti-down","TVV"),("neutrino","VVV")]:
    q=sum(F(1,3) if c=="T" else 0 for c in comb); print(f"  {name:10s} = {comb}: charge {q}")
print("  -> one generation's charges from 2 pieces: HOLDS. Compositeness limits (LHC contact terms) push binding above ~20-30 TeV,")
print("     while the electron weighs 0.5 MeV: binding 10^7-10^8 times its weight. Generations unexplained. Dynamics FAIL.")

print("="*70); print("7. LIGHT SPECIES IN THE EARLY UNIVERSE")
print("  Planck+BAO N_eff = 2.99 +/- 0.17 (SM 3.044); BBN helium agrees. Fewer light species than the SM: excluded at ~6 sigma per missing neutrino.")
