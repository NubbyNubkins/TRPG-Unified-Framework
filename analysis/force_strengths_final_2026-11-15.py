"""OPEN ITEM 1, the three force strengths with gravity tied in: final test of every route.
Two-loop running of the three strengths (MS-bar, Standard Model + far-side parts of the E6 27: three spin-1/2 weak pairs L and three
coloured D partners), step thresholds, top decoupled below its mass.
Inputs (PDG 2024): alpha_s(MZ) = 0.1180 +- 0.0009, sin^2 theta_W(MZ, MS-bar) = 0.23122 +- 0.00004, 1/alpha(MZ) = 127.951 +- 0.009,
m_t = 172.57 GeV, y_t(m_t) = 0.9369 (Buttazzo and others 2013). Planck mass: UNREDUCED 1.220890e19 GeV is the main footing
(the 'grain' of an earlier test); the reduced one (2.435e18) is run beside it and labelled.
Credit: Georgi, Quinn & Weinberg 1974; Machacek & Vaughn 1983-85; Jones 1982; Arvanitaki, Davis, Graham & Wacker 2004; Giudice & Romanino 2004;
Gursey, Ramond & Sikivie 1976; Chamseddine & Connes 1996; Kaplunovsky 1988; Cirelli, Fornengo & Strumia 2006; Griest & Kamionkowski 1990;
Super-Kamiokande 2020; LZ 2024; PDG 2024.
Re-run: python3 force_strengths_final_2026-11-15.py > force_strengths_final_2026-11-15.txt   (deterministic, seeded; under a minute)"""
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

# ------------------------------------------------------------------------------------------------------------------
head("1. The Standard Model alone: one number for the miss, and the three older numbers reconciled")
pr("Method A (the standard one): inputs = electric strength and strong strength at the Z; demand that all three lines meet at one energy;")
pr("that fixes the weak mixing number sin^2 at the Z. Compare with the measured 0.23122 +- 0.00004.")
res={}
for two in (False,True):
    for td in (False,True):
        s=predict_s2(two=two,topdec=td); res[(two,td)]=s
        pr(f"  {'two loops' if two else 'one loop '}, top {'decoupled below m_t' if td else 'kept from the Z   '}: sin^2 = {s:.5f}  ({pct(s):+.2f}%)")
s0=res[(True,True)]
ds_a3=(predict_s2(a3=A3+dA3)-predict_s2(a3=A3-dA3))/2
ds_ai=(predict_s2(ainv=AINV+dAINV)-predict_s2(ainv=AINV-dAINV))/2
ds_yt=(predict_s2(yt=YT+dYT)-predict_s2(yt=YT-dYT))/2
ds_top=res[(True,True)]-res[(True,False)]
twoshift=res[(True,True)]-res[(False,True)]
ds_3loop=abs(twoshift)*abs(twoshift)/(0.375-res[(False,True)])*3   # (two-loop step)^2 / (one-loop total), times 3 for safety
ds_3loop=max(ds_3loop,0.0002)
sig_th=math.sqrt(ds_a3**2+ds_ai**2+ds_yt**2+ds_top**2+ds_3loop**2)
sig_tot=math.sqrt(sig_th**2+dS2**2)
pr(f"\nError budget on the predicted sin^2 (two loops):")
pr(f"  alpha_s +-0.0009            : {ds_a3:+.5f}")
pr(f"  1/alpha +-0.009             : {ds_ai:+.5f}")
pr(f"  top link +-0.005            : {ds_yt:+.5f}")
pr(f"  top threshold treatment     : {ds_top:+.5f}  (full difference between the two treatments)")
pr(f"  three loops, not computed   : +-{ds_3loop:.5f}  (estimate: (two-loop step {twoshift:+.4f})^2 / one-loop total, x3, floor 0.0002)")
pr(f"  total theory                : +-{sig_th:.5f};  measurement +-{dS2:.5f};  combined +-{sig_tot:.5f}")
pr(f"\n  ONE NUMBER: Standard Model alone predicts sin^2 = {s0:.4f} +- {sig_th:.4f}; measured 0.23122. Miss = {pct(s0):+.1f}% +- {100*sig_th/S2:.1f}%")
pr(f"              = {(s0-S2)/sig_tot:+.0f} combined sigma. FAIL.")
pr("\nThe same miss seen three other ways (same running, two loops, top decoupled):")
m=mismatch(); pr(f"  Method C (an earlier test's): use the MEASURED sin^2, find where strong = weak ({m[1]:.2e} GeV, 1/alpha = {m[2]:.1f}); there hypercharge is off by {m[0]:+.2f} in 1/alpha.")
m1=mismatch(two=False); pr(f"     one loop: {m1[0]:+.2f} at {m1[1]:.2e} GeV. an earlier test turned this into a sin^2 at the Z of 0.2067 (-10.6%): a mixed method, retired.")
r=run(stop='12'); pr(f"  Triangle: where weak = hypercharge ({r[0]:.2e} GeV, 1/alpha = {ai(r[1][0]):.1f}) the strong line is off by {ai(r[1][2])-ai(r[1][0]):+.2f} in 1/alpha (two loops)")
r=run(stop='12',two=False,topdec=False); pr(f"     one loop, top kept from the Z (the earlier number): {ai(r[1][2])-ai(r[1][0]):+.2f} at {r[0]:.2e} GeV")
def pred_a3(two=True):
    def f(a3):
        M,g=run(a3=a3,stop='12',two=two); return ai(g[2])-ai(g[0])
    return brentq(f,0.04,0.2)
pr(f"  Method B: inputs alpha and sin^2; predict the strong strength at the Z: {pred_a3():.4f} against 0.1180 +- 0.0009 ({100*(pred_a3()/A3-1):+.0f}%)")
pr("  Reconciled: -10.2% = Method A at one loop; -8.9% = Method A at two loops (top kept from the Z; -8.9% also with the top decoupled);")
pr("              -10.6% = the earlier mixed method. The number to quote is Method A at two loops.")

# ------------------------------------------------------------------------------------------------------------------
head("2. Far-side spin-1/2 weak pairs (three, one per placement): the mass that closes the triangle, with its error")
ML0=ML_close(); mm=mismatch(thr=thrL(ML0)); MU0=mm[1]; AU0=mm[2]
pr(f"D partners at (or above) the meeting. Two loops, top decoupled: M_L = {ML0:.3e} GeV; meeting at {MU0:.3e} GeV; 1/alpha there = {AU0:.2f}")
for two,td,lab in ((False,False,"one loop, top from the Z"),(True,False,"two loops, top from the Z"),(False,True,"one loop, top decoupled")):
    x=ML_close(two=two,topdec=td); pr(f"  {lab:36s}: M_L = {x:.3e} GeV")
slope=(predict_s2(thr=thrL(ML0*10**0.5))-predict_s2(thr=thrL(ML0/10**0.5)))
pr(f"Slope: the predicted sin^2 changes by {slope:+.5f} for each factor 10 in M_L.")
dl=[]
for lab,kwp,kwm in (("alpha_s +-0.0009",dict(a3=A3+dA3),dict(a3=A3-dA3)),("sin^2 +-0.00004",dict(s2=S2+dS2),dict(s2=S2-dS2)),
                    ("1/alpha +-0.009",dict(ainv=AINV+dAINV),dict(ainv=AINV-dAINV)),("top link +-0.005",dict(yt=YT+dYT),dict(yt=YT-dYT))):
    d=(math.log10(ML_close(**kwp))-math.log10(ML_close(**kwm)))/2; dl.append(d); pr(f"  {lab:20s}: {d:+.3f} decades in M_L")
d=math.log10(ML_close(topdec=True))-math.log10(ML_close(topdec=False)); dl.append(d); pr(f"  top threshold       : {d:+.3f} decades")
d=ds_3loop/abs(slope); dl.append(d); pr(f"  three loops (est.)  : +-{d:.3f} decades")
sigL=math.sqrt(sum(x*x for x in dl))
pr(f"  TOTAL: log10(M_L) = {math.log10(ML0):.2f} +- {sigL:.2f}  ->  M_L = {ML0:.2e} GeV, range {ML0/10**sigL:.2e} to {ML0*10**sigL:.2e} GeV (1 sigma)")
pr(f"         2 sigma: {ML0/10**(2*sigL):.2e} to {ML0*10**(2*sigL):.2e} GeV")
pr("\nFewer light pairs (the other placements' pairs at the meeting):")
for n in (1,2):
    x=ML_close(n=n); pr(f"  {n} light pair(s): " + (f"M_L = {x:.2e} GeV" if x else "no mass above 3 GeV closes it"))
pr("\nTHE HIDDEN SECOND NUMBER: where the D partners sit. If they are lighter than the meeting, the weak pairs must be lighter too:")
pr(f"  {'M_D (GeV)':>10s} {'M_L needed':>11s} {'M_D / M_L':>10s} {'meeting (GeV)':>14s} {'1/alpha there':>13s}")
ratios=[]
for MD in (3e13,1e13,1e12,1e11,1e10,1e9,1e8,4.2e7,2.1e7,1e7,1e6,5.28e5,2.5e5):
    x=ML_close(MD=MD)
    if x:
        q=mismatch(thr=thrLD(x,MD)); ratios.append(MD/x); pr(f"  {MD:10.2e} {x:11.3e} {MD/x:10.3e} {q[1]:14.3e} {q[2]:13.2f}")
    else: pr(f"  {MD:10.2e}   none above 3 GeV")
RAT=np.exp(np.mean(np.log(ratios)))
pr(f"  The strengths fix only the RATIO: M_D / M_L = {RAT:.2e} (5.3 decades), anywhere from M_L = about 100 GeV to {ML0:.1e} GeV.")
pr(f"  (If the D partners are at or above the meeting, M_L = {ML0:.2e} GeV and the 'ratio' is meeting/M_L = {MU0/ML0:.2e}.)")
pr("  Bounds: charged weak-pair partners above 103.5 GeV (LEP), roughly 200 GeV for pure pairs (LHC disappearing tracks);")
pr(f"  so M_D must be above about {RAT*103.5:.1e} to {RAT*200:.1e} GeV for the triangle to close.")
v=246.22; p0=2.873/3
for k,lab in ((0.25,"k = 0.25, link = balance"),(0.26,"phased fold k = 0.26"),(0.12,"real fold k = 0.12")):
    MD=v**2/(k*p0); pr(f"  The fold's D with link = Higgs value ({lab}): M_D = {MD:.2e} GeV -> weak pairs would need M_L = {MD/RAT:.2f} GeV: EXCLUDED (LEP). 'Both jobs' FAILS for every M_L.")

# ------------------------------------------------------------------------------------------------------------------
head("3. Checks at the meeting point and on the pairs themselves")
hbar=6.582119569e-25; yr=3.15576e7; mp=0.93827
tau=lambda M,ainvU: M**4*ainvU**2/mp**5*hbar/yr
pr(f"Meeting: {MU0:.2e} GeV, 1/alpha = {AU0:.1f}. Above it the three lines part again unless the theory changes there (no merged force in the main line).")
r=run(thr=thrLD(ML0,MU0),stop=MPL); g=r[1]
pr(f"  Carried on to the grain (unreduced Planck mass) with the whole far side counted: 1/alpha = {ai(g[0]):.1f} (hypercharge), {ai(g[1]):.1f} (weak), {ai(g[2]):.1f} (strong)")
pr(f"If carriers that turn quarks into leptons sat at the meeting (as in Georgi-Glashow): proton lifetime about M^4/(alpha^2 m_p^5) = {tau(MU0,AU0):.1e} years")
pr(f"  (rough to a factor 10). Super-Kamiokande 2020: p -> e+ pi0 longer than 2.4e34 years. FAIL by a factor {2.4e34/tau(MU0,AU0):.0e}.")
pr(f"  A meeting with such carriers must sit above about {(2.4e34/tau(1.0,AU0))**0.25:.1e} GeV. So the main line survives ONLY because it has no such carriers (43.4).")
pr(f"  Gravity alone (Planck-suppressed): about {tau(MPL,1.0):.0e} years. No signal expected at Hyper-Kamiokande.")
pr("\nThe pairs: a charged member and a neutral member. Loops make the charged one heavier by about 0.34 GeV (Thomas & Wells 1998;")
pr("  Cirelli, Fornengo & Strumia 2006), so it decays to the neutral one and a pion in about 1e-11 s. Does the neutral one decay?")
pr("  It has the same charges as the ordinary lepton pair, so it can mix with it, exactly as the D mixes with the down quarks (the fold).")
for y in (1e-6,1e-10,1e-14):
    G=y*y*ML0/(16*math.pi); pr(f"    mixing link {y:.0e}: lifetime {hbar/G:.1e} s")
ymin=math.sqrt(hbar/0.1*16*math.pi/ML0)
pr(f"  To be gone before the light elements form (0.1 s) the link need only exceed {ymin:.0e}. Any fold-like mixing does it.")
om=0.12*(ML0/1100)**2
pr(f"  If the neutral one were exactly stable and the early universe was ever hotter than M_L: leftover = 0.12 x (M_L / 1.1 TeV)^2 = {om:.0e}")
pr(f"    against the measured dark matter 0.120: too much by {om/0.12:.0e}. (Any stable thermal relic above about 1e5 GeV overcloses: Griest & Kamionkowski 1990.)")
pr(f"    And a Z-coupled relic scatters on nuclei at about 1e-39 cm^2; LZ 2024 allows about {2.2e-48*ML0/40*1.3:.0e} cm^2 at this mass (scaled linearly). Excluded as the dark matter.")
pr(f"  So: either the neutral member decays (fold-like mixing), or the start was colder than about M_L/20 = {ML0/20:.0e} GeV.")
pr("  The cold, non-thermal start of an earlier test (combination 3) already satisfies the second; the fold-like mixing satisfies the first. No conflict.")
pr("  Light-element formation (BBN) and the microwave background: no effect; the pairs are gone or were never made. Colliders: out of reach by 4 decades.")

# ------------------------------------------------------------------------------------------------------------------
head("4. Every fixed rule for the far-side pair mass (D partners at the meeting, the most favourable case)")
me=0.51099895e-3; mmu=0.1056584; mtau=1.77686; mH=125.20; mW=80.3692
mnu=[0.00036e-9,0.00866e-9,0.05034e-9]; N=[me**2/mnu[0],mmu**2/mnu[1],mtau**2/mnu[2]]
def row(name,M=None,thr=None):
    if thr is None:
        if M<3 or M>MU0*0.999: pr(f"  {name:58s} {M:10.3e} GeV   outside the usable range: FAIL"); return None
        thr=thrL(M)
    s=predict_s2(thr=thr); pull=(s-S2)/sig_tot
    lab="closes" if abs(pull)<1 else ("within 2 sigma" if abs(pull)<2 else "FAIL")
    pr(f"  {name:58s} {('%10.3e GeV'%M) if M else '   (three masses)'}   sin^2 = {s:.5f} ({pct(s):+.2f}%, {pull:+.1f} sigma)  {lab}")
    return s
pr(f"Combined error on sin^2: +-{sig_tot:.5f} ({100*sig_tot/S2:.2f}%). 'sigma' below uses it. Exact closure needs {ML0:.2e} GeV.")
pr("\n(i) Electron-gravity balance and its relatives")
row("sqrt(m_e x M_Planck)  [unreduced]",math.sqrt(me*MPL))
row("sqrt(m_e x M_Planck)  [reduced]",math.sqrt(me*MRED))
row("v^2/m_e (v = 246.22)",v**2/me); row("(v/sqrt2)^2/m_e",v**2/2/me); row("top^2/m_e",MT**2/me); row("Higgs^2/m_e",mH**2/me); row("Z^2/m_e",mZ**2/me)
pr("(ii) Other balance points with gravity")
row("sqrt(m_mu x M_Planck)",math.sqrt(mmu*MPL)); row("sqrt(proton x M_Planck)",math.sqrt(mp*MPL)); row("sqrt(v x M_Planck)",math.sqrt(v*MPL))
row("sqrt(m_nu3 x M_Planck)",math.sqrt(mnu[2]*MPL)); row("v^(2/3) x M_Planck^(1/3)",v**(2/3)*MPL**(1/3)); row("v^2/M_Planck",v**2/MPL)
pr("(iii) The hidden neutrino partners (earlier tests: 7.2e5, 1.3e9, 6.3e10 GeV)")
for i in range(3): row(f"all three pairs at hidden partner {i+1}",N[i])
row("geometric mean of the three hidden partners",(N[0]*N[1]*N[2])**(1/3))
sN=row("each pair at its own placement's hidden partner",thr=[(N[i],1,dbL,dBL) for i in range(3)])
for f in (0.5,2.0):
    Nv=[me**2/(mnu[0]*f),N[1],N[2]]; s=predict_s2(thr=[(Nv[i],1,dbL,dBL) for i in range(3)])
    pr(f"     the same with the lightest neutrino x{f} (it is TRPG's prediction, not measured): sin^2 = {s:.5f} ({pct(s):+.2f}%, {(s-S2)/sig_tot:+.1f} sigma)")
pr("(iv) Seesaw with the D partner, and the D's own mass from the fold")
for k in (0.25,0.12):
    MD=v**2/(k*p0); row(f"v^2/M_D with the fold's D (k = {k})",v**2/MD)
    s=predict_s2(thr=thrLD(MD,MD)); pr(f"  pairs AT the fold's D mass {MD:.2e} GeV (whole set together, k = {k}):            sin^2 = {s:.5f} ({pct(s):+.2f}%): no better than the Standard Model. FAIL")
    row(f"pairs at the fold's D mass (k = {k}), D's moved to the meeting",MD)
pr(f"  v^2/M_D together with the ratio rule M_D = {RAT:.1e} M_L: M_L = v/sqrt(ratio) = {v/math.sqrt(RAT):.2f} GeV: EXCLUDED (LEP). FAIL")
pr("(v) Gravity seesaw of the meeting itself: M_L = (meeting)^2 / M_Planck, solved together with the running (no input mass at all)")
def meetA(ML):
    s_=predict_s2(thr=thrL(ML)); q=mismatch(s2=s_,thr=thrL(ML)); return q[1]
def fixed(MP):
    f=lambda x: math.log10(meetA(10**x)**2/MP)-x     # Method A: the lines are made to meet, so 'the meeting' is unambiguous
    return 10**brentq(f,5,12,xtol=1e-6)
FIX={}
for lab,MP in (("unreduced",MPL),("reduced",MRED)):
    M=fixed(MP); FIX[lab]=M; row(f"M_L = M_meet^2/M_Planck [{lab}], self-consistent",M)
pr(f"  The same thing read the other way: at exact closure, meeting = {MU0:.3e}; sqrt(M_L x M_Planck) = {math.sqrt(ML0*MPL):.3e} (unreduced), {math.sqrt(ML0*MRED):.3e} (reduced)")
a_=math.log(MU0/ML0)/math.log(MPL/ML0); ar=math.log(MU0/ML0)/math.log(MRED/ML0)
da=abs(a_-math.log(meetA(ML0*10**sigL)/(ML0*10**sigL))/math.log(MPL/(ML0*10**sigL)))
pr(f"  Position of the meeting between M_L and the grain, on a log scale: {a_:.3f} +- {da:.3f} (unreduced), {ar:.3f} (reduced). One half would be exact.")
pr("\n(v-b) Rules that land BELOW the closing mass can be made exact by moving the D partners below the meeting (M_D = ratio x M_L):")
for nm,M in (("sqrt(m_e x M_Planck) [unreduced]",math.sqrt(me*MPL)),("v^2/m_e",v**2/me),("M_meet^2/M_Planck [unreduced]",FIX["unreduced"]),("top^2/m_e",MT**2/me)):
    f=lambda x: mismatch(thr=thrLD(M,10**x))[0]
    MD=10**brentq(f,math.log10(M)+1,13.9); pr(f"    {nm:34s}: exact if M_D = {MD:.2e} GeV (the meeting is then at {mismatch(thr=thrLD(M,MD))[1]:.2e})")
pr("  Rules that land ABOVE it (the hidden partners) cannot be rescued this way. So no rule for M_L alone is a test until M_D has a rule too.")
pr("\n(v-c) Unequal pair masses: only the product matters (one loop exactly, two loops nearly)")
tot=3*math.log10(MU0/ML0)
for M1,M2 in ((1e6,1e8),(N[0],N[1]),(1e3,1e9)):
    f=lambda x: mismatch(thr=[(M1,1,dbL,dBL),(M2,1,dbL,dBL),(10**x,1,dbL,dBL)])[0]
    try:
        M3=10**brentq(f,2,13.6); q=mismatch(thr=[(M1,1,dbL,dBL),(M2,1,dbL,dBL),(M3,1,dbL,dBL)])
        pr(f"    pairs at {M1:.2e}, {M2:.2e} -> third at {M3:.2e} GeV; sum of log10(meeting/M) = {sum(math.log10(q[1]/m_) for m_ in (M1,M2,M3)):.2f} (equal masses: {tot:.2f})")
    except ValueError: pr(f"    pairs at {M1:.2e}, {M2:.2e}: no third mass closes it")
pr(f"  General statement: the triangle closes when the pair-only stretches add up to {tot:.1f} decades:")
pr(f"     sum over the three placements of log10( min(M_D, meeting) / M_L ) = {tot:.1f}.  One condition on up to six far-side masses.")
pr("\n(v-d) Rules for the RATIO M_D/M_L = %.2e (x/÷ %.2f): ratios of the theory's scales" % (RAT,10**sigL))
Sr={"m_e":me,"m_mu":mmu,"m_tau":mtau,"proton":mp,"Z":mZ,"Higgs":mH,"top":MT,"v":v,"nu2":mnu[1],"nu3":mnu[2],"N1":N[0],"N2":N[1],"N3":N[2],"meeting":MU0,"Planck":MPL}
rat=[(f"{a}/{b}",Sr[a]/Sr[b]) for a in Sr for b in Sr if Sr[a]>Sr[b]]
lr=np.array([math.log10(x) for _,x in rat]); t0=math.log10(RAT)
hits=[(n,x) for n,x in rat if abs(math.log10(x)-t0)<sigL]
rng0=np.random.default_rng(7); T=rng0.uniform(1,12,20000); frac=np.mean([(np.abs(lr-t)<sigL).any() for t in T])
pr(f"    {len(rat)} ratios; within 1 sigma: " + "; ".join(f"{n} = {x:.2e}" for n,x in hits))
pr(f"    luck: a target ratio placed anywhere from 1e1 to 1e12 has at least one such hit {100*frac:.0f}% of the time. NOT EVIDENCE.")
pr("\n(vi) LUCK CHECK over a defined family")
def family(MP):
    S={"m_e":me,"m_mu":mmu,"m_tau":mtau,"proton":mp,"Z":mZ,"Higgs":mH,"top":MT,"v":v,"nu2":mnu[1],"nu3":mnu[2],
       "N1":N[0],"N2":N[1],"N3":N[2],"meeting":MU0,"Planck":MP}
    out=[]; ks=list(S)
    for i,a in enumerate(ks):
        for j,b in enumerate(ks):
            if i<j: out.append((f"sqrt({a} x {b})",math.sqrt(S[a]*S[b])))
            if i!=j: out.append((f"{a}^2/{b}",S[a]**2/S[b]))
    return S,out
rng=np.random.default_rng(20261115)
for lab,MP in (("UNREDUCED Planck mass (main footing)",MPL),("reduced Planck mass",MRED)):
    S,fam=family(MP); lg=np.array([math.log10(x) for _,x in fam]); inr=(lg>2)&(lg<16)
    pr(f"  {lab}: {len(S)} scales, {len(fam)} balance points (geometric means and seesaws), {inr.sum()} of them between 1e2 and 1e16 GeV")
    t0=math.log10(ML0)
    for ns in (1,2):
        hits=[(n,x) for n,x in fam if abs(math.log10(x)-t0)<ns*sigL]
        T=rng.uniform(3,15,20000); frac=np.mean([(np.abs(lg-t)<ns*sigL).any() for t in T])
        pr(f"    within {ns} sigma of the closing mass (x/÷ {10**(ns*sigL):.2f}): {len(hits)} hit(s): " + "; ".join(f"{n} = {x:.2e}" for n,x in hits))
        pr(f"      luck: a target placed anywhere from 1e3 to 1e15 GeV has at least one such hit {100*frac:.0f}% of the time")
    best=sorted(fam,key=lambda p:abs(math.log10(p[1])-t0))[:6]
    pr("    six nearest: " + "; ".join(f"{n} = {x:.2e} ({math.log10(x)-t0:+.2f} dec)" for n,x in best))
pr("  Note: 'meeting' in the family is itself an output of the fit, so hits that use it are relations to test, not free predictions.")

# ------------------------------------------------------------------------------------------------------------------
head("5. Gravity: does the meeting know about G?")
aU=1/AU0
pr(f"Meeting: M = {MU0:.3e} GeV, alpha = 1/{AU0:.2f} = {aU:.5f}.")
for lab,MP in (("unreduced",MPL),("reduced",MRED)):
    pr(f" [{lab} Planck mass {MP:.3e} GeV]")
    pr(f"  (a) gravity's strength at the meeting, (M/M_Pl)^2 = {(MU0/MP)**2:.2e}; the forces' strength there = {aU:.2e}: gravity is {aU/(MU0/MP)**2:.1e} times weaker. FAIL")
    pr(f"  (b) 'gravity equals the others at the meeting' would put the grain at M/sqrt(alpha) = {MU0/math.sqrt(aU):.2e} GeV: {MP/(MU0/math.sqrt(aU)):.1e} times too low. FAIL")
    pr(f"  (c) energy where gravity would equal 1/{AU0:.0f}: {MP*math.sqrt(aU):.2e} GeV, {math.log10(MP*math.sqrt(aU)/MU0):.1f} decades above the meeting. FAIL")
    L=math.log(MP/MU0); pr(f"  (d) ln(grain/meeting) = {L:.2f}. As exp(-2 pi/(b alpha)) this needs b = {2*math.pi*AU0/L:.2f}. As (M_Pl/M)=alpha^-n: n = {L/math.log(AU0):.2f}.")
    pr(f"  (e) state counting (Sakharov 1967; Dvali 2007): grain^2 = N x meeting^2 needs N = {(MP/MU0)**2:.1e} states; TRPG has 124 (plus 30 far-side). FAIL")
gU=math.sqrt(4*math.pi*aU)
pr(f"  (f) string-type tie (Kaplunovsky 1988): meeting = 5.27e17 x g = {5.27e17*gU:.2e} GeV against {MU0:.2e}: {math.log10(5.27e17*gU/MU0):.1f} decades off. FAIL")
pr(f"  (g) proton from the grain (the earlier tie): ln(grain/proton) = {math.log(MPL/mp):.1f}; 1/alpha at the meeting = {AU0:.1f}; ratio {math.log(MPL/mp)/AU0:.3f}. No simple rule.")
pr("  (h) simple counts for 1/alpha at the meeting (38.97 if D at the meeting; 34.7 to 38.5 if D is lighter, see section 2):")
cands={"12 pi":12*math.pi,"27 + 12":39,"124/pi":124/math.pi,"4 pi^2":4*math.pi**2,"36":36,"40":40,"27 sqrt2":27*math.sqrt(2)}
for n,x in cands.items(): pr(f"      {n:14s} = {x:6.2f}  ({100*(x/AU0-1):+.1f}%)")
dAU=abs(mismatch(thr=thrL(ML0*10**sigL))[2]-AU0)+0.3
pr(f"      Error on 1/alpha at the meeting: about +-{dAU:.1f} (inputs, M_L range, three loops), and the D's position moves it by up to 4.8.")
pr("      Whole numbers alone are spaced 1 apart, so some 'simple count' always lands inside the band. NOT EVIDENCE.")
pr(f"  (i) the one pattern that lands: meeting = sqrt(M_L x M_Planck): off by {100*(math.sqrt(ML0*MPL)/MU0-1):+.0f}% (unreduced), {100*(math.sqrt(ML0*MRED)/MU0-1):+.0f}% (reduced). Tested in 4(v).")
fr=[1/4,1/3,1/2,2/3,3/4]; tol=max(abs(a_-0.5),da)
T=rng.uniform(0.15,0.85,20000)
p1=np.mean([min(abs(t-f) for f in fr)<tol for t in T])
pr(f"      luck: position {a_:.3f} +- {da:.3f}; a random position between 0.15 and 0.85 lies within {tol:.3f} of one of 1/4, 1/3, 1/2, 2/3, 3/4")
pr(f"      {100*p1:.0f}% of the time for one Planck convention, about {100*(1-(1-p1)**2):.0f}% allowing either convention; and it needs the D partners exactly at the meeting.")
pr("\nEnd.")
