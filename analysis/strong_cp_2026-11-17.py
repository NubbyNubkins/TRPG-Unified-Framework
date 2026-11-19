"""Borrowed test 6: the strong-CP angle (measured |theta_bar| < 1e-10, neutron EDM; Abel et al. 2020).
A. Nelson-Barr (Nelson 1984; Barr 1984) with TRPG's own D partners.
   Structure: CP exact, broken only by a complex scalar vev <S>; mass matrix [[m, 0],[B, mu]] with B = a<S> + b<S*> (D d^c S coupling),
   mu = D D^c mass, and the corner Q H D^c FORBIDDEN by a Z2 (D, D^c, S odd). Then arg det = 0: theta_bar = 0 at tree level.
   TEST FIXED BEFORE RUNNING (all three must pass):
   A1 fields: TRPG must already hold a complex scalar whose CP-breaking phase is physical (the Higgs pair: no, its phase is a gauge
      direction; the majoron field sigma: its phase is the majoron itself, removable by the lepton-number turn).
   A2 scale: with sigma used anyway (<sigma> = f/sqrt2, f = 6.3e10 GeV, Part 58) and D at F = 4.6e13 GeV (record), the CP phase
      reaching the quark tilts is ~ |B/mu|^2. PASS needs it to reach order one (measured CKM phase 1.14 rad) with couplings below sqrt(4 pi).
   A3 quality: the grain-suppressed operator Q H D^c sigma / M_grain (gravity is not expected to respect the Z2) must leave theta_bar < 1e-10.
B. Break check, fixed before running: the 'majoron = axion' route (KSVZ: Kim 1979, Shifman, Vainshtein & Zakharov 1980; majoron-axion:
   Langacker, Peccei & Yanagida 1986; SMASH: Ballesteros, Redondo, Ringwald & Tamarit 2017) uses only sigma, D, D^c, nu^c.
   B1 strengths: D then weighs y f/sqrt2 ~ 4.5e10 GeV (y = 1), not F. Do the force strengths still close (alpha_s within 2 sigma of 0.1180)?
   B2 quality: TRPG's keV majoron mass is itself an explicit breaking; it must not move theta_bar above 1e-10.
   B3 dark matter: can an axion with f/N (N = 3 D partners) be all the dark matter (Omega h^2 = 0.12)?
"""
import math, numpy as np
f=6.3e10; s=f/math.sqrt(2); F=4.6e13; MG=8.633e18; v=246.22; md=4.7e-3
print("A. NELSON-BARR")
print("  A1: Higgs phase = gauge direction (removable); majoron-field phase = the majoron (removable by the lepton-number turn).")
print("      -> no physical CP-breaking scalar in TRPG: Nelson-Barr needs a NEW complex scalar. FAIL under rule 6.")
ymax=math.sqrt(4*math.pi); r=2*ymax*s/F
rng=np.random.default_rng(1); best=0
for _ in range(20000):
    a=rng.uniform(-ymax,ymax,(3,3)); b=rng.uniform(-ymax,ymax,(3,3)); O1=np.linalg.qr(rng.normal(size=(3,3)))[0]; O2=np.linalg.qr(rng.normal(size=(3,3)))[0]; mu=F*O1@np.diag(rng.uniform(0.3,3,3))@O2   # D masses between 0.3F and 3F, as fixed
    ph=np.exp(1j*rng.uniform(0,2*np.pi)); B=a*s*ph+b*s*np.conj(ph)
    K=B.conj().T@np.linalg.inv(B@B.conj().T+mu@mu.T)@B
    best=max(best,np.abs(K.imag).max())
print(f"  A2: largest |B/mu| = {r:.1e}; largest imaginary part of the correction over 20,000 random choices: {best:.1e} (needs ~1). FAIL by {math.log10(1/best):.0f} powers of ten.")
eps=v/math.sqrt(2)*s/MG; th=r/2*eps/md
print(f"  A3: corner from the grain = {eps:.1e} GeV -> theta_bar ~ {th:.0e} (limit 1e-10) -> FAIL unless the Z2 and CP are exact against gravity (a further assumption).")
print()
print("B. MAJORON = AXION (break check)")
src=open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'leak_establish_three_loop_2026-11-16.py')).read()
ns={}; exec(src.split('print("A. grain')[0],ns); exec(src[src.index('nG=3; nt=1'):src.index('# sanity')],ns)
exec('''
def solveD(c, MD, **kw):
    def meet(ML):
        s_=predict_s2(thr=thrLD(ML,MD),**kw); return mismatch(s2=s_,thr=thrLD(ML,MD),**kw)[1]
    f_=lambda x: math.log10(c*meet(10**x)**2/MPL)-x
    M=10**brentq(f_,3,10.6,xtol=1e-6); return M,meet(M),predict_s2(thr=thrLD(M,MD),**kw)
def a3D(c,MD): return brentq(lambda a3: solveD(c,MD,a3=a3)[2]-S2,0.105,0.130,xtol=1e-6)
''',ns)
for MD in (4.5e10,):
    ML=10**ns['brentq'](lambda x: ns['mismatch'](thr=ns['thrLD'](10**x,MD))[0],4,6,xtol=1e-6); Fm=ns['mismatch'](thr=ns['thrLD'](ML,MD))[1]
    print(f"  B1: D at {MD:.1e} GeV: strengths close with the MMG at {ML:.2e} GeV (meeting {Fm:.2e}), but the gravity leak needs sqrt2 F^2/M_Planck = {math.sqrt(2)*Fm**2/ns['MPL']:.2e}:")
    print("      off by a factor ~1000. Scan of the strong strength: impose the gravity leak and solve, for each alpha_s:")
    for a3s in (0.105, 0.110, 0.115, 0.118, 0.120, 0.125, 0.130):
        try:
            Mx, Fx, s2x = ns['solveD'](math.sqrt(2), MD, a3=a3s)
            print(f"        alpha_s {a3s:.3f}: MMG {Mx:.2e}, meeting {Fx:.2e}, weak angle {s2x:.5f} vs measured {ns['S2']:.5f}")
        except ValueError:
            print(f"        alpha_s {a3s:.3f}: no solution with the gravity leak")
    print("      no strong strength between 0.105 and 0.130 rescues it -> FAIL (it would undo the MMG = gravity-leak lead).")
chi=0.0755**4; mmaj=1.4e-6
dth=mmaj**2*f**2/chi
print(f"  B2: keV majoron mass as explicit breaking shifts theta_bar by ~{dth:.0e} (limit 1e-10) -> the keV majoron and the axion cannot both stand.")
fa=f/3; ma=5.7e-6*(1e12/fa); om=0.12*(fa/9e11)**1.165*math.pi**2
print(f"  B3: f/N = {fa:.1e} GeV -> axion mass {ma*1e6:.0f} micro-eV; misalignment with the harmonic estimate at start angle pi: Omega h^2 ~ {om:.3f} = {om/0.12:.0%} of dark matter.")
# near the hilltop the cosine potential adds an anharmonic factor (Lyth 1992; Visinelli & Gondolo 2009; Kobayashi et al. 2013):
# Omega ~ theta^2 F(theta), F = [ln(e / (1 - theta^4/pi^4))]^(7/6)
Fan=lambda th: math.log(math.e/(1-(th/math.pi)**4))**(7/6)
need=0.12/(om/math.pi**2)
thn=__import__('scipy.optimize',fromlist=['brentq']).brentq(lambda th: th*th*Fan(th)-need, 1.0, math.pi*(1-1e-15))
print(f"      with the anharmonic factor the axion is all of the dark matter only for a start angle within {math.pi-thn:.1e} of pi "
      f"(theta = {thn:.6f}): not excluded, but it needs the start angle tuned to about {(math.pi-thn)/math.pi:.0e} of its range -> PARTIAL (tuned).")
print("      (N = 3 also leaves domain walls if the symmetry breaks after inflation.)")
