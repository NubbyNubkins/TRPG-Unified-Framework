"""Five further tests of open questions (18 November 2026)."""
import numpy as np
print("1. PROTON HOME = REDUCED MASS")
mp=938.272; me=0.51099895
home=mp/(1+np.exp(-2*0.5*np.log(me/mp))); red=me*mp/(me+mp)
print(f"   electron: logistic home {home:.6f} MeV, reduced mass {red:.6f} MeV (identical by algebra)")
print("   link value (m1+m2)/sqrt(m1 m2) = 2cosh S: a function of the mass ratio only;")
print("   the 171-pair matrix of 19 particles holds exactly 18 independent numbers (the mass ratios).")

print("\n2. WHAT PUSHES QUARKS OFF THE 45 DEG LEPTON CONE?")
op={'lepton':45.000,'down':48.150,'up':52.147}; err={'down':0.015,'up':0.018}
dd=op['down']-45; du=op['up']-45
r=du/dd; sr=r*np.hypot(err['down']/dd,err['up']/du)
print(f"   shift down {dd:.3f} deg, up {du:.3f} deg, ratio {r:.3f} +- {sr:.3f}")
for name,pred in [("|charge| (1/3 : 2/3)",2.0),("charge^2",4.0)]:
    print(f"   shift proportional to {name}: needs {pred}, measured {r:.3f} -> {(r-pred)/sr:+.0f} sigma -> FAIL")
Q=lambda th: (1+np.tan(np.radians(th))**2)/3
qd,qu=Q(op['down'])-2/3,Q(op['up'])-2/3
print(f"   in Koide-Q terms: excess down {qd:.4f}, up {qu:.4f}, ratio {qu/qd:.2f} (neither 2 nor 4)")
print("   The D fold presses DOWN quarks only, yet UP quarks sit farther from 45 deg:")
print("   so the fold is not what moves quarks off the lepton cone. FAIL for that reading.")

print("\n3. DOES THE FAR-SIDE CHAIN CONTINUE DOWN TO THE HIGGS VALUE?")
MPl=1.22089e19; grain=MPl/np.sqrt(2); F=4.62e13; MMG=F**2/grain; step=grain/F
nxt=MMG/step
print(f"   grain {grain:.3e}, F {F:.3e}, MMG = F^2/grain {MMG:.3e} GeV, step {step:.3e}")
print(f"   next rung MMG/step = {nxt:.0f} GeV vs Higgs value 246.2 ({nxt/246.22:.1f}x), top 172.6 ({nxt/172.57:.1f}x) -> FAIL")

print("\n4. MATTER PARITY (-1)^{3(B-L)} vs MMG DECAY")
reps={'16 quarks (B-L=1/3)':1/3,'16 leptons (B-L=-1)':-1,'10 weak pairs, MMG (B-L=0)':0,'10 coloured D (B-L=+-2/3)':2/3,'1 singlet (B-L=0)':0}
for k,v in reps.items(): print(f"   {k:30s} parity {'odd' if round(3*v)%2 else 'even'}")
print("   A fermion that is EVEN cannot decay to an odd number of SM fermions (all odd) plus bosons:")
print("   with the parity exact, the lightest even fermion (neutral MMG) is STABLE, and MMG-lepton or D-down")
print("   mixing (needed for MMG decay and for the D fold) is FORBIDDEN. The exact parity and those two pieces contradict.")
# stable MMG as dark matter: direct detection (Dirac doublet, Z exchange) rough
GF=1.1664e-5; mn=0.9396; s2=0.2312; Z,A=54,131; N=A-Z
mu=mn  # m_chi >> m_n
gV=0.5  # neutral component of a vector-like doublet (T3=-1/2 both chiralities), coupling to Z ~ g/cw * 1/2 * 2 chiralities /2
sigma_n=GF**2*mu**2/(2*np.pi)*((N-(1-4*s2)*Z)/A)**2*(2*gV)**2*0.389e-27  # GeV^-2 -> cm^2
lim_1TeV=1e-46; lim=lim_1TeV*(MMG/1e3)
print(f"   stable neutral MMG as dark matter: sigma_n ~ {sigma_n:.1e} cm^2 vs direct-detection limit ~ {lim:.0e} at {MMG:.1e} GeV")
print(f"   -> excluded by ~{sigma_n/lim:.0f}x (order of magnitude). Mass splitting from v^2/F = {246.22**2/F*1e9:.1f} eV, far below ~100 keV: stays Dirac.")

print("\n5. MATTER EXCESS WITHOUT RULE 5 (heavy-partner masses free)")
v=174.0; m3=0.0503e-9  # GeV
for M1 in (1e9,3e9,1e10,3e10,1e11):
    eps=3*M1*m3/(16*np.pi*v**2)
    for kappa in (0.01,0.1):
        eta=0.96e-2*eps*kappa
    e01=0.96e-2*eps*0.01; e1=0.96e-2*eps*0.1
    print(f"   M1 {M1:.0e} GeV: max asymmetry eps {eps:.1e}; eta {e01:.1e} (kappa 0.01) to {e1:.1e} (kappa 0.1); need 6.1e-10")
print("   -> allowed for M1 above ~1e9-1e10 GeV with a start hotter than that (Davidson & Ibarra 2002). Allowed, not predicted.")
