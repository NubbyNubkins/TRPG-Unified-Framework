# (A) tie the handover number to the start level; (B) quark mirror hint at each quark's own energy
import numpy as np
ln=np.log
print("(A) Handover number X = (share of the start's energy through the partner) x (start level / weight of one piece of the source)")
X=1.68e-3; M2=1.63e10; Ttop=2.41e8; Tbot=1e6     # small-tilt corner; start band
msrc_min=2*M2
print(f"   needed X = {X:.2e}; source piece must outweigh two partners: > {msrc_min:.2e} GeV; start level < {Ttop:.2e} GeV (MMG never made)")
print(f"   old ceiling on X: 0.05 -> room 30x.  New ceiling: {Ttop/msrc_min:.2e} -> room {Ttop/msrc_min/X:.1f}x")
print(f"   share of the start's energy that must pass through partner 2: at least {X*msrc_min/Ttop:.2f}")
print(f"   start level: at least {X*msrc_min:.2e} GeV (share 1) -> band {X*msrc_min:.1e} to {Ttop:.2e} GeV (was {Tbot:.0e} to {Ttop:.2e}): {Ttop/(X*msrc_min):.1f}x wide instead of {Ttop/Tbot:.0f}x")
print(f"   weight of one piece of the source: {msrc_min:.2e} to {Ttop/X:.2e} GeV")
g=106.75; MP=1.22e19
age=lambda T: 0.301*MP/(np.sqrt(g)*T*T)*6.58e-25
print(f"   ages: {age(Ttop):.1e} s to {age(X*msrc_min):.1e} s")
print("   things in the main line inside the source's weight range: settling value 6.3e10, partner 3 at 1.01e11")
ms=3.16*M2; br=0.093/1.093
print(f"   settling piece as the source (weighs {ms:.2e}, at most {br:.3f} to partners): X at most {br*Ttop/ms:.1e}: short by {X/(br*Ttop/ms):.1f}x -> FAIL (and 91% hot majorons)")
print(f"   D (4.6e13) as the source: X at most {Ttop/4.6e13:.1e}: short by {X/(Ttop/4.6e13):.0f}x -> FAIL")
for th in (0.5,1,1.81,3.14):
    print(f"   majoron starts to swing at {1.45e6/th**2:.1e} GeV (angle {th}): below the new band's floor -> amount stays as an earlier test (no new tie)")
print(f"   the swing would begin inside the band only for angle < {np.sqrt(1.45e6/(X*msrc_min)):.2f}, i.e. majoron weight > {8.8/(1.45e6/(X*msrc_min))**2:.0f} keV: above the 52 keV ceiling -> excluded")

print("\n(B) Quark mirror hint: middle-rung fractions, up + down = 1 ?")
f=lambda m: ln(m[1]/m[0])/ln(m[2]/m[0])
cases={"as weighed, top by its pole weight":([2.16,1273,172570],[4.70,93.5,4183]),
       "each at its own energy (top 162.7 GeV at its own level)":([2.16,1273,162700],[4.70,93.5,4183]),
       "light three at 1 GeV instead of 2 GeV (x1.35), heavy at own":([2.16*1.35,1273,162700],[4.70*1.35,93.5*1.35,4183]),}
for k,(u,d) in cases.items(): print(f"   {k:62s} up {f(u):.4f} down {f(d):.4f} sum {f(u)+f(d):.4f}")
rng=np.random.default_rng(3); S=[]
for _ in range(200000):
    u=[rng.normal(2.16,0.21),rng.normal(1273,4.6),rng.normal(162700,700)]; d=[rng.normal(4.70,0.07),rng.normal(93.5,0.8),rng.normal(4183,7)]
    S.append(f(u)+f(d))
S=np.array(S); print(f"   own-energy sum with errors: {S.mean():.4f} +- {S.std():.4f}: {(S.mean()-1)/S.std():+.1f} sigma from 1")
