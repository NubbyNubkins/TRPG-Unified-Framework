# Opposites, pulls and inputs: numbers
import math
def kq(m):
    Q=sum(m)/sum(x**.5 for x in m)**2; e=math.sqrt(2*(3*Q-1)); return math.degrees(math.atan(e/math.sqrt(2)))
d21,d31=7.42e-5,2.511e-3
nuop=lambda m1: kq([m1,(m1*m1+d21)**.5,(m1*m1+d31)**.5])
print("1) Grain squish passes energy to the MMG (counter relation): MMG = F^2/grain")
F=4.56e13
for g in (8.63e18/2,8.63e18,2*8.63e18): print(f"   grain {g:.2e} -> MMG {F*F/g:.2e} GeV  (heavier grain, lighter MMG)")
print("\n2) Opposite press on the leptons? Lepton opening = 45 deg to within the twist's 7.9e-7;")
print("   D's press changes the down ladder at the 1e-3 level and more (1-kappa = 8.4e-4): an equal opposite press on leptons is ruled out by ~1000x")
print("\n3) Up and down on opposite sides of 45?")
op={"neutrinos (rule 5)":nuop(0.36e-3),"charged leptons":45.0,"down":48.150,"up":52.147}
for k,v in op.items(): print(f"   {k:20s} {v:7.3f}  ({v-45:+.3f} from 45)")
print("   up and down are on the SAME side (+7.15, +3.15): FAIL for up/down as opposites")
target=90-52.147; lo,hi=0.0,2e-3
for _ in range(80):
    mid=(lo+hi)/2; lo,hi=(mid,hi) if nuop(mid)>target else (lo,mid)
print(f"   up (+7.147) and neutrinos ({nuop(0.36e-3)-45:+.3f}) ARE on opposite sides; exact mirror (37.853) needs m1 = {lo*1e3:.3f} meV vs rule 5 0.36: HINT, conflicts with rule 5 by 0.72 deg")
print("\n4) Gravity's pull at each rung: (weight / Planck)^2  (the other forces are ~1/137 to ~1)")
MPl=1.22089e19
for nm,m in (("electron",0.000511),("proton",0.938272),("top",172.57),("MMG",2.41e8),("D",4.56e13),("grain 2",8.63e18)):
    print(f"   {nm:9s} {m:10.3g} GeV  gravity strength {(m/MPl)**2:.2e}")
print("\n5) Force strengths along the height (one loop, with three MMGs at 2.6e7 GeV, the one-loop closing weight)")
MZ=91.1876; a_em=1/127.951; s2=0.23122; a3=0.1180
A=[(3/5)*(1-s2)/a_em,s2/a_em,1/a3]; b=[41/10,-19/6,-7]; db=[1.2,2.0,0.0]; M=2.63e7
def inv(mu):
    t=math.log(mu/MZ); out=[A[i]-b[i]/(2*math.pi)*t for i in range(3)]
    if mu>M: out=[out[i]-db[i]/(2*math.pi)*math.log(mu/M) for i in range(3)]
    return out
for mu in (91.19,1e3,1e6,2.63e7,1e10,1e12,5.07e13):
    v=inv(mu); print(f"   {mu:9.3g} GeV: 1/strength  hypercharge {v[0]:6.2f}  weak {v[1]:6.2f}  strong {v[2]:6.2f}")
print("\n6) Where the antimatter went: eta = 6.1e-10 baryons per photon ->", f"{1/6.1e-10:.2e} photons per surviving baryon (Planck 2018; BBN)")
print("   early bag: about 1,640,000,001 quarks for every 1,640,000,000 antiquarks; the pairs became light (today's 411 photons per cm^3)")
print("\n7) If all matter met equal antimatter today: ordinary matter 4.9% of the bag -> light + neutrinos;")
print("   dark matter 26.5% (majoron: its own antiparticle) and dark energy 68.6% untouched (Planck 2018). Total energy unchanged.")
