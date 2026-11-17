# Top-bottom symmetry: neutrinos, hidden partners, MMG
import math
def kq(m):
    Q=sum(m)/sum(x**.5 for x in m)**2; e=math.sqrt(2*(3*Q-1)); return Q,math.degrees(math.atan(e/math.sqrt(2)))
d21,d31=7.42e-5,2.511e-3; m1=0.36e-3
nu=[m1,(m1*m1+d21)**.5,(m1*m1+d31)**.5]           # eV
lep=[0.51099895e6,105.6584e6,1776.86e6]            # eV
N=[l*l/n for l,n in zip(lep,nu)]                    # eV, diagonal seesaw estimate
print("1) Top three circles: hidden partners = neutrinos mirrored through the charged leptons (seesaw, diagonal)")
for i in range(3):
    print(f"   family {i+1}: neutrino {nu[i]*1e3:7.3f} meV  <- mirror at {lep[i]/1e6:8.3f} MeV ->  partner {N[i]/1e9:.3g} GeV ; steps of ten {math.log10(lep[i]/nu[i]):.2f} each way")
for nm,m in (("neutrinos (bottom)",nu),("charged leptons (middle, the pivots)",lep),("hidden partners (top)",N)):
    Q,o=kq(m); print(f"   {nm:38s} opening {o:6.2f} deg  (Q={Q:.4f})")
print("   spread of each ladder, heaviest/lightest:", f"bottom {nu[2]/nu[0]:.0f}x, middle {lep[2]/lep[0]:.0f}x, top {N[2]/N[0]:.0f}x")
print("   top spread = middle spread^2 / bottom spread:", f"{(lep[2]/lep[0])**2/(nu[2]/nu[0]):.0f}x  (mirror rule)")
print("\n2) Does the top ladder go back down to the bottom's shape?")
inv=[1/x for x in nu]; Q,o=kq(inv); print(f"   pure inversion of neutrinos (no pivot): opening {o:.2f} deg")
print("   with the lepton pivots the top ladder is much wider than the bottom one: the top is NOT a copy of the bottom")
print("\n3) MMG: the strengths fix only the average (geometric mean) of the three MMG weights, 2.2e8 GeV (two loop)")
gm=2.2e8
for nm,r in (("equal (cube pinches: top goes back together)",(1,1,1)),("split like the neutrinos",tuple(x/nu[0] for x in nu)),("split like the charged leptons",tuple(x/lep[0] for x in lep)),("split like the hidden partners",tuple(x/N[0] for x in N))):
    g=(r[0]*r[1]*r[2])**(1/3); w=[gm*x/g for x in r]
    print(f"   {nm:46s}: " + ", ".join(f"{x:.2g}" for x in w)+" GeV")
print("   all four give the same strengths at one loop: no measurement decides the MMG split (FREE)")
print("\n4) Neutrino <-> MMG link")
v=174.0; M=2.41e8
print("   shared: weak force (W, Z); both are colourless; the MMG neutral member is a heavy both-handed neutrino-like state")
ymin=math.sqrt(6.6e-25*16*math.pi/M); ymax=math.sqrt(0.05e-9*M)/v
print(f"   a direct link y (MMG-neutrino-Higgs) breaks kept matter number; decay before 1 s needs y > {ymin:.1e}")
print(f"   neutrino weights stay below 0.05 eV if y < {ymax:.1e}  (rough, Majorana-type estimate)")
print("   charged MMG -> neutral MMG + soft pion (weak-force split ~0.34 GeV, as for Higgsinos; Thomas & Wells 1998)")
