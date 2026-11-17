# Neutrinos near the centre; time both ways. Plain numbers only.
import math, itertools
print("1) Ladder openings vs heaviest member (is the neutrino the closest ladder to the centre line?)")
lad={"neutrinos":(37.13,0.0502),"charged leptons":(45.00,1.77686e9),"down quarks":(48.15,4.18e9),"up quarks":(52.147,172.57e9)}
for k,(o,m) in sorted(lad.items(),key=lambda x:x[1][0]): print(f"   {k:16s} opening {o:6.2f} deg   heaviest {m:.3g} eV")
print("   order of openings = order of heaviest masses: 4 of 4; chance of this by luck = 1/24 = 4.2% (hint, no law)")
print("\n2) If the neutrinos' smallness were a cancellation (pushes adding to almost zero), how fine must it be?")
m3=0.0502
for nm,mD in (("tau",1.77686e9),("muon",0.1057e9),("top",172.57e9)):
    print(f"   cancel a {nm}-size push down to 0.050 eV: precision 1 part in {mD/m3:.1e}")
print("   three separate cancellations at 1e-10..1e-12 with no reason: FAIL as the mechanism")
print("\n3) Dilution instead of cancellation (seesaw: Minkowski 1977, Gell-Mann-Ramond-Slansky 1979, Yanagida 1979, Mohapatra-Senjanovic 1980)")
print("   light weight = (Higgs link)^2 / (hidden partner far up the centre line). Partner weight needed per family:")
mnu=[0.36e-3,8.66e-3,50.2e-3]  # eV, rule-5 m1 and measured splittings (normal order)
for lab,mD in (("charged-lepton-size links",[0.511e6,105.66e6,1776.86e6]),("up-quark-size links",[2.16e6,1.27e9,172.57e9])):
    M=[d*d/n/1e9 for d,n in zip(mD,mnu)]
    print(f"   {lab:26s}: " + ", ".join(f"{x:.1e}" for x in M)+" GeV")
print("   (diagonal estimate, mixing ignored; the hidden-partner window 1e9..1e11 GeV was itself set this way, so landing in it is not new evidence)")
print("\n4) Time: the one part that does not reverse by itself")
J=3.08e-5; print(f"   J (area of the quark mixing triangle) = {J:.2e}; BaBar 2012 saw time-reversal broken directly at 14 sigma")
print("   T alone: FAIL (J != 0).  C, P, T all together (exact opposite side + antimatter): HOLDS to every test so far")
print("   K0 vs anti-K0 mass difference / mass < 6e-19 (PDG 2024): the reversed film weighs the same")
