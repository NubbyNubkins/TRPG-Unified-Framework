# Antimatter on the top half? numbers used
import math, itertools
print("1) Is antimatter lighter? Lab tests (matter vs antimatter weight)")
for k,v in (("proton vs antiproton, charge-to-mass (BASE 2022)","equal to 1.6e-11"),
            ("hydrogen vs antihydrogen 1S-2S light (ALPHA 2018)","equal to 2e-12"),
            ("neutral kaon vs antikaon mass (PDG 2024)","equal to 6e-19"),
            ("electron vs positron mass (PDG 2024)","equal to 8e-9"),
            ("antihydrogen falls (ALPHA-g 2023)","0.75 +- 0.21 g, down, not up")):
    print(f"   {k:52s}: {v}")
print("   => antimatter lighter: FAIL (any lightness would have to be below these levels)")
print("\n2) Equal antimatter hidden nearby? Annihilation glow at matter/antimatter borders not seen;")
print("   domains must be at least about the size of the visible universe (Cohen, De Rujula & Glashow 1998); no antihelium confirmed (AMS-02)")
print("\n3) Equal and opposite matter number stored in weakly-linked top-half particles (Dirac leptogenesis, Dick, Lindner, Ratz & Wright 2000)")
print("   needs neutrinos WITHOUT their own heavy (Majorana) weight; main line has the seesaw planks -> conflicts with main line")
d21,d31=7.42e-5,2.511e-3; m1=0.36e-3; m=[m1,(m1*m1+d21)**.5,(m1*m1+d31)**.5]; U=[0.681,0.297,0.022]
t=[u*x for u,x in zip(U,m)]; mx=sum(t); mn=max(0,max(t)-(sum(t)-max(t)))
print(f"   deciding test: neutrinoless double beta decay. Main line (Majorana): m_bb = {mn*1e3:.1f} to {mx*1e3:.1f} meV; Dirac picture: exactly 0")
print("   next experiments reach about 10-20 meV (LEGEND-1000, nEXO): cannot decide before ~2035+")
print("\n4) Energies made of matter-antimatter pairings (known physics, HOLDS):")
for k,v in (("pion (quark + antiquark)","140 MeV; every meson is a matter-antimatter pair"),
            ("positronium para (e- e+, spins opposite)","lives 1.25e-10 s"),
            ("positronium ortho (spins same)","lives 1.42e-7 s: 1,100x longer, a 'stage of stability'"),
            ("neutral kaon","switches matter<->antimatter 5.6e9 times a second"),
            ("proton sea","quark-antiquark pairs carry a few % of its momentum"),
            ("force strengths running","made entirely of virtual particle-antiparticle pairs; MMG pairs join above 2.4e8 GeV")):
    print(f"   {k:44s}: {v}")
J=3.08e-5; print(f"\n5) Matter side vs antimatter side differ only by the mixing-triangle area J = {J:.2e} (Jarlskog 1985); weights identical")
