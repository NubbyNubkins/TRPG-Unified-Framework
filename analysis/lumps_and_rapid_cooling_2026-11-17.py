# Does rapid cooling make the lumps that would raise the early Hubble value? what in the theory could?
import numpy as np
hbar=6.582e-16; MPL=1.2209e28; T0=2.349e-4; G=6.674e-8; c=2.998e10; Msun=1.989e33; ly=9.461e17
def Hrad(T,g): return 1.66*np.sqrt(g)*T*T/MPL/hbar            # 1/s, T in eV
print("(1) How rapid is the cooling, compared with the thing that is happening? (process rate / cooling rate; big = slow, gentle cooling)")
rows=[("the surface's resistance switches on (160 GeV)",0.03*160e9/hbar,Hrad(160e9,106.75)),
      ("quarks lock into protons (155 MeV)",155e6/hbar,Hrad(155e6,17.25)),
      ("neutrinos let go (1.5 MeV)",None,None),
      ("first nuclei (0.066 MeV)",None,None)]
for n,r,h in rows[:2]: print(f"   {n:48s} {r/h:.1e}")
# deuterium formation rate p+n -> D + light: <sigma v> ~ 4.6e-20 cm^3/s ; n_b at T
nb=lambda T: 6.1e-10*0.2436*(T/1.9733e-5)**3      # cm^-3, T in eV (n_gamma = 0.2436 T^3)
print(f"   {'neutrinos let go (1.5 MeV)':48s} 1 (by definition: this is why they let go)")
print(f"   {'first nuclei (0.066 MeV)':48s} {4.6e-20*nb(0.066e6)/Hrad(0.066e6,3.36):.1e}")
# atoms: recombination at stretch 1100, x_e ~ 0.1, with the bottleneck factor ~0.03
z=1100; T=2.7255*(1+z); al=4.309e-13*(T/1e4)**-0.6166/(1+0.6703*(T/1e4)**0.53); nH=1.9e-7*(1+z)**3; H=2.184e-18*np.sqrt(0.315*(1+z)**3*(1+(1+z)/3400))
print(f"   {'atoms form (0.26 eV)':48s} {al*0.1*nH/H:.0f} before the bottleneck, about {0.03*al*0.1*nH/H:.0f} with it; falls to 1 as the last electrons are left over")
print("   -> atom formation is the ONE step where the cooling is rapid compared with the process. That is why its timing depends on how lumpy matter is.")
print("\n(2) Atoms forming is what lets matter go lumpy (pressure drops): smallest lump that can pull itself together")
rho_m=0.143*1.878e-29*1101**3; rho_b=0.02237/0.143*rho_m
for lab,cs in (("just before (matter tied to the light)",c/np.sqrt(3*(1+0.63))),("just after (atoms, light set free)",np.sqrt(5*8.617e-5*3000/(3*938.27e6))*c)):
    lam=cs*np.sqrt(np.pi/(G*rho_m)); M=np.pi/6*rho_b*lam**3/Msun; print(f"   {lab:42s} sound speed {cs/1e5:,.0f} km/s; smallest lump {M:.1e} suns")
print("   -> true, and large (a factor of 1e13), but these lumps come AFTER atoms form, so they cannot change WHEN atoms form.")
print("\n(3) When could lumps have been made? Test: lumps already there at 4 minutes would spoil the deuterium (D/H goes as density^-1.6)")
for lab,z3 in (("mild 0.6/1/1.4",(0.6,1,1.4)),("medium 0.3/1/1.7",(0.3,1,1.7)),("strong 0.1/1/1.9",(0.1,1,1.9))):
    print(f"   {lab:20s} deuterium x{np.mean([d**-1.6 for d in z3]):.2f}   (measured to 1%: Cooke, Pettini, Steidel 2018)")
print("   -> lumps must be made AFTER the first nuclei (4 minutes) and BEFORE atoms (363 thousand years): light between 0.07 MeV and 0.3 eV.")
print("\n(4) Size of one horizon then, stretched to today (a lump cannot be made bigger than this at that time)")
def comov(T,g): t=2.42/np.sqrt(g)/(T/1e6)**2; return 2*c*t*(T/T0)*(g/3.91)**(1/3)/ly
for lab,T,g in (("the light began (1e6 GeV)",1e15,106.75),("the light began (2.4e8 GeV)",2.4e17,106.75),("quarks lock (155 MeV)",155e6,17.25),("first nuclei (0.066 MeV)",0.066e6,3.36),("majoron level (8.8 keV)",8.8e3,3.36),("light at 1 keV",1e3,3.36),("light at 10 eV",10,3.36)):
    print(f"   {lab:30s} {comov(T,g):.2e} light-years today")
print("   needed size of the lumps (Jedamzik & Pogosian): roughly 3 to 3,000 light-years today")
