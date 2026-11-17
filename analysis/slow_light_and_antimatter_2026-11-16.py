# Checks: Higgs at the middle?, majoron wave size, light stretched since annihilation,
# where the "missing" matter number went (sphaleron split), confined light weighs.
import math
MPl=1.22089e28; grain=MPl/math.sqrt(2)   # eV
tip=1e-5
for nm,top in (("Planck",MPl),("grain",grain)):
    mid=math.sqrt(tip*top)
    print(f"middle height (tip 1e-5 eV, top {nm}): {mid/1e9:.0f} GeV")
    for p,m in (("Higgs value",246.22e9),("Higgs mass",125.2e9),("top",172.6e9)):
        print(f"   {p}: {math.log10(m/mid):+.2f} steps from the middle")
mid_nu=math.sqrt(3.6e-4*grain); print(f"middle if tip pinned at nu1 0.36 meV: {mid_nu/1e12:.2f} TeV")
# majoron de Broglie length in a galaxy
m=1.4e3; v=200/2.998e5; lam=1239.84e-9/(m*v)   # m (hc=1239.84 eV nm)
print(f"majoron 1.4 keV at 200 km/s: wavelength {lam*1e6:.2f} micrometres")
mfuzzy=1239.84e-9/(1*3.086e19*v)  # mass whose wavelength is ~1 kpc
print(f"mass needed for 1 kpc waves: {mfuzzy:.1e} eV")
# light stretched since proton-antiproton annihilation (~40 MeV temperature) to today 2.725 K
kT0=8.617e-5*2.725; print(f"stretch since annihilation (T~40 MeV): {40e6/kT0:.1e}")
print(f"typical photon then: ~{2.7*40e6/1e6:.0f} MeV gamma; now {2.7*kT0*1e3:.2f} meV (microwave, peak 160 GHz)")
Og=5.38e-5; Ob=0.049; print(f"energy in background light today / ordinary matter: {Og/Ob:.1e}")
# sphaleron split (Khlebnikov-Shaposhnikov 1988; Harvey-Turner 1990), SM with 3 families, 1 Higgs doublet
Nf,NH=3,1; cB=(8*Nf+4*NH)/(22*Nf+13*NH); cL=cB-1
print(f"B = {cB:.4f} (B-L) = 28/79; L = {cL:.4f} (B-L) = -51/79; L/B = {cL/cB:.3f}")
ne_per_B=0.88  # electrons per baryon (helium carries neutrons)
print(f"neutrino lepton number L_nu = L - n_e = {cL/cB - ne_per_B:.2f} per baryon (a surplus of about 2.7 antineutrinos per proton)")
nnu=336; nb=0.25e-6
print(f"background neutrinos per baryon: {nnu/nb:.1e}; surplus fraction {abs(cL/cB-ne_per_B)/(nnu/nb):.1e}")
# confined light weighs: proton mass mostly field energy; a mirror box with 1 J of light gains
print(f"1 joule of trapped light adds {1/2.998e8**2:.2e} kg")
