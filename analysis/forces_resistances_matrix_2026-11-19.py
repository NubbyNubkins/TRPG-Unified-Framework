"""Forces, resistances, matrix rebuild, universe analogies (19 Nov 2026). Exploration only."""
import math, numpy as np
pi=math.pi; rng=np.random.default_rng(1)
G=6.674e-11; c=2.998e8; Msun=1.989e30; Mpc=3.086e22; pc=3.086e16
print("="*72); print("A. RELATIONSHIP MATRIX: matter on the diagonal, forces on the links")
# Sites = places in the grid; each site holds a colour triplet (matter). Link U_ij in SU(3) connects sites.
def su3():
    z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)); q,r=np.linalg.qr(z); q=q@np.diag(np.diag(r)/abs(np.diag(r)))
    return q/np.linalg.det(q)**(1/3)
N=4; psi=[rng.normal(size=3)+1j*rng.normal(size=3) for _ in range(N)]
U={(i,(i+1)%N):su3() for i in range(N)}
def energy(psi,U):  # sum of 'held-to-linked-to-held' terms + closed loops of links (the force's own energy)
    hop=sum(np.vdot(psi[i],U[(i,j)]@psi[j]) for (i,j) in U)
    loop=np.trace(U[(0,1)]@U[(1,2)]@U[(2,3)]@U[(3,0)])
    return hop,loop
h0,l0=energy(psi,U)
g=[su3() for _ in range(N)]   # change each site's private 'colour axes' independently
psi2=[g[i]@psi[i] for i in range(N)]; U2={(i,j):g[i]@U[(i,j)]@g[j].conj().T for (i,j) in U}
h1,l1=energy(psi2,U2)
print(f"  held-link-held energy before/after relabelling every spot's colours: {h0:.6f} / {h1:.6f}")
print(f"  closed loop of links (force field energy) before/after:              {l0:.6f} / {l1:.6f}")
print("  -> identical: matter on the diagonal + forces on the links is exactly self-consistent (Wilson 1974).")
print("  Layout of one generation's relationship matrix (diagonal = held states, off-diagonal = links):")
print("    diagonal: Q(3 colours x 2), u_R(3), d_R(3), L(2), e_R(1), nu_R(1) = 16 held states")
print("    links between DIFFERENT spots: 8 gluon + 3 weak + 1 hypercharge = 12 link types (the pulls)")
print("    links between the TWO SIDES of the SAME spot (left <-> right): the Higgs (4 real parts) = the resistance")
print("    the spacing of the grid itself: gravity / geometry (not a link at all)")
print("    link strengths across the Higgs link: the 13 masses+mixings (where all TRPG's open numbers live)")

print("="*72); print("B. HOW MANY FORCES DO THE NUMBERS NEED?")
aem=1/127.951; s2=0.23122; a3=0.1180; MZ=91.1876
a=[(5/3)*aem/(1-s2), aem/s2, a3]; b=[41/10,-19/6,-7.]
inv=lambda i,mu: 1/a[i]-b[i]/(2*pi)*math.log(mu/MZ)
for mu in (1e3,1e13,1e16):
    print(f"  1/strength at {mu:.0e} GeV: {inv(0,mu):5.1f} {inv(1,mu):5.1f} {inv(2,mu):5.1f}")
t12=2*pi*(1/a[0]-1/a[1])/(b[0]-b[1]); M12=MZ*math.exp(t12)
print(f"  pulls 1 and 2 meet at {M12:.1e} GeV where pull 3 is at 1/{inv(2,M12):.1f} vs 1/{inv(0,M12):.1f}: SM alone = 3 separate pulls (miss {abs(inv(2,M12)-inv(0,M12))/inv(0,M12)*100:.0f}%)")
print("  TRPG main line (with far-side MMG, fitted choice 2): the three meet near 5e13 GeV -> ONE pull, alpha_s 0.1178 vs 0.1180+/-0.0009")
# Strong 'resistance' derived from the pull: dimensional transmutation
b0=11-2*5/3; L=MZ*math.exp(-2*pi/(b0*a3))
print(f"  Strong-force resistance scale from the pull alone (one loop): Lambda_QCD ~ {L*1000:.0f} MeV (higher loops ~210 MeV); the proton's 938 MeV follows from it (lattice QCD)")
# Weak force = EM-family pull + Higgs resistance (Yukawa range)
mW=80.37; hbarc=0.1973e-15  # GeV*m
print(f"  Weak force range = hbar/(m_W c) = {hbarc/mW:.2e} m: the photon-family pull, cut short by the Higgs resistance (Anderson 1963, superconductor)")
# Gravity vs EM: can gravity be the same pull?
mp=1.6726e-27; e=1.602e-19; k=8.988e9
print(f"  gravity / electric between two protons = {G*mp**2/(k*e**2):.2e}")
print("  MICROSCOPE: everything falls the same to 1 part in 10^15; electric pull differs 100% (neutron 0, proton 1)")
print("  -> gravity is NOT the same pull as the other three (it acts on all energy, only attracts, spin 2).")
# Expansion: pull at the edge vs tension everywhere
print("  Expansion as a pull from the EDGE: needs a centre. CMB is the same in every direction to 1 part in 100,000;")
print("     rotation/anisotropic expansion disfavoured 121,000:1 (Saadeh 2016). FAILS.")
print("  Expansion as a pull EVERYWHERE (tension): dark energy has pressure = -energy density (w = -1.03+/-0.03):")
print("     negative pressure IS tension, a pull, like a stretched rubber band. In gravity, tension drives FASTER expansion")
print("     (acceleration ~ -(rho + 3p) = +2 rho_L > 0). HOLDS.")
Om,OL=0.315,0.685
z_acc=(2*OL/Om)**(1/3)-1; z_eq=(OL/Om)**(1/3)-1
print(f"  Cosmic pull (tension, dark energy {OL}) vs resistance (matter's gravity {Om}): equal at z = {z_eq:.2f}; speed-up begins z = {z_acc:.2f} (measured ~0.6-0.7)")

print("="*72); print("C. RESISTANCES: list and try to merge")
Mpl=1.221e19; v=246.22; rhoL4=2.24e-12  # GeV (2.24 meV)
scales=[("Planck (space's stiffness, gravity)",Mpl),("far side MMG (TRPG)",2.4e8),("Higgs (weak/mass resistance)",v),
        ("QCD vacuum (strong resistance)",0.21),("neutrino weight m3",0.050e-9),("dark energy (rho_L^1/4)",rhoL4)]
for n,s in scales: print(f"  {n:38s} {s:9.2e} GeV")
print("  Merge tests:")
print(f"   QCD resistance from the strong pull: HOLDS (Lambda from alpha_s, above). One fewer independent resistance.")
print(f"   Higgs = dark energy?  Higgs field energy ~ (100 GeV)^4 vs ({rhoL4*1e12:.2f} meV)^4: off by 10^{4*math.log10(100/rhoL4):.0f}. FAILS.")
print(f"   Zero-point 'vacuum' energy up to Planck vs dark energy: off by 10^{4*math.log10(Mpl/rhoL4):.0f} (cosmological constant problem). FAILS (unsolved for everyone).")
m=(0.36e-3,8.6e-3,50.2e-3)
print(f"   Dark energy ~ neutrino weights? rho_L^1/4 = 2.24 meV vs m1 {m[0]*1e3:.2f}, m2 {m[1]*1e3:.1f}, m3 {m[2]*1e3:.1f} meV; m2^2/m3 = {m[1]**2/m[2]*1e3:.2f} meV. Same ballpark, no exact rule: HINT only.")
print(f"   (TeV)^2/M_Planck = {1e3**2/Mpl*1e12:.2f} meV vs 2.24 meV: off ~{2.24/(1e3**2/Mpl*1e12):.0f}x (needs ~5.2 TeV instead of 1 TeV, nothing there). FAILS as a rule.")
H0=67.4e3/Mpc; a0=1.2e-10
print(f"   Dark-matter-like 'resistance' tied to expansion: Milgrom a0 = 1.2e-10 m/s^2 vs c*H0/2pi = {c*H0/(2*pi):.2e}. Coincidence HOLDS (Milgrom 1983);")
print("     resistance alone FAILS the Bullet Cluster and the CMB third peak, so some extra matter/wave is still needed (TRPG: majoron, wave DM).")
print("   Higgs shape from space's stiffness: TRPG's grain rule (Shaposhnikov-Wetterich) gives Higgs 127.3 vs 125.13 (1.4 sigma): PENDING, links 2 resistances.")
print("   Dark energy + dark matter as one fluid (Chaplygin gas): FAILS structure growth (Sandvik et al. 2004) in its simple form.")

print("="*72); print("D. UNIVERSE AS SOLAR SYSTEM / RUBBER BAND")
Lam=3*OL*H0**2/c**2
print("  Newtonian form of the expansion (Milne-McCrea 1934): acceleration = -G M / r^2  +  (Lambda c^2 / 3) r")
print("     first term = a sun's pull, fading with distance; second = a spring (Hooke) force growing with distance, pointing OUT.")
for name,M in [("Sun",1*Msun),("Milky Way",1.5e12*Msun),("Local Group",3e12*Msun),("Virgo cluster",1.2e15*Msun)]:
    R=(3*G*M/(Lam*c**2))**(1/3)
    print(f"   {name:13s}: pull = outward spring at {R/Mpc:7.3f} Mpc ({R/pc:,.0f} pc)")
print("   Local Group measured zero-velocity edge ~0.9-1.0 Mpc, zero-gravity ~1.2-1.6 Mpc (Chernin et al.): HOLDS.")
print("   Sun's balance radius ~100 pc vs solar system edge (Oort cloud) ~1 pc: the solar system is held 100x inside its balance point.")
sigma=0.18  # GeV^2 string tension
print(f"  Quark rubber band: constant pull sigma = {sigma} GeV^2 = {sigma/0.1973:.2f} GeV/fm (~14 tonnes); energy = sigma x length.")
print("  Dark energy: constant energy per VOLUME; stretching adds energy in proportion = same constant-tension rule in 3D. HOLDS (structure).")
print("  Differences: band pulls its ends together; cosmic tension speeds expansion (gravity of tension is repulsive). Band force is constant;")
print("     cosmic 'spring' grows with distance. Band snaps into new quark pairs; cosmic stretching makes particles only negligibly today (Parker 1968).")
print("  'Expansion shows decays': supernova Ni-56 light curves at z~1 run slow by exactly (1+z) (Goldhaber 2001; Blondin 2008): local decay rates")
print("     unchanged by expansion. Oklo reactor: constants fixed to ~1e-8 over 2 billion years. Expansion causing decays FAILS.")
print("  One central sun for the universe: FAILS (no centre, see B). Every galaxy sees the same expansion around itself.")
