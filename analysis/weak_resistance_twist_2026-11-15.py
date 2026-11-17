"""Numbers: the weak force as a change of perspective against the Higgs field's resistance; motion of the needles; the twist.
Credit: Fermi 1934; Glashow 1961; Weinberg 1967; Salam 1968; Higgs 1964; Englert & Brout 1964; Gell-Mann & Nishijima (Q = T3 + Y);
Pontecorvo 1957; Gell-Mann & Pais 1955; Christenson, Cronin, Fitch & Turlay 1964; Kobayashi & Maskawa 1973; ATLAS & CMS 2022; PDG 2024."""
import math
def pr(*a): print(*a, flush=True)
v = 246.22; MZ = 91.1876; MW = 80.369; s2 = 0.23122; aem = 1/127.95
g = math.sqrt(4*math.pi*aem/s2); gp = g*math.sqrt(s2/(1 - s2))
pr("1. Weak force against the Higgs field's resistance")
pr(f"  W mass = g v / 2 = {g*v/2:.2f} GeV (measured {MW}); Z mass = W / cos(theta_W) = {g*v/2/math.sqrt(1-s2):.2f} GeV (measured {MZ}); photon 0 (unresisted mix)")
pr(f"  Reach of the weak force = hbar c / M_W = {0.1973269804/MW:.5f} fm (a proton is ~0.84 fm across)")
GF = math.sqrt(2)*g**2/(8*MW**2); pr(f"  Fermi strength = sqrt2 g^2 / (8 M_W^2) = {GF:.4e} GeV^-2 (measured 1.1664e-5): the heavier the resistance, the weaker the force at low energy")
pr(f"  Beta decay has ~1 MeV to spend but the W needs {MW*1e3:.0f} MeV: suppression ~ (1 MeV / M_W)^2 = {(1e-3/MW)**2:.1e}")
pr(f"  rho = M_W^2 / (M_Z^2 cos^2) = {MW**2/(MZ**2*(1-s2)):.4f} (tree level; loop effects from the top's mass shift it by ~1%; TRPG main line: rho = 1 at tree level)")
pr("\n2. Motion of the needles (measured)")
hbar_eVs = 6.582119569e-16
for name, dm in (("K0 (s-d)", 3.484e-6), ("B0 (b-d)", 3.334e-4), ("Bs (b-s)", 1.1693e-2)):
    pr(f"  {name:10s} swing between matter and antimatter: {dm/hbar_eVs/1e12:.4f} per ps (mass split {dm:.3e} eV)")
for name, dm2, E in (("atmospheric nu", 2.5e-3, 1.0), ("solar nu", 7.4e-5, 0.005)):
    L = 2.48*E/dm2   # km, E in GeV, dm2 in eV^2 (full oscillation length)
    pr(f"  {name:15s} full swing length {L:,.0f} km at {E} GeV")
pr("  Higgs vibration: the Higgs boson (125.2 GeV) is a vibration of the field that sets every axis length; each particle's coupling to it is")
pr("  proportional to its mass (measured from muon to top over a factor 1600: ATLAS & CMS, Nature 2022). All axes vibrate together, in proportion.")
pr("\n3. The twist and polarity")
pr("  Lepton twist 2/9: electron, muon and tau all have charge -1, so a charge-sign reading cannot tell them apart. Positron and electron masses")
pr("  agree to about 8 parts in a billion (PDG 2024), so the twist is the same for antimatter: it is not about polarity.")
pr("  Quark turn's twist (phase ~66 deg, J = 3.0e-5): this IS a matter/antimatter difference (CP violation: Christenson, Cronin, Fitch & Turlay")
pr("  1964; Kobayashi & Maskawa 1973). Known size is ~1e10 too small to make the matter excess by itself (Gavela et al. 1994; Huet & Sather 1995).")
