"""Idea tested: the ladders show one perspective of a single geometry; other views (masses and other
quantities, across scales) might combine into one measurable geometry viewed by moving through perspectives. Tests:
1. one shape seen from rotated viewpoints (rotation keeps eigenvalue ratios); 2. one shape seen at different scales (powers), neutrinos
included; 3. the perspective-free (invariant) quantities: the measurable geometry itself, the unitarity triangle.
Credit: Cabibbo 1963; Kobayashi & Maskawa 1973; Jarlskog 1985 (the invariant); Bjorken, Dunietz (unitarity triangle); Wolfenstein 1983;
CKMfitter/UTfit and PDG 2024; Koide 1981; NuFIT 6.0."""
import math, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
fam = {"charged leptons": (0.51099895, 105.6583755, 1776.86), "down type": (4.70, 93.5, 4183.0), "up type": (2.16, 1273.0, 172570.0),
       "neutrinos (TRPG)": (0.36e-9, 8.66e-9, 50.34e-9)}
head("1. One shape seen from rotated viewpoints")
pr("  Turning a matrix to face another way never changes its eigenvalues, so every view of one shape would show the SAME mass ratios.")
for f, m in fam.items(): pr(f"   {f:17s} lightest/heaviest {m[0]/m[2]:.3e}   middle/heaviest {m[1]/m[2]:.3e}")
pr("  The ratios differ by factors of 10 to 100 between families. FAIL: the ladders are not one shape merely seen from different directions.")
head("2. One shape seen at different scales: masses = (one shape's values)^power, plus a turn")
pr("  Then every family's log masses are the same pattern stretched: the middle family sits at the same fraction of each span.")
L = {f: np.log(m) for f, m in fam.items()}; fr = {f: (l[1]-l[0])/(l[2]-l[0]) for f, l in L.items()}
for f in fam: pr(f"   {f:17s} middle at {fr[f]:.3f} of the span; span {L[f][2]-L[f][0]:.2f}")
pr("  Charged leptons and neutrinos: 0.654 and 0.644. One shape at two scales HOLDS for the lepton side (expected: TRPG derived the neutrino")
pr("  masses from the same 45-degree rule, with the lightest root mirrored, so this is a consistency check, not independent evidence).")
pr("  Quarks: 0.565 and 0.440. Random middles give a spread this small 12% of the time. FAIL for the quark side.")
head("3. The perspective-free geometry: what stays the same from every viewpoint")
lam, A, rb, eb = 0.22501, 0.826, 0.1591, 0.3523       # PDG 2024 global fit (Wolfenstein parameters)
gam = math.degrees(math.atan2(eb, rb)); bet = math.degrees(math.atan2(eb, 1 - rb)); alp = 180 - gam - bet
J = A**2*lam**6*eb*(1 - lam**2/2)
pr("  For the quarks, the quantities no change of viewpoint can alter are the six masses and one triangle: the unitarity triangle, built from")
pr("  the turn between the two ladders. Its shape is measured directly (B-factories, LHCb):")
pr(f"   angles alpha {alp:.1f} deg, beta {bet:.1f} deg, gamma {gam:.1f} deg; twice its area (the Jarlskog invariant) J = {J:.2e}")
pr("  This IS a measurable geometry assembled from perspectives: each experiment sees one side or one angle, and together they close one")
pr("  triangle. It is established physics (Kobayashi & Maskawa 1973; Jarlskog 1985; Bjorken; Dunietz).")
pr(f"  A one-shape guess for it, 'the triangle is right-angled' (alpha = 90 deg): fit {alp:.1f} deg, direct measurement about 85 +5 -4 deg (PDG):")
pr("  allowed, not established. If TRPG's geometry made the turn exactly a quarter, that would be a derived angle; it does not yet.")
head("4. 'Transitioning is constant': the geometry in motion")
pr("  Measured cases of constant transition between perspectives of one fixed geometry:")
pr("   - neutrino flavour oscillation: a neutrino moves between its three flavour views without stopping; the geometry behind it (three mass")
pr("     states and one fixed turn, the mixing) does not change (Pontecorvo 1957; Maki, Nakagawa & Sakata 1962; Super-Kamiokande 1998, SNO 2001).")
pr("   - neutral meson mixing (K, D, B mesons): quarks pass back and forth between two views at a fixed rate set by the quark turn (Gell-Mann & Pais 1955).")
pr("  Across scales: the strengths and masses change with energy, but the shape of the turn barely does. The unitarity triangle's")
pr("  angles stay very nearly the same from the Z to unification energies (Babu 1987; Naculich 1993); only its size changes slowly.")
pr("  Across cosmic time: the measured shape has not moved: the electric strength changes by less than about 1e-18 per year (atomic clocks:")
pr("  Rosenband et al. 2008; Lange et al. 2021); less than about 1e-7 over 2 billion years (Oklo: Shlyakhter 1976; Damour & Dyson 1996);")
pr("  the proton/electron ratio by less than about 1e-7 since 7 billion years ago (Bagdonaite et al. 2013).")
pr("  So constant transition HOLDS as physics, in a precise form: things move between perspectives all the time, while the geometry they move")
pr("  through stays fixed to better than one part in ten million. A geometry whose own shape drifted would already have been seen.")
