"""Inputs-versus-predictions ledger for TRPG, version 1.3.2 (18 November 2026).

Every quantity the framework uses or produces is listed once, with one label:
  INPUT     taken from measurement
  FITTED    a choice (number or discrete option) made with today's data in view
  DERIVED   follows from the framework's rules plus inputs, and agrees with a measurement already made
  PREDICTED follows from the framework's rules plus inputs, and its deciding measurement has not returned
  EFFECTIVE an effective number of known physics (like the constants of the nuclear-mass formula), counted on neither side

The script recounts the inputs three ways, and re-runs the break tests for the results that move the count
(tau mass, down quark, lightest neutrino mass, Higgs mass) plus the strong-force strength, the JUNO lead and TM1 against JUNO 2026,
so that a reader can try to break them. It does not re-check the other DERIVED rows of the ledger (for example
the up quark or rho = 1): those are standard results or come from the scripts named in the paper and README. Nothing here is new physics; it is bookkeeping that anyone can rerun.
"""
import math

# ---------------------------------------------------------------------------
# 1. The like-for-like count. Standard Model + gravity + cosmology, counted the usual way.
SM = {
    "force strengths (3)": 3, "Higgs value and Higgs mass": 2, "charged-fermion masses (9)": 9,
    "quark mixing (3 angles + 1 phase)": 4, "strong-CP angle": 1,
    "neutrino masses (3)": 3, "neutrino mixing (3 angles, 1 Dirac phase, 2 Majorana phases)": 6,
    "gravity's strength": 1, "dark energy": 1, "dark matter amount": 1, "matter excess": 1,
    "seeds of galaxies (size, tilt)": 2,
}
n_sm = sum(SM.values())

# What TRPG removes from that list, and on what evidence.
SAVINGS = [
    # (quantity, how, status, deciding measurement)
    ("tau mass", "45 deg rule for leptons (Koide 1983; reason: family symmetry + maximum uncertainty)", "confirmed", "already measured"),
    ("down-quark mass", "Cabibbo relation m_s/m_d from the mixing angle (Gatto, Sartori, Tonin 1968)", "confirmed", "already measured"),
    ("Higgs mass", "self-coupling zero at the grain (Shaposhnikov and Wetterich 2010) with the far-side states", "pending", "sharper top mass and alpha_s"),
    ("lightest neutrino mass", "45 deg rule with one negative root (as Brannen 2006 proposed for neutrinos)", "pending", "JUNO ordering; cosmological mass sum"),
    ("neutrino mixing number 1 (theta_12 given theta_13)", "TM1 pattern (Xing and Zhou 2007; Albright and Rodejohann 2009)", "pending", "JUNO precision theta_12"),
    ("neutrino mixing number 2 (Dirac phase)", "TM1 pattern read as perspectives", "pending", "DUNE, Hyper-Kamiokande"),
]
n_conf = sum(1 for s in SAVINGS if s[2] == "confirmed")
n_all = len(SAVINGS)

# Choices made with today's data in view (not continuous knobs, but they are choices).
FITTED = [
    ("grain of space = Planck energy / sqrt2", "one of five equally natural readings (see derive_45_and_sqrt2); feeds alpha_s and the Higgs mass"),
    ("weight rule of the far-side weak pairs", "meeting energy squared / grain; one stated assumption (no direct link to the background value)"),
    ("coloured far-side (D) states where the strengths cross", "placed at 4.6e13 GeV with the three strengths taken as equal there (a boundary condition, not derived)"),
    ("TM1 rather than another trimaximal pattern", "exact tribimaximal and TM2 are disfavoured by data"),
    ("sign pattern of the neutrino roots (-,+,+)", "the all-positive pattern gives no solution with the measured splittings"),
    ("maximum uncertainty as the reason for the 45 deg rule", "a starting principle chosen knowing Koide holds for the charged leptons; it fails for quarks"),
]

EFFECTIVE = [
    ("containment pressure B = 60 (57-64) MeV/fm^3", "from the proton and Delta masses; an effective number of the strong force"),
    ("stiffness of the containment field", "from the lattice glueball mass 1.71-1.73 GeV"),
]

print("1. THE COUNT")
print(f"   Standard Model + gravity + cosmology, counted like for like: {n_sm}")
for k, v in SM.items(): print(f"      {v:2d}  {k}")
print(f"   Savings already confirmed by measurement: {n_conf}  -> TRPG takes {n_sm - n_conf} measured numbers today")
print(f"   Savings waiting for a measurement:        {n_all - n_conf}  -> {n_sm - n_all} if every one of them lands")
for q, how, st, dec in SAVINGS: print(f"      [{st:9s}] {q}: {how}; decided by {dec}")
print(f"   Choices made with today's data in view (FITTED, listed, not counted as numbers): {len(FITTED)}")
for q, why in FITTED: print(f"      {q}: {why}")
print("   Effective numbers of the strong force, counted on neither side:")
for q, why in EFFECTIVE: print(f"      {q}: {why}")
print("   Note: the grain choice feeds both alpha_s and the Higgs mass. If alpha_s lands at 0.1178 one force strength")
print("   is saved too, but the grain was chosen with data in view, so that saving only counts once both land.")
print("   Honest summary: 32 today, 28 if every pending measurement agrees, plus 6 stated choices.")
assert n_sm == 34 and n_sm - n_conf == 32 and n_sm - n_all == 28

# ---------------------------------------------------------------------------
# 2. Try to break each DERIVED label.
print("\n2. BREAK TESTS ON THE DERIVED LABELS")
me, mmu = 0.51099895, 105.6583755
# tau from Q = 2/3 with e and mu: solve for sqrt(m_tau)
a, b = math.sqrt(me), math.sqrt(mmu)
# Q(x) = (me+mmu+x^2)/(a+b+x)^2 = 2/3  ->  3(me+mmu+x^2) = 2(a+b+x)^2
A, B, C = 1.0, -4 * (a + b), 3 * (me + mmu) - 2 * (a + b) ** 2
x = (-B + math.sqrt(B * B - 4 * A * C)) / (2 * A)
mtau = x * x
for lab, val, err in (("PDG 2022", 1776.86, 0.12), ("PDG 2024", 1776.93, 0.09)):
    print(f"   tau mass from the 45 deg rule: {mtau:.3f} MeV vs {lab} {val} +- {err}: {(mtau - val) / err:+.1f} sigma")
print("   -> DERIVED holds (from a principle counted among the fitted choices). Not unique to TRPG: it is Koide's relation (1983). TRPG adds a reason for it.")
print("   Lepton masses: three masses from TWO numbers (scale and twist) plus the rule, not from one.")

# Down quark from the Cabibbo relation: sin(theta_C) = sqrt(m_d/m_s)  (Gatto, Sartori, Tonin 1968)
ratio = 1 / 0.2243 ** 2; md = 93.5 / ratio
print(f"   m_s/m_d from the Cabibbo angle: {ratio:.1f} vs lattice 20.4 +- 0.3 -> {(ratio - 20.4) / 0.3:+.1f} sigma; m_d = {md:.2f} MeV (lattice 4.70 +- 0.07)")
print("   (the tangent form, tan(theta_C) = sqrt(m_d/m_s), would give 18.9, about -5 sigma; the sine form is the one used)")
print("   -> DERIVED holds as a check. The same relation cannot also be counted as 'Cabibbo angle from the masses';")
print("      it is one relation, used once.")

# Neutrinos: Q = 2/3 with the lightest root negative, using the measured splittings
d21, d31 = 7.49e-5, 2.513e-3
def Q(m1):
    m2, m3 = math.sqrt(m1 * m1 + d21), math.sqrt(m1 * m1 + d31)
    s = -math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return (m1 + m2 + m3) / s ** 2 - 2 / 3
lo, hi = 1e-6, 2e-3
for _ in range(200):
    mid = (lo + hi) / 2
    if Q(lo) * Q(mid) <= 0: hi = mid
    else: lo = mid
m1 = (lo + hi) / 2; m2 = math.sqrt(m1 * m1 + d21); m3 = math.sqrt(m1 * m1 + d31)
print(f"   neutrinos, 45 deg rule with roots (-,+,+): m1 {m1*1e3:.2f}, m2 {m2*1e3:.2f}, m3 {m3*1e3:.2f} meV, sum {1e3*(m1+m2+m3):.1f} meV")
print("   -> PREDICTED (not measured). Brannen (2006) proposed the same relation for neutrinos; the numbers are shared.")

# Higgs mass and the top mass it needs
need, mt, emt = 171.63, 172.57, 0.29
print(f"   Higgs rule needs top pole mass {need} vs {mt} +- {emt}: {(mt - need) / emt:.1f} sigma on the top mass alone")
print("   -> PREDICTED, under strain. Not unique: any extra weak matter at high energy moves it the same way.")

# alpha_s
pa, lat, elat, wa, ewa = 0.1178, 0.11873, 0.00056, 0.1180, 0.0009
print(f"   alpha_s 0.1178 vs most precise lattice {lat} +- {elat}: {(lat - pa) / elat:+.1f} sigma; vs world average {wa} +- {ewa}: {(wa - pa) / ewa:+.1f} sigma")
print("   -> PREDICTED, rests on a FITTED grain choice; a 0.3% lattice result decides it.")

# Small neutrino gap (a lead, not core)
r = 3 - 2 * math.sqrt(2)
p21 = r * r * (0.36e-3 ** 2 + 2.519e-3) - 0.36e-3 ** 2
print(f"   lead: small gap from m2/m3 = 3-2sqrt2: {p21*1e5:.3f}e-5 eV^2 (7.37 +- 0.08 with JUNO's own large splitting)")
print("   vs global fits with JUNO's first data 7.48-7.53e-5 (+0.7 to +1.1 sigma); vs JUNO 2026 7.388 +- 0.078e-5 (within 0.2 sigma)")
print("   -> LEAD, holds against JUNO 2026; JUNO's full precision decides it.")

# theta_12 from TM1 against JUNO 2026 (arXiv:2609.37895)
s13 = 0.0220
t12 = 1 - 2 / (3 * (1 - s13))
print(f"   TM1 sin^2 theta12 = {t12:.3f} vs JUNO 2026 0.3036 +- 0.0064: {(t12 - 0.3036) / 0.0064:+.1f} sigma")
print("   -> PREDICTED, under pressure: 3 sigma would fail TM1 and with it the theta12 and CP-phase savings.")

print("\n3. LABELS THAT DID NOT SURVIVE AS 'DERIVED' (kept, relabelled in the ledger)")
print("   Cabibbo angle from the masses: same relation as the down-quark check; counted once.")
print("   Age of the universe, light bending, ringdown: standard results reproduced from measured inputs (H0, G); not TRPG predictions.")
print("   Heaviest neutron star 2.0-2.15 Msun: the observed window the pairing gap must respect, not a separate prediction.")
print("   The 45 deg rule for leptons: derived from a starting principle (maximum uncertainty) that was chosen knowing Koide holds.")
