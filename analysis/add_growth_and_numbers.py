"""Adds growth (up-the-ladder), per-quark relational profiles, scale visibility and a decimals table
to TRPG_Quark_Gluon_Matrices.xlsx (14 Nov 2026). Values are taken only from the matrix sheets plus PDG 2024."""
import math, pandas as pd
from fractions import Fraction as Fr
F = "TRPG_Quark_Gluon_Matrices.xlsx"
mes = pd.read_excel(F, "Mesons_qqbar"); bar = pd.read_excel(F, "Baryons_qqq"); ckm = pd.read_excel(F, "CKM", index_col=0)
mass = dict(up=2.16, down=4.70, strange=93.5, charm=1273.0, bottom=4183.0, top=172570.0)
sym = dict(up="u", down="d", strange="s", charm="c", bottom="b", top="t")
v = 246219.65
# ---- per-quark profile from the matrix entries that contain it ----
prof = []
for q, s in sym.items():
    m_rows = mes[(mes.quark == q) | (mes.antiquark == "anti-" + q)].dropna(subset=["spin-0 mass (MeV)"])
    b_rows = bar[bar.content.str.split().apply(lambda L: s in L) & bar["status"].eq("observed")]
    hadron_masses = list(m_rows["spin-0 mass (MeV)"]) + list(b_rows["spin-1/2 mass (MeV)"].dropna()) + list(b_rows["spin-3/2 mass (MeV)"].dropna())
    shares = list(m_rows["field-energy share of spin-0 mass"].dropna()) + list(b_rows["field-energy share of mass"].dropna())
    lightest = min(hadron_masses) if hadron_masses else None
    row = {"quark": q, "mass (MeV)": mass[q], "observed hadrons containing it (in these sheets)": len(set(hadron_masses)),
           "lightest hadron (MeV)": lightest,
           "growth factor = lightest hadron / quark mass": round(lightest / mass[q], 3) if lightest else None,
           "mass added by binding to make lightest hadron (MeV)": round(lightest - mass[q], 1) if lightest else None,
           "average binding share over its hadrons": round(sum(shares) / len(shares), 3) if shares else None,
           "Higgs coupling y": round(math.sqrt(2) * mass[q] / v, 6)}
    if q in ckm.index: r = ckm.loc[q]; row["stays in own generation |V|^2"] = round(max(r) ** 2, 4); row["leaks to other generations"] = round(1 - max(r) ** 2 / sum(r ** 2), 4)
    else:
        c = ckm[q]; row["stays in own generation |V|^2"] = round(max(c) ** 2, 4); row["leaks to other generations"] = round(1 - max(c) ** 2 / sum(c ** 2), 4)
    prof.append(row)
prof = pd.DataFrame(prof)
# ---- growth: moving UP the ladder (production, binding, freeze-out) ----
growth = pd.DataFrame([
 ("photon -> electron + positron", "electron", "1.022 MeV", "pair production (Anderson 1932); two photons -> matter seen 2021 (STAR, Breit-Wheeler)"),
 ("quarks + gluons -> protons, neutrons, pions", "proton / muon-pion", "~150 MeV (1e-5 s after the big bang)", "quark-gluon plasma made at RHIC and LHC; transition 156 MeV (lattice QCD)"),
 ("proton + neutron -> deuterium -> helium", "proton (nuclei)", "~0.07 MeV (3 minutes)", "big bang nucleosynthesis: 25% helium by mass, matches measurement"),
 ("p + p -> deuterium + positron + neutrino (in the Sun)", "proton + weak", "~1.3 keV Sun core, tunnelling", "needs the weak force, which is why the Sun burns slowly (Bethe 1939)"),
 ("cosmic-ray proton + air -> pions -> muons", "muon-pion", "> 0.3 GeV", "about 1 muon per cm^2 per minute at sea level"),
 ("electron + nucleus -> atom", "electron (atoms)", "13.6 eV, 0.26 eV temperature (380,000 years)", "the cosmic microwave background was released then"),
 ("gluon + gluon -> charm, bottom, top pairs", "proton -> weak", "2.5 / 8.4 / 345 GeV", "LHC makes ~1 top pair per second"),
 ("gluon + gluon -> Higgs", "weak", "125 GeV", "main way the Higgs is made (LHC 2012)"),
 ("Higgs field switches on", "weak", "~160 GeV (1e-11 s)", "gives W, Z, quarks and leptons their mass (electroweak crossover)"),
 ("neutron capture -> heavy elements", "proton (nuclei)", "neutron star mergers, supernovae", "gold seen made in GW170817 kilonova (2017)"),
], columns=["growth step (going UP the ladder)", "scale group", "energy or temperature needed", "where it is seen"])
freeze = pd.DataFrame([
 ("weak scale", "~100 GeV", "1e-11 s", "Higgs on; W, Z, top get mass"),
 ("proton scale", "~1 GeV -> 150 MeV", "1e-5 s", "quarks bind into protons and neutrons"),
 ("muon-pion scale", "~100 MeV", "1e-4 s", "muons and pions disappear"),
 ("electron scale", "~1 MeV -> 0.5 MeV", "1 s -> 10 s", "neutrinos let go; electrons and positrons annihilate"),
 ("nuclear", "~0.07 MeV", "3 min", "helium forms"),
 ("atomic (electron binding)", "0.26 eV", "380,000 years", "atoms form, light set free"),
 ("neutrino / dark-energy scale", "~meV", "today (13.8 billion years)", "dark energy takes over the expansion"),
], columns=["scale group", "temperature (energy)", "time after big bang", "what grew or froze"])
# ---- visibility by scale ----
vis = pd.DataFrame([
 ("neutrino (meV)", "stable, weak only", "hardly at all", "dark energy (expansion speeding up); neutrino detectors; cosmic microwave glow is also at this energy"),
 ("electron (~1 MeV)", "stable, EM", "EVERYWHERE, as light and chemistry", "lightning, electricity, colour, light, magnetism, every chemical bond, the size of atoms"),
 ("muon-pion (~100 MeV)", "unstable (us to ns)", "indirectly", "the range of the nuclear force (pion -> 1.4 fm, so the size of nuclei); cosmic-ray muons in cloud chambers"),
 ("proton (~1 GeV)", "stable, strong + EM", "EVERYWHERE, as mass", "99% of the weight of everything; sunlight from fusion; nuclear power"),
 ("weak (~100 GeV)", "unstable (1e-25 s)", "only through effects", "radioactivity, why the Sun burns for billions of years, why anything has mass (Higgs)"),
 ("Planck / stars / black holes", "gravity", "as gravity", "orbits, stars, black holes (containment ladder)"),
], columns=["scale group", "stability and force", "how visible to us", "what its relationships produce that we can see"])
# ---- reasonable numbers as decimals ----
def fr(x): return "1/sqrt2" if abs(x-1/math.sqrt(2))<1e-12 else str(Fr(x).limit_denominator(50))
nums = [
 ("Koide balance (leptons)", 2/3, 0.666661, "leptons"), ("twist", 2/9, 0.222239, "leptons"), ("generation link", 1/math.sqrt(2), None, "leptons"),
 ("Koide value, down-type quarks", None, 0.731, "quarks"), ("Koide value, up-type quarks", None, 0.849, "quarks"),
 ("quark colour charge C_F", 4/3, None, "strong"), ("gluon colour charge C_A", 3, None, "strong"), ("gluon / quark C_A/C_F", 9/4, None, "strong"),
 ("quark / gluon C_F/C_A", 4/9, None, "strong"), ("quark-gluon pull vs quark-antiquark", 9/8, None, "strong"),
 ("quark-quark pull vs quark-antiquark", 1/2, None, "strong"), ("settled gluons' share", 1/4, None, "strong"), ("restless gluons' share", 3/4, None, "strong"),
 ("weak carrier / matter", 8/3, None, "weak"), ("weak matter charge C_F(SU2)", 3/4, None, "weak"),
 ("up-quark charge squared", 4/9, None, "EM"), ("down-quark charge squared", 1/9, None, "EM"),
 ("Buchdahl compactness 2GM/Rc^2", 8/9, None, "stars"), ("Buchdahl GM/Rc^2", 4/9, None, "stars"), ("universe compactness", 1, 1.0000, "cosmos"),
 ("Cabibbo |Vus|", None, 0.22431, "weak"), ("sqrt(m_d/m_s)", None, math.sqrt(4.70/93.5), "quarks"), ("sin^2 theta_W", None, 0.23121, "weak"),
 ("3M^2 / proton mass", None, 941.55/938.272, "leptons-proton"), ("spin constant (GeV^2)", None, 0.55, "mesons"),
 ("strong baseline / pion mass", None, 150/139.57, "strong"), ("alpha (EM)", None, 1/137.036, "EM"), ("alpha_s at Z", None, 0.1180, "strong"), ("alpha_W", None, 0.0338, "weak"),
]
num = pd.DataFrame([{"number": n, "exact fraction": fr(e) if e is not None else "", "exact decimal": round(e, 6) if e is not None else None,
                     "measured": round(m, 6) if m is not None else None, "where": w} for n, e, m, w in nums])
lineups = pd.DataFrame([
 ("0.4444 = 4/9", "quark/gluon charge ratio (strong); up-quark charge squared (EM); Buchdahl GM/Rc^2 (heaviest star)", "three unrelated derivations; COINCIDENCE as far as known physics goes. Worth noting, not a result."),
 ("1.125 = 9/8 and 0.8889 = 8/9", "quark-gluon pull vs quark-antiquark (strong) and Buchdahl switch point (gravity)", "reciprocals from different physics; coincidence unless a derivation links them"),
 ("0.2222 = 2/9 vs 0.2243 (Cabibbo) vs 0.2312 (weak angle)", "lepton twist; quark mixing; weak mixing", "close but NOT equal: Cabibbo differs from 2/9 by ~1%, ruled out as equal by the data"),
 ("0.2243 = Cabibbo vs 0.2242 = sqrt(m_d/m_s)", "quark mixing vs quark masses", "REAL relation, known since 1968 (Gatto, Sartori, Tonin); in the record as K4"),
 ("1.0035 = 3M^2 / proton", "lepton family scale vs proton", "unexplained, noticed by Brannen; hint for the proton scale view"),
 ("0.9912 = top quark Higgs coupling", "top quark vs the Higgs field", "almost exactly 1: the top is the one particle fully linked to the Higgs. Known, unexplained puzzle"),
 ("1.07 = strong baseline / pion", "proton layer baseline vs muon-pion scale", "real: the pion is light because of how the strong force's vacuum breaks (Nambu 1960)"),
], columns=["number", "where it appears", "honest reading"])
with pd.ExcelWriter(F, mode="a", engine="openpyxl", if_sheet_exists="replace") as w:
    for name, df in [("Quark_Relational_Profile", prof), ("Growth_Up_The_Ladder", growth), ("Big_Bang_Growth_Timeline", freeze),
                     ("Visibility_By_Scale", vis), ("Numbers_As_Decimals", num), ("Number_Lineups", lineups)]:
        df.to_excel(w, sheet_name=name, index=False)
print(prof.to_string(index=False)); print(num.to_string(index=False))
