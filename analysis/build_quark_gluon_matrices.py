"""Quark and gluon relationship matrices (TRPG), built 14 Nov 2026.
Companion to TRPG_Every_Particle_to_Every_Particle.xlsx; same PDG quark masses.
Data: PDG Review of Particle Physics 2024 (masses, lifetimes, widths, PMNS); exotic hadrons from the LHCb, Belle,
BESIII, CMS and WASA-at-COSY papers named in the sheets; lattice QCD predictions from Brown, Detmold, Meinel & Orginos
(2014, baryons) and Morningstar & Peardon (1999) / Chen et al. (2006) (glueballs). Gluon algebra: Gell-Mann (1962).
Every number that is a prediction rather than a measurement is labelled as such.
"""
import itertools, math
import numpy as np, pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment

OUT = "TRPG_Quark_Gluon_Matrices.xlsx"
Q = ["up", "down", "strange", "charm", "bottom", "top"]
sym = dict(up="u", down="d", strange="s", charm="c", bottom="b", top="t")
S2Q = {v: k for k, v in sym.items()}
mass = dict(up=2.16, down=4.70, strange=93.5, charm=1273.0, bottom=4183.0, top=172570.0)  # MeV, PDG (same as main workbook)
charge = dict(up=2/3, down=-1/3, strange=-1/3, charm=2/3, bottom=-1/3, top=2/3)
gen = dict(up=1, down=1, strange=2, charm=2, bottom=3, top=3)
v_higgs = 246219.65

# ---------------- 1. The six quarks ----------------
quarks = pd.DataFrame([{
    "quark": q, "symbol": sym[q], "mass_MeV (PDG)": mass[q], "charge": round(charge[q], 4), "generation": gen[q],
    "colour states": 3,
    "left-handed: weak isospin T3": 0.5 if charge[q] > 0 else -0.5,
    "left-handed: feels the W?": "yes",
    "right-handed: weak isospin T3": 0,
    "right-handed: feels the W?": "no",
    "hypercharge Y (left / right)": f"1/6 / {'2/3' if charge[q] > 0 else '-1/3'}",
    "Higgs coupling y = sqrt2 m / v": math.sqrt(2) * mass[q] / v_higgs,
    "weak partner (same generation)": {"up": "down", "down": "up", "strange": "charm", "charm": "strange", "bottom": "top", "top": "bottom"}[q],
    "makes hadrons?": "no: decays in 5e-25 s, before a hadron can form (~3e-24 s)" if q == "top" else "yes",
} for q in Q])

# ---------------- 2. Quark + antiquark (mesons), 6 x 6 ----------------
# lightest spin-0 (pseudoscalar) and spin-1 (vector) state for each quark / antiquark pair; PDG 2024
# (name, mass MeV, lifetime s or None, width MeV or None, main decay, force that ends it)
M = {
 ("u","u"): [("pi0",134.9768,8.43e-17,None,"photon photon","EM"),("rho0",775.26,None,149.1,"pi+ pi-","strong")],
 ("d","d"): [("pi0 (same mixed state as u u-bar)",134.9768,8.43e-17,None,"photon photon","EM"),("omega",782.66,None,8.68,"pi+ pi- pi0","strong")],
 ("u","d"): [("pi+",139.57039,2.6033e-8,None,"mu+ nu_mu","weak"),("rho+",775.11,None,149.1,"pi+ pi0","strong")],
 ("d","u"): [("pi-",139.57039,2.6033e-8,None,"mu- anti-nu_mu","weak"),("rho-",775.11,None,149.1,"pi- pi0","strong")],
 ("u","s"): [("K+",493.677,1.2380e-8,None,"mu+ nu_mu","weak"),("K*+",891.67,None,51.4,"K pi","strong")],
 ("s","u"): [("K-",493.677,1.2380e-8,None,"mu- anti-nu_mu","weak"),("K*-",891.67,None,51.4,"K pi","strong")],
 ("d","s"): [("K0 (K-short 8.95e-11 s / K-long 5.12e-8 s)",497.611,8.954e-11,None,"pi pi","weak"),("K*0",895.55,None,47.3,"K pi","strong")],
 ("s","d"): [("anti-K0",497.611,8.954e-11,None,"pi pi","weak"),("anti-K*0",895.55,None,47.3,"K pi","strong")],
 ("s","s"): [("eta' (mixed with u,d)",957.78,None,0.188,"pi pi eta","strong"),("phi",1019.461,None,4.249,"K+ K-","strong")],
 ("c","u"): [("D0",1864.84,4.103e-13,None,"K- + anything","weak"),("D*0",2006.85,None,None,"D0 pi0","strong")],
 ("u","c"): [("anti-D0",1864.84,4.103e-13,None,"K+ + anything","weak"),("anti-D*0",2006.85,None,None,"anti-D0 pi0","strong")],
 ("c","d"): [("D+",1869.66,1.033e-12,None,"K- + anything","weak"),("D*+",2010.26,None,0.0834,"D0 pi+","strong")],
 ("d","c"): [("D-",1869.66,1.033e-12,None,"K+ + anything","weak"),("D*-",2010.26,None,0.0834,"anti-D0 pi-","strong")],
 ("c","s"): [("Ds+",1968.35,5.04e-13,None,"K K pi, tau nu","weak"),("Ds*+",2112.2,None,None,"Ds+ photon","EM")],
 ("s","c"): [("Ds-",1968.35,5.04e-13,None,"K K pi, tau nu","weak"),("Ds*-",2112.2,None,None,"Ds- photon","EM")],
 ("c","c"): [("eta_c",2983.9,None,30.5,"hadrons (via 2 gluons)","strong"),("J/psi",3096.900,None,0.0926,"hadrons (via 3 gluons), e+e-, mu+mu-","strong")],
 ("u","b"): [("B+",5279.34,1.638e-12,None,"anti-D0 + anything","weak"),("B*+",5324.71,None,None,"B+ photon","EM")],
 ("b","u"): [("B-",5279.34,1.638e-12,None,"D0 + anything","weak"),("B*-",5324.71,None,None,"B- photon","EM")],
 ("d","b"): [("B0",5279.66,1.519e-12,None,"D- + anything","weak"),("B*0",5324.71,None,None,"B0 photon","EM")],
 ("b","d"): [("anti-B0",5279.66,1.519e-12,None,"D+ + anything","weak"),("anti-B*0",5324.71,None,None,"anti-B0 photon","EM")],
 ("s","b"): [("Bs0",5366.92,1.520e-12,None,"Ds- + anything","weak"),("Bs*0",5415.8,None,None,"Bs0 photon","EM")],
 ("b","s"): [("anti-Bs0",5366.92,1.520e-12,None,"Ds+ + anything","weak"),("anti-Bs*0",5415.8,None,None,"anti-Bs0 photon","EM")],
 ("c","b"): [("Bc+",6274.47,5.10e-13,None,"J/psi + anything","weak"),("Bc*+ (not yet measured)",None,None,None,"Bc+ photon","EM")],
 ("b","c"): [("Bc-",6274.47,5.10e-13,None,"J/psi + anything","weak"),("Bc*- (not yet measured)",None,None,None,"Bc- photon","EM")],
 ("b","b"): [("eta_b",9398.7,None,10.0,"hadrons (via 2 gluons)","strong"),("Upsilon(1S)",9460.30,None,0.05402,"hadrons (via 3 gluons), lepton pairs","strong")],
}
NOTOP = "none: the top decays (to W + b) before it can bind"
mes_rows, grid = [], pd.DataFrame(index=[f"{q} (quark)" for q in Q], columns=[f"anti-{q}" for q in Q])
for a in Q:
    for b in Q:
        A, B = sym[a], sym[b]
        r = {"quark": a, "antiquark": "anti-" + b, "content": f"{A} {B}-bar", "net charge": round(charge[a] - charge[b], 4),
             "quark masses added (MeV)": mass[a] + mass[b], "M_ij = 2cosh S (TRPG pair number)": (mass[a] + mass[b]) / math.sqrt(mass[a] * mass[b]),
             "same flavour (can annihilate)?": "yes" if a == b else "no (flavour must change by the W)"}
        if "t" in (A, B):
            if a == b:
                r.update({"spin-0 state": "toponium: no lasting bound state; LHC (CMS 2025, ATLAS 2025) saw a threshold bump near 345 GeV"})
            else:
                r.update({"spin-0 state": NOTOP})
            r["stable?"] = "n/a (never forms)"; grid.loc[f"{a} (quark)", f"anti-{b}"] = "none (top decays first)"
        else:
            P0, V1 = M[(A, B)]
            r.update({"spin-0 state": P0[0], "spin-0 mass (MeV)": P0[1], "spin-0 lifetime (s)": P0[2], "spin-0 width (MeV)": P0[3],
                      "spin-0 main decay": P0[4], "spin-0 ended by": P0[5],
                      "spin-1 state": V1[0], "spin-1 mass (MeV)": V1[1], "spin-1 width (MeV)": V1[3], "spin-1 ended by": V1[5],
                      "field-energy share of spin-0 mass": (1 - (mass[a] + mass[b]) / P0[1]),
                      "stable?": "no (all mesons decay)"})
            if P0[1] and V1[1]:
                r["spin-1 minus spin-0 (MeV)"] = V1[1] - P0[1]
                r["M(spin-1)^2 - M(spin-0)^2 (GeV^2)"] = (V1[1]**2 - P0[1]**2) / 1e6
            grid.loc[f"{a} (quark)", f"anti-{b}"] = f"{P0[0].split(' (')[0]} {P0[1]:.0f}"
        mes_rows.append(r)
mesons = pd.DataFrame(mes_rows)

# ---------------- 3. Three quarks (baryons), all 56 flavour contents ----------------
# (spin-1/2 name, mass, lifetime/width text) , (spin-3/2 name, mass, text) ; None = not allowed or not seen
Bdat = {
 "uuu": (None, ("Delta++", 1232, "width 117 MeV, strong")),
 "uud": (("proton", 938.27208816, "stable (> 1e34 years)"), ("Delta+", 1232, "width 117 MeV, strong")),
 "udd": (("neutron", 939.56542052, "878.4 s, weak"), ("Delta0", 1232, "width 117 MeV, strong")),
 "ddd": (None, ("Delta-", 1232, "width 117 MeV, strong")),
 "uus": (("Sigma+", 1189.37, "8.02e-11 s, weak"), ("Sigma*+", 1382.83, "width 36 MeV, strong")),
 "uds": (("Lambda (and Sigma0 at 1192.6, 7.4e-20 s, EM)", 1115.683, "2.62e-10 s, weak"), ("Sigma*0", 1383.7, "width 36 MeV, strong")),
 "dds": (("Sigma-", 1197.449, "1.48e-10 s, weak"), ("Sigma*-", 1387.2, "width 39 MeV, strong")),
 "uss": (("Xi0", 1314.86, "2.90e-10 s, weak"), ("Xi*0", 1531.80, "width 9 MeV, strong")),
 "dss": (("Xi-", 1321.71, "1.64e-10 s, weak"), ("Xi*-", 1535.0, "width 10 MeV, strong")),
 "sss": (None, ("Omega-", 1672.45, "8.21e-11 s, weak (predicted by Gell-Mann 1962, found 1964)")),
 "udc": (("Lambda_c+", 2286.46, "2.02e-13 s, weak"), ("Sigma_c*+", 2517.5, "strong")),
 "uuc": (("Sigma_c++", 2453.97, "width 1.9 MeV, strong"), ("Sigma_c*++", 2518.41, "width 15 MeV, strong")),
 "ddc": (("Sigma_c0", 2453.75, "width 1.8 MeV, strong"), ("Sigma_c*0", 2518.48, "width 15 MeV, strong")),
 "usc": (("Xi_c+", 2467.71, "4.53e-13 s, weak"), ("Xi_c*+", 2645.1, "strong")),
 "dsc": (("Xi_c0", 2470.44, "1.51e-13 s, weak"), ("Xi_c*0", 2646.2, "strong")),
 "ssc": (("Omega_c0", 2695.2, "2.73e-13 s, weak"), ("Omega_c*0", 2765.9, "EM")),
 "ucc": (("Xi_cc++", 3621.6, "2.56e-13 s, weak (LHCb 2017)"), None),
 "dcc": (("Xi_cc+ (not confirmed; same mass as Xi_cc++ expected)", None, ""), None),
 "scc": (("Omega_cc+ (not seen; lattice ~3738)", None, ""), None),
 "ccc": (None, ("Omega_ccc++ (not seen; lattice ~4796)", None, "")),
 "udb": (("Lambda_b0", 5619.60, "1.47e-12 s, weak"), None),
 "uub": (("Sigma_b+", 5810.56, "width 5 MeV, strong"), ("Sigma_b*+", 5830.32, "width 9 MeV, strong")),
 "ddb": (("Sigma_b-", 5815.64, "width 5 MeV, strong"), ("Sigma_b*-", 5834.74, "width 10 MeV, strong")),
 "usb": (("Xi_b0", 5791.9, "1.48e-12 s, weak"), ("Xi_b*0", 5952.3, "strong")),
 "dsb": (("Xi_b-", 5797.0, "1.57e-12 s, weak"), ("Xi_b*-", 5955.3, "strong")),
 "ssb": (("Omega_b-", 6045.2, "1.64e-12 s, weak"), None),
 "ucb": (("Xi_cb+ (not seen; lattice ~6943)", None, ""), None), "dcb": (("Xi_cb0 (not seen; lattice ~6943)", None, ""), None),
 "scb": (("Omega_cb0 (not seen; lattice ~6998)", None, ""), None), "ccb": (("Omega_ccb+ (not seen; lattice ~8007)", None, ""), None),
 "ubb": (("Xi_bb0 (not seen; lattice ~10143)", None, ""), None), "dbb": (("Xi_bb- (not seen; lattice ~10143)", None, ""), None),
 "sbb": (("Omega_bb- (not seen; lattice ~10273)", None, ""), None), "cbb": (("Omega_cbb0 (not seen; lattice ~11195)", None, ""), None),
 "bbb": (None, ("Omega_bbb- (not seen; lattice ~14366)", None, "")),
}
order = "udscbt"
bar_rows = []
for combo in itertools.combinations_with_replacement(order, 3):
    key = "".join(combo)
    qs = [S2Q[c] for c in combo]
    ch = sum(charge[q] for q in qs)
    r = {"content": " ".join(combo), "net charge": round(ch, 4), "quark masses added (MeV)": sum(mass[q] for q in qs),
         "identical quarks": max(combo.count(c) for c in set(combo)),
         "colour": "always red+green+blue (antisymmetric) = colour neutral"}
    if "t" in key:
        r.update({"spin-1/2 state": NOTOP, "status": "never forms"})
    else:
        h, f = Bdat[key]
        if key in ("uuu", "ddd", "sss", "ccc", "bbb"):
            r["spin-1/2 state"] = "forbidden: three identical quarks only make spin-3/2 (the colour rule)"
        elif h:
            r["spin-1/2 state"], r["spin-1/2 mass (MeV)"], r["spin-1/2 lifetime / width"] = h
        if f:
            r["spin-3/2 state"], r["spin-3/2 mass (MeV)"], r["spin-3/2 lifetime / width"] = f
        seen = (h and h[1]) or (f and f[1])
        r["status"] = "observed" if seen else "predicted, not yet observed"
        light = (h[1] if h and h[1] else None) or (f[1] if f and f[1] else None)
        if light: r["field-energy share of mass"] = 1 - r["quark masses added (MeV)"] / light
        r["stable?"] = "YES (only stable hadron)" if key == "uud" else "no"
    bar_rows.append(r)
baryons = pd.DataFrame(bar_rows)

# ---------------- 4. Bigger combinations: counts and the observed exotic states ----------------
def n_multi(k, n=6): return math.comb(n + k - 1, k)
counts = pd.DataFrame([
 {"type": "single quark", "make-up": "q", "flavour contents possible (6 flavours)": 6, "contents possible without top": 5, "observed": "never alone (confined); top decays as a bare quark"},
 {"type": "meson", "make-up": "q + anti-q", "flavour contents possible (6 flavours)": 36, "contents possible without top": 25, "observed": "all 25 non-top contents seen"},
 {"type": "baryon", "make-up": "q q q", "flavour contents possible (6 flavours)": n_multi(3), "contents possible without top": n_multi(3, 5), "observed": f"{int((baryons.status=='observed').sum())} of {n_multi(3,5)} seen"},
 {"type": "tetraquark", "make-up": "q q + anti-q anti-q", "flavour contents possible (6 flavours)": n_multi(2)**2, "contents possible without top": n_multi(2, 5)**2, "observed": "about 10 well-established states (below)"},
 {"type": "pentaquark", "make-up": "q q q q + anti-q", "flavour contents possible (6 flavours)": n_multi(4) * 6, "contents possible without top": n_multi(4, 5) * 5, "observed": "about 6 states, all with a c c-bar pair (below)"},
 {"type": "hexaquark / dibaryon", "make-up": "6 q", "flavour contents possible (6 flavours)": n_multi(6), "contents possible without top": n_multi(6, 5), "observed": "deuteron (bound p+n, stable); d*(2380) resonance"},
])
exotic = pd.DataFrame([
 ("chi_c1(3872) a.k.a. X(3872)", "c c-bar u u-bar (tetraquark or D D* molecule)", 3871.64, "width 1.2 MeV", "Belle 2003"),
 ("Zc(3900)+", "c c-bar u d-bar", 3887.1, "width 28 MeV", "BESIII, Belle 2013"),
 ("Tcc+(3875)", "c c anti-u anti-d", 3874.83, "width 0.05 MeV (very long-lived for a strong state)", "LHCb 2021"),
 ("Zcs(3985)-", "c c-bar s u-bar", 3985.2, "width 13.8 MeV", "BESIII 2021"),
 ("Tcs0(2900)", "c s anti-u anti-d", 2866, "width 57 MeV", "LHCb 2020"),
 ("X(6900)", "c c anti-c anti-c (all-charm)", 6886, "width ~170 MeV", "LHCb 2020; ATLAS, CMS 2023"),
 ("Zb(10610), Zb(10650)", "b b-bar u d-bar", 10607.2, "width 18 MeV", "Belle 2011"),
 ("Pc(4312)+", "u u d c c-bar (pentaquark)", 4311.9, "width 10 MeV", "LHCb 2019"),
 ("Pc(4440)+", "u u d c c-bar", 4440.3, "width 21 MeV", "LHCb 2019"),
 ("Pc(4457)+", "u u d c c-bar", 4457.3, "width 6.4 MeV", "LHCb 2019"),
 ("Pcs(4459)0", "u d s c c-bar", 4458.8, "width 17 MeV", "LHCb 2021"),
 ("Pcs(4338)0", "u d s c c-bar", 4338.2, "width 7 MeV", "LHCb 2022"),
 ("deuteron", "u u d + u d d (p n bound)", 1875.61, "STABLE", "Urey 1931"),
 ("d*(2380)", "6 quarks, u u u d d d", 2380, "width 70 MeV", "WASA-at-COSY 2011/2014"),
 ("H-dibaryon", "u u d d s s", None, "predicted (Jaffe 1977), not found", "searches by Belle 2013 and others"),
], columns=["state", "quark content", "mass (MeV)", "lifetime / width", "first seen (experiment, year)"])

# ---------------- 5. Patterns the quark matrices show (known physics; tested here) ----------------
mp = {k: v[0][1] for k, v in M.items()}
pat = []
for lab, (a, b) in [("u d-bar  (rho - pi)", ("u","d")), ("u s-bar  (K* - K)", ("u","s")), ("c d-bar  (D* - D)", ("c","d")),
                    ("c s-bar  (Ds* - Ds)", ("c","s")), ("b u-bar  (B* - B)", ("b","u")), ("b s-bar  (Bs* - Bs)", ("b","s")),
                    ("c c-bar  (J/psi - eta_c)", ("c","c")), ("b b-bar  (Upsilon - eta_b)", ("b","b"))]:
    P0, V1 = M[(a, b)]
    pat.append({"pattern": "spin-spin constant: M(spin-1)^2 - M(spin-0)^2", "pair": lab, "value (GeV^2)": round((V1[1]**2 - P0[1]**2) / 1e6, 3),
                "reading": "about 0.55 whenever a light quark is present; grows for two heavy quarks"})
N, L, Sg, Xi = (938.272 + 939.565) / 2, 1115.683, (1189.37 + 1192.642 + 1197.449) / 3, (1314.86 + 1321.71) / 2
pat.append({"pattern": "Gell-Mann-Okubo (1961-62), spin-1/2 baryons: 2(N + Xi) = 3 Lambda + Sigma", "pair": "left side vs right side",
            "value (GeV^2)": None, "reading": f"{2*(N+Xi):.1f} vs {3*L+Sg:.1f} MeV, off by {abs(2*(N+Xi)-(3*L+Sg))/(3*L+Sg)*100:.2f}%"})
dec = [1232, 1384.6, 1533.4, 1672.45]
pat.append({"pattern": "Equal spacing, spin-3/2 baryons (Delta, Sigma*, Xi*, Omega)", "pair": "gaps in MeV",
            "value (GeV^2)": None, "reading": f"gaps {dec[1]-dec[0]:.0f}, {dec[2]-dec[1]:.0f}, {dec[3]-dec[2]:.0f}. Gell-Mann used the empty slot to predict the Omega- (1962); it was found in 1964 (Barnes et al.)."})
pat.append({"pattern": "Identical-quark rule", "pair": "uuu, ddd, sss, ccc, bbb",
            "value (GeV^2)": None, "reading": "three identical quarks never make spin-1/2. Explaining the Delta++ and Omega- is why colour was invented (Greenberg 1964; Han & Nambu 1965)."})
pat.append({"pattern": "Field (containment) energy share", "pair": "light hadrons vs heavy",
            "value (GeV^2)": None, "reading": "proton ~99% field energy; B meson ~21%; Upsilon ~12%. The lighter the quarks, the more the mass is the binding itself."})
pat.append({"pattern": "Top row is empty", "pair": "every combination with a top",
            "value (GeV^2)": None, "reading": "top lives ~5e-25 s, shorter than the time the strong force needs to bind (~3e-24 s). An empty row is a result, not a gap."})
patterns = pd.DataFrame(pat)

# ---------------- 6. Gluons: two groups ----------------
l = np.zeros((8, 3, 3), dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1; l[6][1,2]=-1j; l[6][2,1]=1j
l[7]=np.diag([1,1,-2])/math.sqrt(3)
T = l / 2
f = np.zeros((8,8,8)); d = np.zeros((8,8,8))
for a, b, c in itertools.product(range(8), repeat=3):
    f[a,b,c] = (-2j * np.trace((T[a]@T[b]-T[b]@T[a]) @ T[c])).real
    d[a,b,c] = (2 * np.trace((T[a]@T[b]+T[b]@T[a]) @ T[c])).real
grpA, grpB = [2, 7], [0, 1, 3, 4, 5, 6]     # lambda3, lambda8  |  the six colour-changing ones
gname = {i: f"lambda{i+1}" for i in range(8)}
# colour-flow labels (the six off-diagonal come as complex pairs of the real Gell-Mann ones)
flow = {0: "r g-bar + g r-bar", 1: "-i(r g-bar - g r-bar)", 2: "(r r-bar - g g-bar)/sqrt2", 3: "r b-bar + b r-bar",
        4: "-i(r b-bar - b r-bar)", 5: "g b-bar + b g-bar", 6: "-i(g b-bar - b g-bar)", 7: "(r r-bar + g g-bar - 2 b b-bar)/sqrt6"}
glu = pd.DataFrame([{"gluon (Gell-Mann label)": gname[i], "colour make-up": flow[i],
    "group": "A: settled (colour-neutral, diagonal)" if i in grpA else "B: restless (colour-changing, off-diagonal)",
    "changes a quark's colour?": "no" if i in grpA else "yes",
    "its own antiparticle?": "yes" if i in grpA else "pairs with its partner (r g-bar <-> g r-bar)",
    "share of a quark's colour charge it carries": None} for i in range(8)])
CF_part = {i: (T[i] @ T[i]).real[0,0] + (T[i] @ T[i]).real[1,1] + (T[i] @ T[i]).real[2,2] for i in range(8)}
tot = sum(CF_part.values())
glu["share of a quark's colour charge it carries"] = [CF_part[i] / tot for i in range(8)]
# complex colour-flow basis for pair/net-colour work
wq = {"r": np.array([0.5, 1/(2*math.sqrt(3))]), "g": np.array([-0.5, 1/(2*math.sqrt(3))]), "b": np.array([0, -1/math.sqrt(3)])}
cf = [("r g-bar","r","g"),("g r-bar","g","r"),("r b-bar","r","b"),("b r-bar","b","r"),("g b-bar","g","b"),("b g-bar","b","g"),
      ("diag-3 (r r-bar - g g-bar)",None,None),("diag-8 (r r-bar + g g-bar - 2 b b-bar)",None,None)]
def wt(x): return np.zeros(2) if x[1] is None else wq[x[1]] - wq[x[2]]
gA = [c for c in cf if c[1] is None]; gB = [c for c in cf if c[1] is not None]
pair_rows = []; pgrid = pd.DataFrame(index=[c[0] for c in cf], columns=[c[0] for c in cf])
for x, y in itertools.combinations_with_replacement(cf, 2):
    w = wt(x) + wt(y); neutral = bool(np.allclose(w, 0))
    gx = "A" if x[1] is None else "B"; gy = "A" if y[1] is None else "B"
    pair_rows.append({"gluon 1": x[0], "gluon 2": y[0], "groups": "".join(sorted(gx + gy)),
                      "net colour (T3, T8)": f"({w[0]:+.3f}, {w[1]:+.3f})",
                      "can be part of a colour-neutral pair (glueball)?": "yes" if neutral else "no (net colour left over)"})
    pgrid.loc[x[0], y[0]] = pgrid.loc[y[0], x[0]] = "neutral" if neutral else "coloured"
gpairs = pd.DataFrame(pair_rows)
# 3-gluon couplings
trip = []
for a, b, c in itertools.combinations(range(8), 3):
    if abs(f[a,b,c]) > 1e-9 or abs(d[a,b,c]) > 1e-9:
        nA = sum(i in grpA for i in (a,b,c))
        trip.append({"gluons": f"{gname[a]}, {gname[b]}, {gname[c]}", "how many from group A": nA,
                     "f_abc (3-gluon vertex: the gluons pull on each other)": round(f[a,b,c], 4) if abs(f[a,b,c]) > 1e-9 else 0,
                     "d_abc (makes a 3-gluon colour-neutral state)": round(d[a,b,c], 4) if abs(d[a,b,c]) > 1e-9 else 0})
gtrip = pd.DataFrame(trip)
def closes(idx):
    """does the set close under commutators? returns list of generators produced outside the set"""
    out = set()
    for a, b in itertools.combinations(idx, 2):
        for c in range(8):
            if abs(f[a,b,c]) > 1e-9 and c not in idx: out.add(gname[c])
    return sorted(out)
def nvert(idx): return sum(1 for a,b,c in itertools.combinations(idx,3) if abs(f[a,b,c])>1e-9)
def npairs_neutral(g, h=None):
    if h is None: P = list(itertools.combinations_with_replacement(g, 2))
    else: P = [(x, y) for x in g for y in h]
    return sum(np.allclose(wt(x) + wt(y), 0) for x, y in P), len(P)
nAA, tAA = npairs_neutral(gA); nBB, tBB = npairs_neutral(gB); nAB, tAB = npairs_neutral(gA, gB); nall, tall = npairs_neutral(cf)
gres = pd.DataFrame([
 {"set": "Group A alone (2 settled gluons)", "members": 2,
  "share of quark colour charge": round(sum(CF_part[i] for i in grpA)/tot, 4),
  "self-interaction (3-gluon vertices inside the set)": nvert(grpA),
  "closed by itself?": "yes" if not closes(grpA) else "no",
  "colour-neutral pairs / all pairs": f"{nAA} / {tAA}",
  "what it would be on its own": "two copies of an electromagnetism-like force: no self-pull, so no confinement; long range"},
 {"set": "Group B alone (6 restless gluons)", "members": 6,
  "share of quark colour charge": round(sum(CF_part[i] for i in grpB)/tot, 4),
  "self-interaction (3-gluon vertices inside the set)": nvert(grpB),
  "closed by itself?": "no: their pulls create " + ", ".join(closes(grpB)),
  "colour-neutral pairs / all pairs": f"{nBB} / {tBB}",
  "what it would be on its own": "cannot exist alone: two restless gluons meeting produce a settled one"},
 {"set": "A with B (cross pairs only)", "members": "2 x 6",
  "share of quark colour charge": None, "self-interaction (3-gluon vertices inside the set)": "every vertex has >= 2 from B",
  "closed by itself?": "", "colour-neutral pairs / all pairs": f"{nAB} / {tAB}",
  "what it would be on its own": "a settled gluon never neutralises a restless one"},
 {"set": "Both groups combined (all 8)", "members": 8, "share of quark colour charge": 1.0,
  "self-interaction (3-gluon vertices inside the set)": nvert(range(8)),
  "closed by itself?": "yes: this is SU(3), the strong force",
  "colour-neutral pairs / all pairs": f"{nall} / {tall}",
  "what it would be on its own": "confinement, glueballs, and quarks bound into hadrons"},
 {"set": "The missing 9th (r r-bar + g g-bar + b b-bar)/sqrt3", "members": 1, "share of quark colour charge": 0,
  "self-interaction (3-gluon vertices inside the set)": 0, "closed by itself?": "yes (it would be its own U(1))",
  "colour-neutral pairs / all pairs": "it is colour-neutral alone",
  "what it would be on its own": "3 colours x 3 anticolours = 9, but nature uses only 8. A 9th would be colour-neutral, escape the proton and give a long-range strong force between all hadrons, which is not seen."},
])
glueballs = pd.DataFrame([
 ("2 gluons", "0++ (scalar)", "~1710 (lattice 1730 +/- 80)", "f0(1710), f0(1500) candidates", "Morningstar & Peardon 1999; Chen et al. 2006"),
 ("2 gluons", "2++ (tensor)", "~2400", "candidates debated", "same"),
 ("2 gluons", "0-+ (pseudoscalar)", "~2590", "X(2370) candidate (BESIII 2024)", "same; BESIII 2024"),
 ("3 gluons", "1-- (vector, needs d_abc)", "~3830", "none", "Chen et al. 2006"),
 ("2 gluons -> Higgs", "fusion", "125200", "Higgs made mostly by gluon fusion; Higgs -> gg is 8.2% of decays", "PDG 2024"),
], columns=["made of", "spin-parity", "mass (MeV)", "experimental status", "source"])

# ---------------- 7. Every unstable particle -> where it sits in the matrices ----------------
U = [
 ("muon","2.197e-6 s","weak","mu -> nu_mu + W-* -> e anti-nu_e","All_Pairs: muon-nu_mu (W amplitude 1)"),
 ("tau","2.903e-13 s","weak","tau -> nu_tau + W-*","All_Pairs: tau-nu_tau"),
 ("neutron","878.4 s","weak","d -> u W- (inside the neutron)","All_Pairs: up-down, |Vud| = 0.974"),
 ("top","~5e-25 s","weak","t -> b W+","All_Pairs: top-bottom, |Vtb| = 0.999"),
 ("W","3.2e-25 s (width 2.085 GeV)","weak","W -> lepton+neutrino or quark pair","every W-amplitude cell, run backwards"),
 ("Z","2.6e-25 s (width 2.495 GeV)","weak (neutral)","Z -> f f-bar","Particle_vs_Antiparticle: every direct-opposite cell"),
 ("Higgs","1.6e-22 s (width ~4.1 MeV)","Higgs coupling","H -> b b-bar 58%","All_Pairs: Higgs_coupling_product; Quarks sheet: y = sqrt2 m/v"),
 ("pi+ / pi-","2.60e-8 s","weak","u d-bar -> W+ -> mu+ nu","Mesons u,anti-d + Particle_vs_Antiparticle up/anti-down 'fuse -> W'"),
 ("pi0","8.4e-17 s","EM","u u-bar, d d-bar -> 2 photons","Mesons u,anti-u (direct opposites annihilate)"),
 ("K+ / K-","1.24e-8 s","weak","u s-bar -> W+ (|Vus| = 0.224)","Mesons u,anti-s + All_Pairs up-strange"),
 ("K0 short / long","8.95e-11 s / 5.12e-8 s","weak","s -> u W- ; K0 <-> anti-K0 mixing","Mesons d,anti-s and s,anti-d (the two cells mix: oscillation)"),
 ("eta, eta'","width 1.3 keV / 0.19 MeV","EM / strong","2 photons, 3 pions","Mesons diagonal cells (u u-bar, d d-bar, s s-bar mixed)"),
 ("rho, omega, phi, K*","widths 4-149 MeV (~1e-23 s)","strong","gluon makes a new q q-bar pair -> 2-3 light mesons","Mesons spin-1 column + Gluon sheets"),
 ("D, Ds","4.1e-13 to 1.0e-12 s","weak","c -> s W+ (|Vcs| = 0.975)","Mesons c rows + All_Pairs charm-strange"),
 ("B, Bs, Bc","5.1e-13 to 1.6e-12 s","weak","b -> c W- (|Vcb| = 0.041); B0 <-> anti-B0 mixing","Mesons b rows + All_Pairs charm-bottom"),
 ("J/psi, Upsilon","widths 93 keV / 54 keV","strong (via 3 gluons) + EM","c c-bar -> 3 gluons or lepton pair","Mesons diagonal cells + Gluon_Triples (d_abc)"),
 ("Lambda, Sigma+-, Xi, Omega","8e-11 to 3e-10 s","weak","s -> u W-","Baryons sheet + All_Pairs up-strange"),
 ("Sigma0","7.4e-20 s","EM","Sigma0 -> Lambda + photon","Baryons uds row (two states with the same quarks)"),
 ("Delta, Sigma*, Xi*","widths 9-117 MeV","strong","-> nucleon or hyperon + pion","Baryons spin-3/2 column + Gluon sheets"),
 ("Lambda_c, Xi_c, Omega_c, Xi_cc","1.5e-13 to 4.5e-13 s","weak","c -> s W+","Baryons c rows"),
 ("Lambda_b, Xi_b, Omega_b","1.5e-12 s","weak","b -> c W-","Baryons b rows"),
 ("tetraquarks, pentaquarks","widths 0.05-170 MeV","strong","fall apart into ordinary mesons and baryons","Exotic_States sheet"),
 ("glueballs (predicted)","~1e-23 s","strong","-> pairs of mesons","Glueballs sheet"),
 ("neutrino flavour (oscillation, not decay)","~ L/E (km per MeV)","weak (mixing)","nu_e <-> nu_mu <-> nu_tau","NEW: PMNS sheet (was missing from the matrices)"),
]
unstable = pd.DataFrame(U, columns=["unstable particle", "lifetime or width", "force that ends it", "what happens", "where it is in the matrices"])
# lifetime ladder by force
ladder = pd.DataFrame([("strong", "1e-24 to 1e-20 s", "rho, Delta, J/psi, exotic states"), ("electromagnetic", "1e-20 to 1e-16 s", "pi0, Sigma0, eta"),
                       ("weak", "1e-25 s (heavy W, Z, top) up to 880 s (neutron)", "muon, pions, kaons, D, B, hyperons, neutron")],
                      columns=["force", "typical lifetime", "examples"])

# PMNS (PDG 2024, normal ordering)
s12, s23, s13 = math.sqrt(0.307), math.sqrt(0.558), math.sqrt(0.0219); dl = 1.19 * math.pi
c12, c23, c13 = [math.sqrt(1 - x*x) for x in (s12, s23, s13)]; e = np.exp(1j*dl)
Upm = np.array([[c12*c13, s12*c13, s13/e],
                [-s12*c23 - c12*s23*s13*e, c12*c23 - s12*s23*s13*e, s23*c13],
                [s12*s23 - c12*c23*s13*e, -c12*s23 - s12*c23*s13*e, c23*c13]])
pmns = pd.DataFrame(np.abs(Upm).round(3), index=["nu_e", "nu_mu", "nu_tau"], columns=["nu_1", "nu_2", "nu_3"])
ckm = pd.DataFrame([[0.97367, 0.22431, 0.00382], [0.221, 0.975, 0.0408], [0.0086, 0.0415, 0.999]],
                   index=["up", "charm", "top"], columns=["down", "strange", "bottom"])

handed = pd.DataFrame([
 ("W boson (weak charged force)", "left-handed particles ONLY", "right-handed antiparticles ONLY", "Lee & Yang 1956; Wu 1957; Garwin, Lederman & Weinrich 1957"),
 ("Z boson (weak neutral force)", "both, but unequally", "both, but unequally", "Glashow 1961, Weinberg 1967, Salam 1968"),
 ("photon (EM), gluons (strong), gravity", "both equally", "both equally", ""),
 ("Higgs", "links left to right (that is how mass works)", "same", "Higgs, Englert & Brout 1964"),
 ("right-handed neutrino", "never seen; would feel no force but gravity and the Higgs", "", "'sterile' neutrino searches, open"),
 ("W_R (a right-handed W)", "hypothetical: would act on right-handed particles only", "", "Pati & Salam 1974; Mohapatra & Senjanovic 1975; LHC limit above ~6 TeV"),
], columns=["force carrier", "acts on particles that are", "acts on antiparticles that are", "source"])

# ---------------- 8. Forces by side: particle vs antiparticle, and one side acting on the other ----------------
forces_side = pd.DataFrame([
 ("strong (gluons)", "colour (r, g, b)", "anticolour (r-bar, g-bar, b-bar)", "acts on both, both handednesses",
  "quark + antiquark in the neutral (singlet) state: pull -4/3, twice the quark-quark pull of -2/3 (see Colour_Pair_Factors)", "gluons are their own anti-set: the force side has no separate mirror"),
 ("electromagnetic (photon)", "charge +q", "charge -q", "acts on both, both handednesses",
  "opposite charges pull, so a particle and its own antiparticle always pull (e- e+ forms positronium, 0.125 ns / 142 ns)", "photon is its own antiparticle"),
 ("weak charged (W+, W-)", "weak isospin +-1/2 (left-handed only)", "weak isospin -+1/2 (right-handed only)", "particles: LEFT only; antiparticles: RIGHT only",
  "W turns a left particle into its partner, and a right antiparticle into its partner; a particle and an antiparticle partner fuse into a W (pi+ = u + anti-d -> W+)", "W+ and W- are each other's antiparticle: the only carrier pair with a mirror"),
 ("weak neutral (Z)", "g = T3 - Q sin^2(theta_W)", "same size, opposite sign", "both, unequally (left stronger)",
  "every particle-antiparticle direct opposite can annihilate to a Z", "Z is its own antiparticle"),
 ("Higgs", "Yukawa y = sqrt2 m / v", "same y", "joins a left to a right state", "f + anti-f -> H strongest for the heaviest (top, bottom)", "Higgs is its own antiparticle"),
 ("gravity", "mass-energy (always +)", "mass-energy (always +): antimatter falls DOWN", "both, equally",
  "always pull, same as matter on matter (ALPHA-g at CERN, 2023: antihydrogen falls at 0.75 +/- 0.29 g)", "graviton (if it exists) is its own antiparticle"),
], columns=["force", "what it grabs on a particle", "what it grabs on the antiparticle", "handedness", "particle side acting on antiparticle side", "the force carrier's own mirror"])

# colour pair factors: <T1.T2> = (C_R - C_1 - C_2)/2, negative = pull, positive = push (Casimir scaling)
Cas = {"1": 0, "3": 4/3, "3bar": 4/3, "6": 10/3, "6bar": 10/3, "8": 3, "10": 6, "10bar": 6, "15": 16/3, "15bar": 16/3, "27": 8}
pairs_c = [("quark + antiquark", "3", "3bar", ["1", "8"]), ("quark + quark", "3", "3", ["3bar", "6"]),
           ("antiquark + antiquark", "3bar", "3bar", ["3", "6bar"]), ("quark + gluon", "3", "8", ["3", "6bar", "15"]),
           ("antiquark + gluon", "3bar", "8", ["3bar", "6", "15bar"]), ("gluon + gluon", "8", "8", ["1", "8", "8", "10", "10bar", "27"])]
cp = []
for lab, a, b, outs in pairs_c:
    for R in dict.fromkeys(outs):
        k = (Cas[R] - Cas[a] - Cas[b]) / 2
        cp.append({"pair": lab, "combined colour state": R, "colour-neutral?": "YES" if R == "1" else "no",
                   "strength factor (T1.T2)": round(k, 4), "pull or push": "pull" if k < 0 else ("push" if k > 0 else "none"),
                   "relative to quark-antiquark neutral pull": round(k / (-4/3), 4)})
colour_pairs = pd.DataFrame(cp)

# second condition: matter side vs force side ("opposite pairing" like antiparticles)
def CF(N): return (N*N - 1) / (2*N)
mvc = pd.DataFrame([
 ("strong SU(3)", "quark: 3 colour states, spin 1/2", "gluon: 8 colour states, spin 1", 4/3, 3.0, 3.0/(4/3), "9/4", "yes (8 = 3 x 3-bar minus 1)",
  "gluon jets spray 9/4 times more than quark jets: measured (LEP: OPAL, DELPHI, ALEPH 1990s, ratio ~2.2)"),
 ("weak SU(2)", "left lepton/quark doublet: 2 states, spin 1/2", "W+, W-, W3: 3 states, spin 1", CF(2), 2.0, 2.0/CF(2), "8/3", "yes (3 = 2 x 2 minus 1)", "hidden: the Higgs gives W, Z mass, so it is not directly seen"),
 ("electromagnetic U(1)", "charged fermion: charge 1, spin 1/2", "photon: charge 0, spin 1", 1.0, 0.0, 0.0, "0", "carrier carries no charge, so no self-pull", "why light passes through light"),
 ("gravity", "mass, spin 1/2 or any", "graviton (spin 2)", None, None, None, "all energy carries the 'charge'", "carrier carries the charge it responds to (gravity pulls on gravity)", "nonlinear GR; graviton not observed"),
], columns=["force", "matter side", "force-carrier side", "matter charge-squared C_F", "carrier charge-squared C_A", "carrier / matter", "as a fraction", "carriers made from matter x anti-matter?", "evidence"])
gen_N = pd.DataFrame([{"N (number of states)": N, "matter C_F = (N^2-1)/2N": CF(N), "carrier C_A = N": N, "ratio C_A/C_F = 2N^2/(N^2-1)": 2*N*N/(N*N-1),
                       "carriers = N^2 - 1": N*N - 1, "N x N-bar (incl. the missing neutral one)": N*N} for N in range(2, 7)])

readme = pd.DataFrame({"TRPG quark and gluon relationship matrices (14 Nov 2026)": [
 "What this is: every quark alone, every quark with an antiquark (36), every three-quark mix (56), bigger mixes, and the 8 gluons in two groups.",
 "Data: PDG 2024 unless a sheet says otherwise. Predictions (lattice QCD) are labelled 'not seen' or 'lattice'.",
 "Gluon groups: no gluon is stable on its own (all are confined). So 'easily stable' is read as: Group A = the 2 colour-neutral gluons that never change a quark's colour and are their own antiparticle; Group B = the 6 colour-changing gluons.",
 "Honest caveat on the groups: which two gluons count as 'diagonal' depends on how you label the colours (a choice of axes). The split into 2 + 6 is always the same; the members can be rotated.",
 "Sheets: Quarks, Meson_Matrix, Mesons_qqbar, Baryons_qqq, Combination_Counts, Exotic_States, Quark_Patterns, Gluons, Gluon_Pair_Matrix, Gluon_Pairs, Gluon_Triples, Gluon_Group_Results, Glueballs, Unstable_Coverage, Lifetime_Ladder, CKM, PMNS, Handedness, Forces_by_Side, Colour_Pair_Factors, Matter_vs_Carrier, Matter_vs_Carrier_SU_N.",
 "Credit: Gell-Mann and Zweig (quarks, 1964), Gell-Mann and Okubo (mass formulas), Greenberg, Han and Nambu (colour), Fritzsch, Gell-Mann and Leutwyler (QCD 1973), Gross, Wilczek and Politzer (1973), De Rujula, Georgi and Glashow (1975, spin-spin splitting), Cabibbo (1963), Kobayashi and Maskawa (1973), Pontecorvo, Maki, Nakagawa and Sakata (neutrino mixing), Lee, Yang and Wu (handedness), and the experiments named in each sheet.",
]})

sheets = [("Read_Me", readme, False), ("Quarks", quarks, False), ("Meson_Matrix", grid, True), ("Mesons_qqbar", mesons, False),
          ("Baryons_qqq", baryons, False), ("Combination_Counts", counts, False), ("Exotic_States", exotic, False),
          ("Quark_Patterns", patterns, False), ("Gluons", glu, False), ("Gluon_Pair_Matrix", pgrid, True), ("Gluon_Pairs", gpairs, False),
          ("Gluon_Triples", gtrip, False), ("Gluon_Group_Results", gres, False), ("Glueballs", glueballs, False),
          ("Unstable_Coverage", unstable, False), ("Lifetime_Ladder", ladder, False), ("CKM", ckm, True), ("PMNS", pmns, True), ("Handedness", handed, False),
          ("Forces_by_Side", forces_side, False), ("Colour_Pair_Factors", colour_pairs, False), ("Matter_vs_Carrier", mvc, False), ("Matter_vs_Carrier_SU_N", gen_N, False)]
with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    for name, df, idx in sheets: df.to_excel(w, sheet_name=name, index=idx)
wb = load_workbook(OUT)
for ws in wb:
    for c in ws[1]: c.font = Font(bold=True, color="FFFFFF"); c.fill = PatternFill("solid", fgColor="1F4E79"); c.alignment = Alignment(wrap_text=True, vertical="top")
    for col in ws.columns:
        L = max(len(str(c.value)) if c.value is not None else 0 for c in col)
        ws.column_dimensions[col[0].column_letter].width = min(max(10, L + 2), 70 if ws.title != "Read_Me" else 160)
    ws.freeze_panes = "B2" if ws.title != "Read_Me" else None
wb.save(OUT)

# ---------------- printed summary ----------------
print("baryon contents (no top):", n_multi(3,5), " observed:", int((baryons.status=='observed').sum()))
print(patterns[["pair","value (GeV^2)","reading"]].to_string(index=False))
print(gres[["set","share of quark colour charge","self-interaction (3-gluon vertices inside the set)","closed by itself?","colour-neutral pairs / all pairs"]].to_string(index=False))
print("3-gluon vertices with 2+ from group A:", int((gtrip[gtrip['f_abc (3-gluon vertex: the gluons pull on each other)']!=0]['how many from group A']>=2).sum()))
print(pmns)
print(mesons[["content","spin-0 state","field-energy share of spin-0 mass"]].dropna().to_string(index=False))

print(colour_pairs.to_string(index=False))
print(mvc[["force","carrier / matter","as a fraction"]].to_string(index=False))
