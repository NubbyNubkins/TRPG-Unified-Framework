"""Hidden-role check (14 Nov 2026): the 18 quiet matrix cells (13 gravity-only + 5 loop-only) and the 18 'no use'
relationships, each marked REAL LINK (could hide roles; what would reveal it) or NUMBER MATCH (failed; nothing hidden)."""
import pandas as pd
F="TRPG_Quark_Gluon_Matrices.xlsx"
cells=[
("neutrinos x photon (3 cells)","REAL LINK, hidden role possible","a tiny neutrino magnetic moment or neutrino decay into light","Borexino, XENONnT limits (< ~6e-12 Bohr magnetons); not seen yet"),
("neutrinos x Higgs (3 cells)","REAL LINK, hidden role CERTAIN","neutrinos have mass, so something links them to mass: the Higgs (Dirac) or a heavy right-handed partner (seesaw)","KATRIN; neutrinoless double-beta decay (LEGEND, nEXO); DESI 2025 cosmology bound sum < ~0.064 eV vs the record's predicted 0.0592 eV: a live test"),
("neutrinos x gluon (3 cells)","REAL LINK only indirectly","neutrinos reach quarks through W and Z, never gluons directly","deep-inelastic neutrino scattering (IceCube, DUNE)"),
("electron/muon/tau x gluon (3 cells)","REAL LINK, hidden role possible","hypothetical leptoquarks would join leptons to colour directly","LHC leptoquark searches; B-meson decay tests (LHCb)"),
("gluon x W (1 cell)","loops only in practice","W and gluon meet only through quarks","LHC W + jets (known physics)"),
("photon x Higgs","LOOP, hidden role ALREADY FOUND","H -> 2 photons (via W and top loops) is how the Higgs was discovered in 2012","ATLAS, CMS"),
("gluon x Higgs","LOOP, hidden role ALREADY FOUND","gluon + gluon -> Higgs makes about 87% of all Higgs at the LHC","ATLAS, CMS"),
("photon x Z","LOOP, hidden role recently seen","Higgs -> Z + photon, first evidence 2023 (3.4 sigma)","ATLAS + CMS 2023"),
("gluon x Z","LOOP","gluon + gluon -> Z Z through quark loops (measured)","ATLAS, CMS"),
("photon x gluon","LOOP","light meets gluons only through quarks (photon-gluon fusion makes heavy quarks)","HERA, LHC"),
("THE WHOLE 'GRAVITY ONLY' CLASS","REAL, and possibly the biggest hidden role","dark matter is by definition something that, to us, does only gravity: it would sit in exactly this class","direct detection (LZ, XENONnT), LHC missing-energy searches, galaxy maps"),
]
rels=[
("Exact twist fractions","NUMBER MATCH, failed test","nothing hidden"),("Up/down opposite weak sides","NUMBER/IDEA, failed","nothing hidden"),
("Lean proportional to charge","NUMBER MATCH, failed","nothing hidden"),("Lean ratio 5/2 or 8/3","UNDECIDED","better light-quark masses (lattice QCD) will decide"),
("Home quark masses for 45 deg","NUMBER test, failed","nothing hidden"),("Charge-type link r^2","NUMBER MATCH, failed","nothing hidden"),
("4/9 three times","COINCIDENCE","nothing hidden unless a derivation links them"),("9/8 vs 8/9","COINCIDENCE","same"),
("2/9 vs Cabibbo vs weak angle","COINCIDENCE, ruled out","nothing hidden"),("Twist difference vs Vcb","NUMBER MATCH, failed","nothing hidden"),
("Band centres","toy numbers, failed","nothing hidden (band RANGES remain a hint)"),("Spiral angle for single masses","failed","nothing hidden"),
("Golden-ratio ladder","failed","nothing hidden"),("Hadron mass bands","none found","nothing hidden"),
("Literal 'right side has strong force'","IDEA, modified version survives","W_R / Z' searches at LHC; right-handed neutrinos"),
("3D space -> 3 colours","REAL IDEA, not shown","octonion / division-algebra derivations (Furey)"),
("3 colours -> 3 generations","REAL IDEA, no link known","same octonion route; CP-violation needs >= 3"),
("Visibility falls with distance","IDEA, replaced","replaced by 'stable + strong/EM = visible' (holds)"),
]
c=pd.DataFrame(cells,columns=["quiet cell(s)","verdict","possible hidden role","experiment or effect that could reveal it"])
r=pd.DataFrame(rels,columns=["'no use' relationship","verdict","what could change it"])
with pd.ExcelWriter(F,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    c.to_excel(w,sheet_name="Quiet_Cells_Hidden_Roles",index=False); r.to_excel(w,sheet_name="No_Use_Rechecked",index=False)
