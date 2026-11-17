"""Relationship ledger (14 Nov 2026): every relationship found in the TRPG matrices, sorted into
'valuable' (tight match or predictive) and 'does nothing (yet)' (coincidence, failed, or no use)."""
import pandas as pd
F="TRPG_Quark_Gluon_Matrices.xlsx"
V=[  # relationship, number, why valuable, source
("Lepton relational matrix (45 deg, 1/sqrt2 link, 2/9 twist)","e, mu, tau to 0.003%","3 masses from 1 input; predicted neutrino masses (sum 0.0592 eV)","Koide 1981, Brannen 2006; record M2"),
("Lepton family scale vs proton","3M^2/m_p = 1.0035","0.35% match; ties leptons to the proton scale view (unexplained)","Brannen; record K6"),
("Cabibbo from quark masses","|Vus| 0.2243 vs sqrt(m_d/m_s) 0.2242","0.05% match; mixing follows from masses","Gatto, Sartori, Tonin 1968; record K4"),
("Gell-Mann-Okubo baryon rule","0.57% off","predicts masses inside the table","Gell-Mann, Okubo 1961-62"),
("Equal spacing, spin-3/2 baryons","gaps 153, 149, 139 MeV","PREDICTED the Omega- before discovery (1962 -> 1964)","Gell-Mann; Barnes et al. 1964"),
("Identical-quark rule","uuu, ddd, sss, ccc, bbb only spin-3/2","forced the discovery of colour","Greenberg 1964; Han & Nambu 1965"),
("Spin-spin constant in mesons","M(1)^2 - M(0)^2 = 0.55 GeV^2 (+/-7%)","holds across light-heavy mesons; predicts unmeasured partners (e.g. Bc*)","De Rujula, Georgi, Glashow 1975"),
("Binding share falls with quark mass","proton 99% -> Upsilon 12%","containment (B1, J2) seen directly","PDG; this workbook"),
("Empty top row","top lives 5e-25 s < binding time","a gap that is a result","PDG"),
("Gluon/quark colour charge","9/4 = 2.25","measured in jets (LEP ~2.2)","Casimir scaling; OPAL, DELPHI, ALEPH"),
("Gluon groups: settled pair never couple","0 of 9 vertices","exact structure: 2 settled + 6 restless; B alone cannot exist","Gell-Mann 1962 algebra"),
("Missing 9th gluon","3 x 3-bar = 9 = 8 + 1","explains why the strong force is short-ranged","SU(3)"),
("3 colours fix quark charges","+2/3, -1/3","forced by consistency once N=3","Geng & Marshak 1989; record M16"),
("Lifetime tells the force","strong 1e-23 s, EM 1e-18 s, weak 1e-13 s and up","sorts every unstable particle to a matrix cell","PDG"),
("Decays fall down the scale ladder via W","t -> b -> c -> s -> u ; tau -> mu -> e","'links above and below' in real data","CKM, PDG"),
("Stability depends on surroundings","neutron, muon, proton all switch","relational stability; Fermi golden rule","Fermi 1934; neutron-star physics"),
("Own-clock rings Q = m/Gamma","rho 5 ... neutron 1e27","every particle completes its own cycle","PDG"),
("Big bang built the ladder in time","100 GeV, 150 MeV, 1 MeV, 0.26 eV","each scale grows at its own temperature (F2)","standard cosmology"),
("Strong baseline vs pion","150 MeV vs 139.6 (1.07)","proton-layer baseline sits at the muon/pion scale","Nambu 1960; record M8"),
("Dark-energy scale in neutrino mass range","2.24 meV in 0.36-50 meV","dark energy = neutrino scale's baseline (hint)","Fardon, Nelson, Weiner 2004; record M8"),
("W acts on right-handed ANTIparticles only","exact","answers the right-spin question","Lee, Yang, Wu 1956-57"),
("Higgs vs conjugate-Higgs side","down +2.5 deg, up +6.2 deg","families on the conjugate side lean most (lead)","Standard Model"),
("Generations: 3 <= N <= 8","CP needs >= 3; Z width 2.996; asymptotic freedom <= 8","brackets the number 3","Kobayashi & Maskawa 1973; LEP; Gross, Wilczek, Politzer"),
("Octonion: fix one direction -> 3 quaternion copies","always 3 (checked all 7)","geometric lead for 3 generations","Hurwitz; Furey 2018"),
("Top Higgs coupling","0.991","top fully linked to the Higgs (known puzzle)","PDG"),
("Every unstable particle has a matrix cell","24 of 24 (neutrino mixing added)","matrices are complete for decays","this workbook"),
]
N=[
("Exact twist fractions 2/9 x (1, 1/2, 1/3)","predicts m_u 2.83 vs 2.13, m_d 4.39 vs 4.68","FAILED prediction test; dropped"),
("Up/down opposite weak sides","both lean the same way","FAILED"),
("Lean proportional to charge","ratio 2.48 vs 2","FAILED (about 4 sigma)"),
("Lean ratio 5/2 or 8/3","both fit 2.48","cannot decide: number-hunting until masses improve"),
("Home (constituent) quark masses for the 45 deg rule","up 49 deg, down 32 deg","FAILED: families move opposite ways"),
("Charge-type link r^2 = 3/6, 4/6, 5/6","7 sigma off","FAILED (record Q5)"),
("0.4444 = 4/9 three times","colour ratio, up charge^2, Buchdahl","coincidence, no shared derivation"),
("9/8 vs 8/9","quark-gluon pull vs Buchdahl","coincidence"),
("2/9 vs Cabibbo vs weak angle","0.2222 / 0.2243 / 0.2312","close but not equal; ruled out"),
("Twist difference vs |Vcb|","0.037 vs 0.041","10% off; no match"),
("Band centres 0.05, 0.35, 1, 2, 4, 10","—","toy estimates, fail (record G1)"),
("Spiral angle for single masses","p = 0.75","no pattern (record F2)"),
("Golden-ratio ladder","3.76% vs random 3.00%","FAILED"),
("Mass bands for hadrons","smooth spread","no banding for composites"),
("Literal 'right side has the strong force'","right quarks feel gluons equally","FAILED as worded (modified version survives)"),
("3D space -> 3 colours","colour is an inner space","not shown (suggestive only)"),
("3 colours -> 3 generations","no known link","nothing"),
("Visibility falls with distance from proton","electron scale very visible","FAILED as worded; replaced by 'stable + strong/EM = visible'"),
]
v=pd.DataFrame(V,columns=["relationship","number","why it is valuable","credit"]); n=pd.DataFrame(N,columns=["relationship","number","why it does nothing (yet)"])
with pd.ExcelWriter(F,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    v.to_excel(w,sheet_name="Valuable_Relationships",index=False); n.to_excel(w,sheet_name="Relationships_No_Use",index=False)
def md(df): return "\n".join("| "+" | ".join(str(x) for x in r)+" |" for r in df.itertuples(index=False))
open("ledger.md","w").write("\n# an earlier test. Relationships that matter vs ones that do nothing (14 Nov)\n\n## Valuable: tight matches or predictive\n\n| relationship | number | why | credit |\n|---|---|---|---|\n"+md(v)+"\n\n## Seem to do nothing: coincidence, failed, or no use (yet)\n\n| relationship | number | why |\n|---|---|---|\n"+md(n)+"\n")
print(len(v),len(n))
