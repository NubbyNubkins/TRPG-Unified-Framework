"""Prior-art check (14 Nov 2026) for every relationship in Valuable_Relationships and Relationships_No_Use."""
import pandas as pd
F="TRPG_Quark_Gluon_Matrices.xlsx"
v=pd.read_excel(F,"Valuable_Relationships"); n=pd.read_excel(F,"Relationships_No_Use")
KV=["KNOWN|Koide 1981; Foot 1994; Brannen 2006","KNOWN, unexplained|Brannen 2006","KNOWN|Gatto, Sartori & Tonin 1968; Cabibbo & Maiani 1968; Fritzsch 1977",
"KNOWN|Gell-Mann 1961; Okubo 1962","KNOWN|Gell-Mann 1962; Barnes et al. 1964","KNOWN|Greenberg 1964; Han & Nambu 1965",
"KNOWN|De Rujula, Georgi & Glashow 1975","KNOWN (standard QCD)|Ji 1995 proton mass decomposition; lattice QCD","KNOWN|Bigi et al. 1986 (top decays before hadronising)",
"KNOWN|Casimir scaling; LEP jet measurements","KNOWN (SU(3) algebra)|Gell-Mann 1962; Cartan","KNOWN|Fritzsch, Gell-Mann & Leutwyler 1973","KNOWN|Geng & Marshak 1989; Minahan, Ramond & Warner 1990",
"KNOWN|textbook (Griffiths, PDG)","KNOWN facts; framed as TRPG's 'links above and below'|CKM, PDG","KNOWN|Fermi 1934; Bethe; neutron-star physics (muons, hyperons)",
"KNOWN quantity (m/Gamma); 'own clock' reading is TRPG's framing|PDG","KNOWN|Kolb & Turner 1990; Weinberg 1977","KNOWN|Nambu 1960; Gell-Mann, Oakes & Renner 1968",
"KNOWN, unexplained|Fardon, Nelson & Weiner 2004","KNOWN|Lee & Yang 1956; Wu 1957","SIMILAR: the Higgs-side fact is known; linking it to the Koide lean appears NEW|this project",
"KNOWN|Kobayashi & Maskawa 1973; LEP Z width; Gross, Wilczek & Politzer 1973","KNOWN maths; used for generations by Furey 2018, Dixon|Hurwitz 1898; Furey 2018",
"KNOWN puzzle|PDG","OUR COMPILATION (facts known)|this workbook"]
KN=["APPARENTLY NEW (and failed)|quark Koide phases studied by Zenczykowski 2012 (arXiv:1210.4125)","APPARENTLY NEW (failed)|—",
"SIMILAR ideas (modified Koide for quarks)|Kartavtsev 2011 (arXiv:1111.0480); arXiv:1512.06349","APPARENTLY NEW (undecided)|—",
"SIMILAR (Koide with constituent masses tried before)|Rivero & Gsponer 2005; Wikipedia 'Koide formula' survey","SIMILAR|Rivero 2011; Rodejohann & Zhang 2011",
"trivial coincidence, not published as a claim|—","coincidence, not published as a claim|—","noticed informally by others; ruled out|—","APPARENTLY NEW (failed)|—",
"Original (toy numbers)|earlier versions","an earlier test|—","golden-ratio mass ladders proposed by several amateur and fringe papers; failed here|—",
"KNOWN: hadrons show no banding|PDG","SIMILAR (left-right and mirror models)|Pati & Salam 1974; Abbott & Farhi 1981",
"SIMILAR (division-algebra programme)|Dixon 1994; Furey 2018","SIMILAR (geometric generations)|Candelas et al. 1985; Furey 2018","Original framing; replaced|this project"]
def split(L): return [x.split("|")[0] for x in L],[x.split("|")[1] for x in L]
v["prior art"],v["citation"]=split(KV); n["prior art"],n["citation"]=split(KN)
extra=pd.DataFrame([("Every-particle-to-every-particle relationship matrices (this project's format)","SIMILAR presentations exist (vertex/interaction charts, the 'Elementary Particle Cube' arXiv:1703.06023); the all-pairs matrix with M_ij = 2cosh S, antiparticle and quark/gluon sheets, and the cell-value sorting appears NEW as a single relational framework","this project")],
                   columns=["relationship","prior art","citation"])
with pd.ExcelWriter(F,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    v.to_excel(w,sheet_name="Valuable_Relationships",index=False); n.to_excel(w,sheet_name="Relationships_No_Use",index=False)
    extra.to_excel(w,sheet_name="Prior_Art_Framework",index=False)
print(v["prior art"].str.split().str[0].value_counts()); print(n["prior art"].str.split().str[0].value_counts())
