"""Antiquark-sector sheets (14 Nov 2026): antibaryons (3 antiquarks) mirrored from Baryons_qqq, plus the up/down
'opposite sides' test. CPT: an antiparticle has exactly the same mass and lifetime, opposite charges
(antiproton/proton charge-to-mass ratio equal to 16 parts per trillion, BASE 2022)."""
import pandas as pd
F="TRPG_Quark_Gluon_Matrices.xlsx"
b=pd.read_excel(F,"Baryons_qqq")
a=b.copy()
a["content"]=["anti-"+" anti-".join(c.split()) for c in b.content]
a["net charge"]=-b["net charge"]
a["colour"]="always anti-red + anti-green + anti-blue = colour neutral"
for col in [c for c in a.columns if "state" in c and "mass" not in c and "lifetime" not in c]:
    a[col]=["anti-"+str(x) if isinstance(x,str) and not x.startswith(("forbidden","none","never")) else x for x in b[col]]
a["stable?"]=["YES while kept away from matter (no antiproton decay ever seen; annihilates on contact with matter)" if c=="u u d" else s for c,s in zip(b.content,b.get("stable?",[""]*len(b)))]
test=pd.DataFrame([
 ("leptons (e, mu, tau)","Higgs H","45.0","0.0","colourless"),
 ("down-type (d, s, b)","Higgs H (same as leptons)","47.5","+2.5 (68%: 2.3-2.6)","coloured, charge -1/3"),
 ("up-type (u, c, t)","CONJUGATE Higgs (the Higgs' 'anti' side)","51.2","+6.2","coloured, charge +2/3"),
],columns=["family","gets its mass from","balance angle (deg)","offset from 45 (deg)","notes"])
test2=pd.DataFrame([
 ("opposite offsets tied to weak isospin T3 (+1/2 / -1/2)","FAILS: both families lean the SAME way (both above 45)"),
 ("offset proportional to electric charge (ratio 2)","FAILS narrowly: measured ratio 2.48 (68%: 2.36-2.74); 2 is about 4 sigma off with standard masses (2.31 at a common scale)"),
 ("offset tied to which side of the Higgs gives the mass","PATTERN HOLDS qualitatively: leptons and down-type share the Higgs and sit at or near 45; up-type uses the conjugate Higgs and leans furthest. Not yet a formula."),
],columns=["idea","result"])
with pd.ExcelWriter(F,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    a.to_excel(w,sheet_name="Antibaryons_antiqqq",index=False); test.to_excel(w,sheet_name="UpDown_Opposite_Sides",index=False)
    test2.to_excel(w,sheet_name="UpDown_Opposite_Tests",index=False)
print(a[["content","net charge"]].head(4))
