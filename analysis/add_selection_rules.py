"""Why the quiet cells are quiet (14 Nov 2026): for each cell in Cell_Value_Map, the conserved charge that blocks a
direct link. A zero in the matrix is a correct prediction, not a wrong table entry."""
import pandas as pd
Q="TRPG_Quark_Gluon_Matrices.xlsx"
c=pd.read_excel(Q,"Cell_Value_Map")
def rule(r):
    s={r.particle_i,r.particle_j}; w=r["what it does"]
    if "gluon" in s and "colourless" in w: return "colour: the gluon only touches colour charge; the partner has none"
    if "photon" in s and ("neutral" in w or "loops" in w or "no direct" in w): return "electric charge: the photon only touches charge; the partner has none"
    if "Higgs" in s and "no Higgs link" in w: return "mass link: the Standard Model gives neutrinos none (yet they have mass: open)"
    if "only the weak neutral" in w: return "no shared colour or charge; only the Z (weak) is left"
    if "repulsion" in w: return "same-sign charges: they push, they can't bind"
    return "active link (charge, weak or Higgs coupling) that scatters but does not build or change anything"
c["blocking rule (why it is quiet)"]=c.apply(rule,axis=1)
q=c[c["class"].str.startswith(("3","4","5"))]
with pd.ExcelWriter(Q,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    q.to_excel(w,sheet_name="Why_Quiet_Cells_Are_Quiet",index=False)
print(q["blocking rule (why it is quiet)"].value_counts().to_string())
