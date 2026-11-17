"""Cell-by-cell value map of the TRPG matrices (14 Nov 2026). Every pair cell in All_Pairs (171) and
Particle_vs_Antiparticle (361) is sorted by what the relationship actually DOES (tree level, Standard Model):
  1 CHANGES  - one turns into another or they merge/annihilate (drives every decay and creation)
  2 BINDS    - holds them together (strong colour, nuclear, or electric attraction)
  3 PUSHES / SCATTERS - they affect each other but nothing is built or changed
  4 LOOPS ONLY - no direct link; only rare indirect (quantum-loop) effects
  5 GRAVITY ONLY - nothing but gravity (in the Standard Model)
"""
import pandas as pd
M="TRPG_Every_Particle_to_Every_Particle.xlsx"
Q="TRPG_Quark_Gluon_Matrices.xlsx"
P=pd.read_excel(M,"Particles").set_index("name"); A=pd.read_excel(M,"All_Pairs"); PA=pd.read_excel(M,"Particle_vs_Antiparticle")
ferm=set(P[P.spin==0.5].index); nu={"nu_e","nu_mu","nu_tau"}; nucl={"proton","neutron"}; quarks=set(P[P.family=="quark"].index)
def col(x): return bool(P.loc[x,"colored"]); 
def ch(x): return P.loc[x,"charge"]!=0
def massive(x): return P.loc[x,"mass_MeV"]>0
def classify(a,b,W):
    s={a,b}
    notes=[]
    if W>0: return 1,"W link: one turns into the other (decay / conversion channel), |amplitude| = %.3g"%W
    if "W" in s:
        o=(s-{"W"}).pop() if len(s)>1 else "W"
        if o in ferm: return 1,"W changes this particle into its weak partner"
        if o in {"photon","Z","Higgs"}: return 2 if o=="photon" else 1, {"photon":"W is charged: photon pulls/pushes it directly","Z":"W W Z vertex: carriers interact (non-abelian)","Higgs":"Higgs gives the W its mass; H -> W W"}[o]
    if "gluon" in s:
        o=(s-{"gluon"}).pop()
        if o in quarks: return 2,"gluon binds quarks (colour)"
        if o in nucl: return 2,"gluons hold the nucleon together (99% of its mass)"
        if o=="Higgs": return 4,"gluon + gluon -> Higgs only through a top-quark loop (but this is how most Higgs are made)"
        return 4 if o in {"photon","Z"} else 5, "gluon has no direct link to colourless %s"%o
    if "photon" in s:
        o=(s-{"photon"}).pop()
        if ch(o): return 3,"photon couples to its charge (light, electricity)"
        if o=="neutron": return 3,"neutron has a magnetic moment from its charged quarks"
        if o=="Higgs": return 4,"only through loops (H -> 2 photons is 0.23% of decays)"
        if o=="Z": return 4,"no direct link; only tiny loop effects"
        return 5,"neutral, no direct link to light"
    if "Z" in s:
        o=(s-{"Z"}).pop()
        if o=="Higgs": return 1,"Higgs gives the Z its mass; H -> Z Z"
        return 3,"Z couples weakly (neutral current)"
    if "Higgs" in s:
        o=(s-{"Higgs"}).pop()
        if o in nu: return 5,"no Higgs link in the Standard Model; but neutrinos DO have mass, so some link must exist (open question)"
        return 3,"Higgs couples in proportion to mass (y = %.3g)"%(1.4142*P.loc[o,"mass_MeV"]/246219.65)
    # fermion-fermion
    if col(a) and col(b): return 2,"strong colour link (they can bind in hadrons)"
    if s<=nucl: return 2,"nuclear force binds them (deuterium, every nucleus)"
    if (a in quarks and b in nucl) or (b in quarks and a in nucl): return 2,"the quark is part of / binds with the nucleon"
    q=P.loc[a,"charge"]*P.loc[b,"charge"]
    if q<0: return 2,"electric attraction (e.g. electron + proton = hydrogen)"
    if q>0: return 3,"electric repulsion"
    return 3,"only the weak neutral force (Z) connects them (very feeble)"
rows=[]
for _,r in A.iterrows():
    c,why=classify(r.particle_i,r.particle_j,r.W_transition_amplitude); rows.append((r.particle_i,r.particle_j,c,why))
cells=pd.DataFrame(rows,columns=["particle_i","particle_j","class","what it does"])
names={1:"1 CHANGES (drives decay/creation)",2:"2 BINDS",3:"3 PUSHES / SCATTERS only",4:"4 LOOPS ONLY",5:"5 GRAVITY ONLY"}
cells["class"]=cells["class"].map(names)
pmap={"can annihilate -> neutral carriers (photons, gluons, Z, Higgs)":names[1],"can fuse -> W boson (weak charged current)":names[1],
      "can annihilate -> pions (strong force; hadronic)":names[1],"flavour-changing: no direct annihilation (only rare loop processes)":names[4],
      "no annihilation; scatter only (colour cannot cancel)":names[3],"no annihilation; scatter only":names[3]}
pa=PA[["i","j","fate"]].copy(); pa["class"]=pa.fate.map(pmap)
summ=pd.concat([cells["class"].value_counts().rename("All_Pairs (171 particle-particle cells)"),
                pa["class"].value_counts().rename("Particle_vs_Antiparticle (361 cells)")],axis=1).fillna(0).astype(int).sort_index()
# who is most connected (valuable hubs)
deg=pd.concat([cells[cells["class"].str.startswith(("1","2"))][c] for c in ("particle_i","particle_j")]).value_counts().rename("links that change or bind (of 18)")
with pd.ExcelWriter(Q,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    summ.to_excel(w,sheet_name="Cell_Value_Summary"); cells.to_excel(w,sheet_name="Cell_Value_Map",index=False)
    pa.to_excel(w,sheet_name="Cell_Value_Map_Anti",index=False); deg.to_frame().to_excel(w,sheet_name="Most_Connected")
print(summ); print(deg.to_string())
print(cells[cells["class"].str.startswith(("4","5"))].to_string(index=False))
