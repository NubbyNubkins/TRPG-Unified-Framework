import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as P
labels=["u r","u g","u b","d r","d g","d b","ν","e",   "u r","u g","u b","d r","d g","d b","ν","e"]
n=16; fig,ax=plt.subplots(figsize=(9,11))
col={"gluon":"#E8C547","weak":"#6BA4E0","hyper":"#B8B8B8","higgs":"#D9534F","held":"#2E7D32"}
for i in range(n):
    for j in range(n):
        Li,Lj=i<8,j<8; qi,qj=(i%8)<6,(j%8)<6
        c="white"
        if i==j: c=col["held"]
        elif Li==Lj and qi and qj and (i%8)//3==(j%8)//3: c=col["gluon"]
        elif Li and Lj and ((i%8)%3==(j%8)%3 and qi and qj or (not qi and not qj)): c=col["weak"]
        elif Li!=Lj and (i%8)==(j%8): c=col["higgs"]
        ax.add_patch(P.Rectangle((j,n-1-i),1,1,fc=c,ec="#ddd",lw=.5))
ax.set_xlim(0,n); ax.set_ylim(0,n); ax.set_aspect("equal")
ax.set_xticks([k+.5 for k in range(n)]); ax.set_xticklabels(labels,fontsize=8)
ax.set_yticks([k+.5 for k in range(n)]); ax.set_yticklabels(labels[::-1],fontsize=8)
ax.axhline(8,color="k",lw=1.5); ax.axvline(8,color="k",lw=1.5)
ax.text(4,16.15,"left-handed side",ha="center"); ax.text(12,16.15,"right-handed side",ha="center")
h=[P.Patch(color=col["held"],label="held state (matter, on the diagonal)"),
   P.Patch(color=col["gluon"],label="gluon link (colour to colour)"),
   P.Patch(color=col["weak"],label="weak link (up↔down, ν↔e; left side only)"),
   P.Patch(color=col["higgs"],label="Higgs link (left ↔ right of the SAME particle) = resistance / weight"),]
ax.legend(handles=h,loc="upper center",bbox_to_anchor=(.5,-.06),ncol=1,fontsize=9,frameon=False)
ax.set_title("One generation's relationship matrix: matter on the diagonal, forces on the links\n(the electric/hypercharge link sits on every diagonal cell as a phase; gravity is the grid itself)",fontsize=10,pad=22)
plt.savefig("figures/relationship_matrix_diagonal_links.png",dpi=150,bbox_inches="tight")
