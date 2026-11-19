"""Draw the TRPG derivation tree: rules -> outputs, each output coloured by its label (version 1.2.2)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

COL = {"rule": "#1F3A5F", "bet": "#6B4E9B", "input": "#7A7A7A", "fitted": "#C0392B",
       "effective": "#8E6C3A", "step": "#2E6F95", "derived": "#2E7D32", "predicted": "#B8860B"}
TXT = {"rule": "white", "bet": "white", "input": "white", "fitted": "white", "effective": "white",
       "step": "white", "derived": "white", "predicted": "black"}

N = {  # name: (x, y, text, kind)
    "R1": (1.5, 9.0, "Rule: one full medium,\none perspective at a time", "rule"),
    "R2": (7.1, 9.0, "Rule: containment\n(push balances pull)", "rule"),
    "R3": (9.9, 9.0, "Rule: maximum uncertainty\n(self vs relationship)", "rule"),
    "BET": (16.0, 9.0, "Bet: the far side of the\nrelationship matrix is real", "bet"),
    "R4": (22.4, 9.0, "Rule: nothing extra", "rule"),
    "OCT": (1.5, 6.6, "Octonions:\n3 colours, strong force", "derived"),
    "FOR": (1.5, 4.4, "Which forces exist:\n8 + 3 + 1 (theorem)", "derived"),
    "BAG": (7.1, 6.6, "B from proton + Delta\n(effective)", "effective"),
    "K45": (9.9, 6.6, "45° rule\nfor leptons", "step"),
    "TM1": (4.3, 6.6, "TM1 mixing read\nas perspectives", "step"),
    "MMG": (16.0, 6.6, "Far-side weak pairs\nat 2.4×10⁸ GeV", "step"),
    "IN": (19.6, 6.6, "28 measured inputs\n(+4 while pending)", "input"),
    "GEN": (1.5, 2.2, "Generation charges,\none-handed matter\n(stated readings)", "derived"),
    "GRV": (-1.3, 6.6, "Form of gravity\n(G measured)", "derived"),
    "FIT": (19.6, 4.4, "3 fitted choices: grain\n= Planck/√2, weight rule,\nequal at the crossing", "fitted"),
    "FIT2": (13.6, 4.4, "3 fitted choices: maximum\nuncertainty, TM1,\nneutrino root signs", "fitted"),
    "GAP": (7.1, 2.2, "Pairing gap\n82–97 MeV", "predicted"),
    "TAU": (9.9, 2.2, "Tau mass 1776.969 MeV\n(Koide)", "derived"),
    "NU": (12.7, 2.2, "Neutrinos: normal,\nsum 0.059 eV\n(Brannen)", "predicted"),
    "DCP": (4.3, 2.2, "sin²θ12 = 0.318,\nCP phase 253–292°", "predicted"),
    "AS": (15.5, 2.2, "α_s = 0.1178", "predicted"),
    "MH": (18.3, 2.2, "Higgs mass\n127.3 GeV", "predicted"),
    "PD": (22.4, 2.2, "No proton decay", "predicted"),
}
E = [("R1", "OCT", 0), ("OCT", "FOR", 0), ("FOR", "GEN", 0), ("R1", "GRV", 0), ("R2", "BAG", 0), ("BAG", "GAP", 0),
     ("R3", "K45", 0), ("K45", "TAU", 0), ("K45", "NU", 0), ("R1", "TM1", 0), ("TM1", "DCP", 0), ("R4", "PD", 0),
     ("BET", "MMG", 0), ("FIT", "MMG", 0), ("FIT", "MH", 0), ("MMG", "AS", 0), ("MMG", "MH", 0),
     ("IN", "MH", 0.35), ("IN", "AS", 0.1), ("IN", "NU", 0.1), ("IN", "TAU", 0.1), ("IN", "GAP", 0.1), ("IN", "DCP", 0.1), ("FIT2", "K45", 0), ("FIT2", "TM1", -0.3), ("FIT2", "NU", 0)]

fig, ax = plt.subplots(figsize=(22, 10.5), dpi=150)
ax.set_xlim(-3.0, 24.0); ax.set_ylim(0.3, 10.3); ax.axis("off")
W, H = 2.6, 1.1
P = {}
for k, (x, y, t, kind) in N.items():
    P[k] = FancyBboxPatch((x - W / 2, y - H / 2), W, H, boxstyle="round,pad=0.04,rounding_size=0.12",
                          fc=COL[kind], ec="white", lw=1.5, zorder=3)
    ax.add_patch(P[k])
    ax.text(x, y, t, ha="center", va="center", fontsize=9.3, color=TXT[kind], zorder=4)
for a, b, rad in E:
    dashed = a == "IN"
    ax.annotate("", xy=N[b][:2], xytext=N[a][:2], zorder=2,
                arrowprops=dict(arrowstyle="-|>", color="#9AA5B1" if dashed else "#34495E", lw=1.0 if dashed else 1.6,
                                linestyle="--" if dashed else "-", patchA=P[a], patchB=P[b], shrinkA=2, shrinkB=2,
                                connectionstyle=f"arc3,rad={rad}"))
leg = [("rule", "Rule"), ("bet", "Bet"), ("step", "Step in the chain"), ("derived", "Derived (agrees)"),
       ("predicted", "Predicted (awaiting)"), ("effective", "Effective number"), ("fitted", "Fitted choice"), ("input", "Measured inputs")]
for i, (k, t) in enumerate(leg):
    x = -2.6 + i * 3.25
    ax.add_patch(FancyBboxPatch((x, 0.55), 0.35, 0.3, boxstyle="round,pad=0.02", fc=COL[k], ec="none"))
    ax.text(x + 0.45, 0.7, t, va="center", fontsize=8.6)
ax.text(10.5, 10.05, "TRPG: from the rules to the outputs", ha="center", fontsize=15, weight="bold", color="#1F3A5F")
ax.text(10.5, 9.65, "Dashed lines: where measured inputs enter. Nothing is counted as established until its measurement returns.",
        ha="center", fontsize=9.5, color="#555")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "derivation_tree.png")
plt.savefig(out, bbox_inches="tight", facecolor="white"); print("saved", out)
