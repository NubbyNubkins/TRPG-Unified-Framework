"""Actual 'sizes' of the needles, how they interact, and a relational picture (no false spatial look).
Masses: PDG 2024 (quarks: u, d, s at 2 GeV; c, b at their own mass; t pole). Neutrinos: TRPG prediction.
Turns: CKM (Wolfenstein, PDG 2024 fit), PMNS (NuFIT 5.2, normal ordering: 33.41, 49.1, 8.54 deg, delta 197 deg)."""
import numpy as np, math, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
def pr(*a): print(*a, flush=True)
hbarc = 197.3269804  # MeV fm
K = [("Up quarks", ["up", "charm", "top"], [2.16, 1273.0, 172570.0], "#d1495b"),
     ("Down quarks", ["down", "strange", "bottom"], [4.70, 93.5, 4183.0], "#2e86ab"),
     ("Charged leptons", ["electron", "muon", "tau"], [0.51100, 105.658, 1776.93], "#3b8b3b"),
     ("Neutrinos (TRPG)", ["nu1", "nu2", "nu3"], [0.36e-9, 8.66e-9, 50.34e-9], "#8e6cc2")]
pr("1. 'Size' of each needle. The axes are masses (energies), not lengths. Two honest ways to turn them into something visible:")
pr("   (a) a drawing at true proportions with 1 mm = the electron's mass; (b) each mass's reach in space, hbar/(m c) (Compton 1923).")
for name, ps, ms, _ in K:
    pr(f"  {name}:")
    for p, m in zip(ps, ms):
        mm = m/0.511; draw = f"{mm:.3g} mm" if mm < 1e3 else f"{mm/1e3:.3g} m"
        if mm < 1e-3: draw = f"{mm*1e6:.3g} nm"
        reach = hbarc/m; rs = f"{reach:.3g} fm" if reach < 1e5 else (f"{reach*1e-9:.3g} micrometres" if reach >= 1e9 else f"{reach*1e-6:.3g} nm")
        pr(f"    {p:9s} mass {m:12.6g} MeV   drawn at 1 mm per electron: {draw:>10s}   reach hbar/mc: {rs}")
# turns
lam, A, rb, eb = 0.22501, 0.826, 0.1591, 0.3523
s12 = lam; s23 = A*lam**2; z = A*lam**3*(rb + 1j*eb)*math.sqrt(1 - A**2*lam**4)/(math.sqrt(1 - lam**2)*(1 - A**2*lam**4*(rb + 1j*eb)))
def mix(s12, s23, s13, d):
    c12, c23, c13 = (math.sqrt(1 - x*x) for x in (s12, s23, s13)); e = np.exp(1j*d)
    return np.array([[c12*c13, s12*c13, s13/e], [-s12*c23 - c12*s23*s13*e, c12*c23 - s12*s23*s13*e, s23*c13],
                     [s12*s23 - c12*c23*s13*e, -c12*s23 - s12*c23*s13*e, c23*c13]])
CKM = mix(s12, s23, abs(z), -np.angle(z))
PM = mix(*(math.sin(math.radians(t)) for t in (33.41, 49.1, 8.54)), math.radians(197))
P1, P2 = np.abs(CKM)**2, np.abs(PM)**2
pr("\n2. How the needles touch: the weak-force turn as chances (squared entries). Rows up/charm/top -> columns down/strange/bottom:")
for r, n in zip(P1, ["up", "charm", "top"]): pr(f"    {n:6s} " + "  ".join(f"{100*x:8.4f}%" for x in r))
pr("   Neutrino flavours e/mu/tau (rows) in mass states nu1/nu2/nu3 (columns):")
for r, n in zip(P2, ["e", "mu", "tau"]): pr(f"    {n:6s} " + "  ".join(f"{100*x:6.1f}%" for x in r))
pr("   Gluons and photons touch every family equally (they do not see the needles' shapes). The Higgs link sets each axis length")
pr("   (mass = link x 174 GeV). The far-side press reaches only the down needle (TRPG, an earlier test).")
# picture
fig, ax = plt.subplots(figsize=(16, 10.5)); ax.set_xlim(-1.2, 12.6); ax.set_ylim(-1.6, 8.6); ax.axis("off")
fig.suptitle("The same geometry drawn as relationships: who touches whom, and how strongly", fontsize=17, weight="bold", y=0.975)
X = [2.0, 4.6, 7.2]; Y = {0: 7.2, 1: 4.9, 2: 2.5, 3: 0.2}
def rad(m_mev): return 0.16 + 0.032*(math.log10(m_mev*1e6) + 4)
for f, x in enumerate(X): ax.text(x, 8.25, ["first family", "second family", "third family"][f], ha="center", fontsize=13, weight="bold", color="#444")
for k, (name, ps, ms, col) in enumerate(K):
    ax.text(-1.1, Y[k], name, va="center", fontsize=12.5, weight="bold", color=col)
    for f, (p, m) in enumerate(zip(ps, ms)):
        ax.add_patch(Circle((X[f], Y[k]), rad(m), color=col, alpha=0.85, zorder=3))
        lab = (f"{m/1e3:.3g} GeV" if m >= 1e3 else f"{m:.3g} MeV") if m >= 1 else (f"{m*1e3:.3g} keV" if m >= 1e-3 else f"{m*1e6:.3g} eV")
        if k in (0, 2): ax.text(X[f] + rad(m) + 0.12, Y[k] + 0.1, f"{p}\n{lab}", ha="left", va="center", fontsize=9.5, zorder=4)
        else: ax.text(X[f] + rad(m) + 0.12, Y[k] - 0.1, f"{p}\n{lab}", ha="left", va="center", fontsize=9.5, zorder=4)
def links(P, ya, yb, col, scale):
    for i in range(3):
        for j in range(3):
            w = P[i, j]
            if w < 1e-6: continue
            ax.plot([X[i], X[j]], [ya, yb], color=col, lw=max(0.6, scale*math.sqrt(w)), alpha=0.55 if w > 1e-3 else 0.9,
                    ls="-" if w > 1e-3 else ":", zorder=1)
links(P1, Y[0] - 0.62, Y[1] + 0.6, "#555", 9)
links(P2, Y[2] - 0.55, Y[3] + 0.35, "#555", 9)
ax.text(8.45, (Y[0] + Y[1])/2, "weak-force turn (W):\nup 94.9% to down,\n5.1% to strange,\n0.0014% to bottom", fontsize=9.5, va="center")
ax.text(8.45, (Y[2] + Y[3])/2, "neutrino turn:\nlarge mixes\n(e = 68% nu1, 30% nu2, 2% nu3)", fontsize=9.5, va="center")
ax.add_patch(FancyBboxPatch((9.9, Y[1] - 0.45), 2.4, 0.9, boxstyle="round,pad=0.08", fc="#fff1d6", ec="#e09f3e", lw=2))
ax.text(11.1, Y[1], "far-side heavy partner\n(10^7 to 10^8 GeV, TRPG)", ha="center", va="center", fontsize=9.5)
ax.annotate("", xy=(X[2] + 1.5, Y[1]), xytext=(9.85, Y[1]), arrowprops=dict(arrowstyle="->", color="#e09f3e", lw=2.5))
ax.text(9.0, Y[1] - 0.75, "presses only this row\n(squeeze ~12%, with a timing)", fontsize=9, color="#a5650f", ha="center")
key = ("How to read: circle size grows with the logarithm of the mass (the axis length of the needle). "
       "Line thickness = the chance of turning from one particle into the other through the weak force (dotted = below 0.1%).\n"
       "Not drawn as lines: the Higgs (it sets every circle's size: mass = link x 174 GeV); gluons and photons (they touch every family equally, "
       "so they do not see the shapes);\nthe hidden twist (a phase; its size is the six equal triangles, J/2 = 1.6e-5). "
       "Neutrino masses are TRPG's prediction; the rest are measured (PDG 2024, NuFIT 5.2).")
ax.text(-1.1, -1.05, key, fontsize=10, va="top", wrap=True)
out = "TRPG_Needle_Relationship_Map.png"
plt.savefig(out, dpi=170, bbox_inches="tight"); pr("\nsaved", out)
