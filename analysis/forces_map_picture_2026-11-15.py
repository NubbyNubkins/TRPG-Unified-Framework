"""Picture 2: where each force pushes and pulls, which carrier (particle or wave) causes it, who feels it, relative strengths."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
fig = plt.figure(figsize=(18, 11))
fig.suptitle("Where each force pushes and pulls, what carries it, and who feels it (one hydrogen atom, a neutron, and a passing neutrino)", fontsize=15.5, weight="bold", y=0.98)
ax = fig.add_axes([0.01, 0.06, 0.56, 0.86]); ax.set_xlim(0, 12); ax.set_ylim(0, 10); ax.set_aspect("equal"); ax.axis("off")
COL = dict(strong="#d1495b", em="#e0a100", weak="#2e86ab", higgs="#6a4c93", grav="#555555", far="#bc4749")
def arrow(a, b, c, style="<->", ls="-", lw=2.2, rad=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=16, color=c, lw=lw, ls=ls, connectionstyle=f"arc3,rad={rad}"))
def spring(a, b, c, n=7, amp=0.12):
    a, b = np.array(a), np.array(b); t = np.linspace(0, 1, 200); d = b - a; L = np.hypot(*d); nrm = np.array([-d[1], d[0]])/L
    P = a[:, None] + d[:, None]*t + nrm[:, None]*amp*np.sin(2*np.pi*n*t); ax.plot(*P, color=c, lw=2)
# proton
ax.add_patch(Circle((3, 6), 1.55, fc="#fdecee", ec="#999", lw=1.2, ls="--")); ax.text(3, 7.75, "proton", ha="center", fontsize=12, weight="bold")
q = {"u1": (2.35, 6.45), "u2": (3.65, 6.45), "d": (3.0, 5.3)}
for k, (x, y) in q.items():
    ax.add_patch(Circle((x, y), 0.32*(1.25 if k == "d" else 1), fc="#ffd6dc" if k != "d" else "#cfe8f3", ec="k", lw=1.2)); ax.text(x, y, k[0], ha="center", va="center", fontsize=12, weight="bold")
for a, b in (("u1", "u2"), ("u2", "d"), ("d", "u1")): spring(q[a], q[b], COL["strong"])
ax.text(1.05, 4.15, "gluon springs (strong): pull harder\nthe farther apart, so quarks never\nget out. Pull ~ 15 tonnes-force at 1 fm", fontsize=9.5, color=COL["strong"])
# far-side press on d
arrow((4.7, 4.2), (3.3, 5.1), COL["far"], style="->", ls=":", lw=2)
ax.text(4.75, 3.7, "far-side press (TRPG):\nacts on down-type quarks only;\nreaches nowhere in space, it\nshows only as the down-quark\ntilt and stretch", fontsize=9, color=COL["far"])
# Higgs halos
for k, (x, y) in q.items(): ax.add_patch(Circle((x, y), 0.45*(1.4 if k == "d" else 1), fc="none", ec=COL["higgs"], lw=1.3, ls=(0, (2, 2))))
# electron
ax.add_patch(Circle((9.6, 7.0), 0.22, fc="#fff3c4", ec="k", lw=1.2)); ax.text(9.6, 7.0, "e", ha="center", va="center", fontsize=11, weight="bold")
ax.add_patch(Circle((9.6, 7.0), 0.36, fc="none", ec=COL["higgs"], lw=1.3, ls=(0, (2, 2))))
ax.text(9.6, 7.55, "electron", ha="center", fontsize=12, weight="bold")
arrow((4.6, 6.3), (9.3, 6.95), COL["em"], style="<->", rad=-0.15)
for i, x in enumerate(np.linspace(5.5, 8.4, 4)):
    t = np.linspace(0, 1, 60); ax.plot(x + 0.35*t, 6.9 + 0.08*np.sin(2*np.pi*3*t) + 0.06*i*0, color=COL["em"], lw=1.4)
ax.text(6.95, 5.75, "photons (electric): opposite charges pull,\nalike push. Holds the atom together.", ha="center", fontsize=9.5, color=COL["em"])
# gravity
arrow((3.0, 7.62), (9.6, 7.85), COL["grav"], style="<->", ls=(0, (1, 3)), lw=1.6, rad=-0.18)
ax.text(6.3, 8.75, "gravity: every mass and energy pulls on every other (TRPG: a byproduct of the\nentangled geometry). Here ~1e-39 of the electric pull: far too weak to draw.", ha="center", fontsize=9.2, color=COL["grav"])
# neutron beta decay
ax.add_patch(Circle((2.4, 1.6), 1.05, fc="#eef4f8", ec="#999", lw=1.2, ls="--")); ax.text(2.4, 2.85, "neutron  (udd)", ha="center", fontsize=11.5, weight="bold")
for (x, y, t) in ((2.0, 1.85, "u"), (2.8, 1.85, "d"), (2.4, 1.15, "d")):
    ax.add_patch(Circle((x, y), 0.25, fc="#cfe8f3" if t == "d" else "#ffd6dc", ec="k")); ax.text(x, y, t, ha="center", va="center", fontsize=11, weight="bold")
arrow((2.65, 1.05), (5.0, 1.1), COL["weak"], style="->", lw=2.5)
ax.text(3.6, 0.45, "W: turns this d into u (the weak turn\nbetween needles). Beta decay", fontsize=9.5, color=COL["weak"])
ax.text(5.1, 1.25, "→  proton + electron + antineutrino", fontsize=10.5, color=COL["weak"])
# neutrino and Z
ax.add_patch(Circle((10.3, 2.6), 0.15, fc="#e8e0f5", ec="k")); ax.text(10.3, 2.95, "neutrino", ha="center", fontsize=11, weight="bold")
arrow((10.2, 2.75), (9.65, 6.7), COL["weak"], style="<->", ls="--", lw=1.8)
ax.text(10.45, 4.4, "Z: a neutrino nudging\nthe electron without\nchanging either. Reach\nonly 0.0025 fm", fontsize=9.2, color=COL["weak"])
ax.text(0.2, 9.55, "Dashed purple halos = the Higgs hold on each particle (bigger = heavier). It gives mass; it does not push or pull anything here.",
        fontsize=9.6, color=COL["higgs"])
ax.text(0.2, 9.15, "NOT TO SCALE: the proton is ~1 fm wide and the electron cloud ~53,000 fm away (the atom would be a stadium with a pea in the middle).", fontsize=9.6, style="italic")
# right: table
tx = fig.add_axes([0.58, 0.06, 0.41, 0.86]); tx.axis("off"); tx.set_xlim(0, 1); tx.set_ylim(0, 1)
rows = [("Force", "Carried by", "Who feels it", "Push or pull?", "Strength\nat ~1 fm"),
        ("Strong", "8 gluons (they\nfeel it too)", "quarks, gluons", "pull, grows with\ndistance; binds", "1"),
        ("Electric\n(magnetic)", "photon\n(a light wave)", "anything charged:\nquarks, e, μ, τ, W", "pull or push;\nfades as 1/r²", "1/137"),
        ("Weak", "W+, W−, Z\n(heavy)", "every matter particle\n(neutrinos ONLY this\nand gravity)", "mostly turns kinds;\nreach 0.0025 fm", "~1e-5"),
        ("Higgs link", "Higgs boson", "every massive\nparticle", "gives mass; binds\nnothing", "∝ mass"),
        ("Gravity", "geometry itself\n(graviton unseen)", "everything with\nmass or energy", "always pull;\nfades as 1/r²", "~6e-39"),
        ("Far-side press\n(TRPG)", "the very heavy\nfar-side partner", "down-type quarks;\nneutrinos (own partner)", "bends needles\n(no reach in space)", "not a force\nin space")]
colx = [0.0, 0.17, 0.37, 0.62, 0.85]; cols = [COL["strong"], COL["em"], COL["weak"], COL["higgs"], COL["grav"], COL["far"]]
y = 0.97
for r, row in enumerate(rows):
    h = 0.075 if r == 0 else 0.095
    if r: tx.add_patch(FancyBboxPatch((0, y - h + 0.008), 1.0, h - 0.012, boxstyle="round,pad=0.003", fc=cols[r - 1], alpha=0.12, ec="none"))
    for c, txt in enumerate(row):
        tx.text(colx[c] + 0.005, y - h/2, txt, va="center", fontsize=10 if r else 10.5, weight="bold" if (r == 0 or c == 0) else None,
                color=cols[r - 1] if (r and c == 0) else "k")
    y -= h
tx.text(0.0, y - 0.02, "WHO TOUCHES WHOM\n"
        "• Gluons pull on each other (that is why the strong pull never lets go).\n"
        "• The photon does not feel itself; the W does feel the photon (it is charged).\n"
        "• The W and Z touch each other and the Higgs; the Higgs gives them their mass.\n"
        "• The neutrino feels only the weak force and gravity: it passes through Earth.\n"
        "• The far-side partner is too heavy to appear; it touches the down quarks through\n"
        "  the Higgs link only, and the neutrinos' own partner touches only them.\n\n"
        "Measured strengths: PDG 2024. Gravity and electric compared for two protons.\n"
        "TRPG main line: the 8|3|1 carriers come from the three-way matrix;\n"
        "the strengths themselves are the one item still open.", va="top", fontsize=10.2,
        bbox=dict(boxstyle="round,pad=0.5", fc="#fffbe6", ec="#F5C518", lw=2))
out = "TRPG_Forces_Map.png"; plt.savefig(out, dpi=160); print("saved", out)
