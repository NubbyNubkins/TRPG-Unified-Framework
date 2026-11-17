"""Picture: the zoom map. Six levels from the full glass down to quarks, each with its size, the force that rules it,
and up (red) / down (blue) quarks marked throughout; neutrinos, W/Z, Higgs hold and far side shown where they act."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, ConnectionPatch
rng = np.random.default_rng(7)
U, D, E, G, N, W, H, F = "#d1495b", "#2e86ab", "#e0a100", "#c0392b", "#8e6cc2", "#2e86ab", "#6a4c93", "#bc4749"
fig = plt.figure(figsize=(20, 13))
fig.suptitle("Zooming in: how everything is built from up and down quarks (red = up, blue = down), and which force rules each level",
             fontsize=16, weight="bold", y=0.985)
pos = [(0.01, 0.50), (0.34, 0.50), (0.67, 0.50), (0.67, 0.05), (0.34, 0.05), (0.01, 0.05)]
axs = []
for (x, y) in pos:
    a = fig.add_axes([x, y, 0.31, 0.42]); a.set_xlim(-1.25, 1.25); a.set_ylim(-1.25, 1.25); a.set_aspect("equal"); a.axis("off"); axs.append(a)
def frame(a, n, title, size, force):
    a.add_patch(plt.Rectangle((-1.24, -1.24), 2.48, 2.48, fill=False, ec="#bbb", lw=1))
    a.text(-1.2, 1.16, f"{n}. {title}", fontsize=12.5, weight="bold", va="top")
    a.text(-1.2, 1.0, f"size: {size}", fontsize=10, va="top", color="#444")
    a.text(-1.2, -1.2, force, fontsize=9.6, va="bottom", color="#222", wrap=True)
def nucleon(a, x, y, r, kind, lw=0.8, label=False):
    a.add_patch(Circle((x, y), r, fc="#fdecee" if kind == "p" else "#e7f1f7", ec="#888", lw=lw))
    cols = [U, U, D] if kind == "p" else [U, D, D]
    for k, c in enumerate(cols):
        t = np.pi/2 + 2*np.pi*k/3; a.add_patch(Circle((x + 0.48*r*np.cos(t), y + 0.48*r*np.sin(t)), 0.24*r, fc=c, ec="none"))
    if label: a.text(x, y - r - 0.08, "proton (uud)" if kind == "p" else "neutron (udd)", ha="center", va="top", fontsize=10)
# 1 the full glass
a = axs[0]; frame(a, 1, "The full glass (what we can see of it)", "about 10^26 m", "Rules here: gravity, and the kept amount of space (dark energy).\nTRPG: one full glass, no gaps; the cosmic web is the drizzle's ripples.")
for _ in range(14):
    x0, y0 = rng.uniform(-0.9, 0.9, 2); t = np.linspace(0, 1, 30); a.plot(x0 + 0.35*rng.normal()*t, y0 + 0.35*rng.normal()*t, color="#9fb3c8", lw=1.2)
pts = rng.uniform(-0.95, 0.85, (260, 2)); a.scatter(pts[:, 0], pts[:, 1], s=3, color="#6b7d91")
a.add_patch(Circle((0.45, -0.25), 0.07, fc="none", ec="k", lw=1.5)); a.text(0.55, -0.15, "a galaxy;\nzoom into one\ndead star in it", fontsize=9)
# 2 neutron star
a = axs[1]; frame(a, 2, "Neutron star", "about 24 km across", "Rules here: gravity crushes atoms; the weak force turns\nproton + electron into neutron + neutrino. TRPG prediction:\nquark pairing in the core caps the mass at 2.0-2.15 Suns (measured 2.07).")
a.add_patch(Circle((0, 0), 0.82, fc="#e7f1f7", ec="#666", lw=1.5)); a.add_patch(Circle((0, 0), 0.42, fc="#f3e6f7", ec="#8e6cc2", lw=1.2, ls="--"))
a.text(0, 0, "core: quarks\nmay pair up\n(u with d)", ha="center", va="center", fontsize=9.5)
for r in np.linspace(0.5, 0.78, 3):
    for t in np.linspace(0, 2*np.pi, int(30*r), endpoint=False): nucleon(a, r*np.cos(t), r*np.sin(t), 0.035, "n" if rng.random() < 0.92 else "p")
a.text(0.86, 0.55, "mostly neutrons:\nthe down quarks\nwin when crushed", fontsize=9, color=D)
for t in (0.5, 1.9, 2.9, 3.7, 5.9):
    a.annotate("", xy=(1.15*np.cos(t), 1.15*np.sin(t)), xytext=(0.8*np.cos(t), 0.8*np.sin(t)), arrowprops=dict(arrowstyle="->", color=N, lw=1.2))
a.text(-1.15, 0.78, "neutrinos stream out\n(they feel only weak + gravity)", fontsize=8.8, color=N)
# 3 atom
a = axs[2]; frame(a, 3, "Atom (carbon)", "about 0.15 nanometres (150,000 fm)", "Rules here: the electric pull (photons) holds electrons to the nucleus.\nThe nucleus would be a grain of sand in a stadium: drawn about 2,000x too big.")
for rr in (0.45, 0.85): a.add_patch(Circle((0, 0), rr, fc="none", ec="#e9c46a", lw=8, alpha=0.35))
for t in np.linspace(0, 2*np.pi, 6, endpoint=False):
    rr = 0.45 if t < 2 else 0.85; a.add_patch(Circle((rr*np.cos(t + 0.3), rr*np.sin(t + 0.3)), 0.05, fc=E, ec="k", lw=0.6))
a.add_patch(Circle((0, 0), 0.12, fc="#f5d0d6", ec="k", lw=0.8)); a.text(0.15, 0.05, "nucleus", fontsize=9)
a.text(0.62, -0.95, "electrons (yellow): feel electric,\nweak and gravity; no colour", fontsize=9, color="#9a6b00")
for k in range(4):
    x = -0.9 + 0.12*k; a.plot(np.linspace(x, x + 0.25, 40), 0.75 + 0.03*np.sin(np.linspace(0, 12, 40)), color=E, lw=1)
a.text(-1.15, 0.86, "photons", fontsize=9, color=E)
# 4 nucleus
a = axs[3]; frame(a, 4, "Nucleus (carbon-12)", "about 5 fm (0.000005 nm)", "Rules here: the leftover strong pull between nucleons (pions:\na u and anti-d passed back and forth) beats the electric push\nbetween protons. The weak force can turn a d into u (beta decay).")
xy = [(0, 0)] + [(0.33*np.cos(t), 0.33*np.sin(t)) for t in np.linspace(0, 2*np.pi, 6, endpoint=False)] + [(0.62*np.cos(t), 0.62*np.sin(t)) for t in np.linspace(0.5, 2*np.pi + 0.5, 5, endpoint=False)]
for i, (x, y) in enumerate(xy): nucleon(a, x, y, 0.16, "p" if i % 2 == 0 else "n")
a.annotate("", xy=(0.62*np.cos(0.5) - 0.05, 0.62*np.sin(0.5)), xytext=(0.33 + 0.12, 0.05), arrowprops=dict(arrowstyle="<->", color="#555", lw=1.5, ls="--"))
a.text(0.55, 0.05, "pion = u + anti-d:\nthe nuclear glue", fontsize=9, color="#555")
# 5 proton and neutron
a = axs[4]; frame(a, 5, "Proton and neutron", "about 1.7 fm across", "Rules here: the strong force (gluons) holds 3 quarks.\n99% of the mass is this binding energy, not the quarks'\nHiggs mass (TRPG: mass = balance of containment and leakage).")
nucleon(a, -0.52, 0.12, 0.42, "p", lw=1.4, label=True); nucleon(a, 0.52, 0.12, 0.42, "n", lw=1.4, label=True)
for cx, kind in ((-0.52, "p"), (0.52, "n")):
    for k in range(3):
        t1, t2 = np.pi/2 + 2*np.pi*k/3, np.pi/2 + 2*np.pi*(k + 1)/3
        p1 = np.array([cx + 0.2*np.cos(t1), 0.12 + 0.2*np.sin(t1)]); p2 = np.array([cx + 0.2*np.cos(t2), 0.12 + 0.2*np.sin(t2)])
        tt = np.linspace(0, 1, 80); d = p2 - p1; nrm = np.array([-d[1], d[0]])/np.hypot(*d)
        P = p1[:, None] + d[:, None]*tt + nrm[:, None]*0.025*np.sin(2*np.pi*6*tt); a.plot(*P, color=G, lw=1.3)
a.text(-0.95, 0.75, "proton: 2 up + 1 down,\ncharge +1, lasts forever", fontsize=9, color=U)
a.text(0.2, 0.75, "neutron: 1 up + 2 down,\ncharge 0, 1.3 MeV heavier", fontsize=9, color=D)
a.annotate("", xy=(-0.12, -0.42), xytext=(0.12, -0.42), arrowprops=dict(arrowstyle="->", color=W, lw=2.2))
a.text(0, -0.5, "W turns one d into u: neutron becomes proton\n(+ electron + antineutrino) when free", ha="center", va="top", fontsize=9, color=W)
# 6 quarks
a = axs[5]; frame(a, 6, "Quarks, gluons and their links", "under 10^-19 m (points)", "Rules here: every force meets the quarks: gluons (colour),\nphotons (charge), W and Z (weak pair), Higgs hold (mass).\nTRPG: the far side presses DOWN quarks only (both are 'singles').")
a.add_patch(Circle((-0.45, 0.15), 0.16, fc=U, ec="k")); a.text(-0.45, 0.15, "u", ha="center", va="center", color="w", fontsize=14, weight="bold")
a.add_patch(Circle((0.45, 0.15), 0.19, fc=D, ec="k")); a.text(0.45, 0.15, "d", ha="center", va="center", color="w", fontsize=14, weight="bold")
for (cx, r) in ((-0.45, 0.28), (0.45, 0.34)): a.add_patch(Circle((cx, 0.15), r, fc="none", ec=H, ls=(0, (2, 2)), lw=1.3))
tt = np.linspace(-0.27, 0.24, 120); a.plot(tt, 0.32 + 0.03*np.sin(60*tt), color=G, lw=1.8); a.text(0, 0.42, "gluon", ha="center", fontsize=9, color=G)
a.annotate("", xy=(0.24, 0.02), xytext=(-0.26, 0.02), arrowprops=dict(arrowstyle="<->", color=W, lw=1.8)); a.text(0, -0.08, "W (turns u <-> d)", ha="center", va="top", fontsize=9, color=W)
a.text(-0.95, -0.38, "dashed purple ring = Higgs hold\n(gives each its mass; d's is bigger)", fontsize=8.8, color=H)
a.add_patch(Circle((0.95, 0.75), 0.12, fc="#f6d6d6", ec=F, ls=":", lw=1.5)); a.text(0.95, 0.75, "D", ha="center", va="center", fontsize=11, color=F, weight="bold")
a.annotate("", xy=(0.55, 0.32), xytext=(0.86, 0.66), arrowprops=dict(arrowstyle="->", color=F, lw=1.6, ls=":"))
a.text(0.3, 0.93, "far side (TRPG): very heavy, never seen;\npresses d only, not u", fontsize=8.8, color=F)
a.text(0.18, -0.4, "also here: photon (charge;\nu feels it 4x more than d), Z", fontsize=8.8, color="#555")
# zoom arrows between panels
order = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
for i, j in order:
    xa, ya = pos[i]; xb, yb = pos[j]
    if ya == yb and xb > xa: s, e = (xa + 0.31, ya + 0.21), (xb, yb + 0.21)
    elif ya == yb: s, e = (xa, ya + 0.21), (xb + 0.31, yb + 0.21)
    else: s, e = (xa + 0.155, ya), (xb + 0.155, yb + 0.42)
    if ya != yb or True:
        fig.patches.append(FancyArrowPatch(s, e, transform=fig.transFigure, arrowstyle="-|>", mutation_scale=22, color="#F5C518", lw=3))
fig.text(0.5, 0.475, "zoom in  (each step is 1,000 to 10^20 times smaller)", ha="center", fontsize=10, color="#715800") if False else None
fig.text(0.995, 0.478, "zoom in \u2193", ha="right", fontsize=11, color="#715800", weight="bold")
fig.text(0.01, 0.006, "Read 1 \u2192 2 \u2192 3, then down and right-to-left 4 \u2192 5 \u2192 6. NOT TO SCALE between panels (the sizes span 45 powers of ten); inside panels 3-6 the layout is schematic.\n"
         "Measured: sizes, charges, which force acts where (PDG 2024; NICER/LIGO for the neutron star). TRPG: the far-side press, the pairing cap, the readings in the captions.", fontsize=9.6, style="italic")
out = "TRPG_Zoom_Map.png"; plt.savefig(out, dpi=150); print("saved", out)
