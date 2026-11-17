"""Up quark vs down quark, every characteristic, two rows (split into three column blocks)."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
hc, v, s2w = 197.327, 246.22, 0.2312
mu, md = 2.16, 4.70
y = lambda m: 2**0.5*m/v/1000
print(f"Compton u {hc/mu:.1f} fm, d {hc/md:.1f} fm; Higgs link u {y(mu):.2e}, d {y(md):.2e}; c {y(1273):.2e} s {y(93.5):.2e} t {y(172570):.3f} b {y(4183):.4f}")
print(f"Z vector coupling u {0.5-4/3*s2w:+.4f}, d {-0.5+2/3*s2w:+.4f}; ladders top/up {172570/2.16:,.0f} bottom/down {4183/4.70:,.0f}")
blocks = [
 [("", "Up quark (u)", "Down quark (d)"),
  ("Mass (weighed at 2 GeV, PDG 2024)", "2.16 MeV (lighter)", "4.70 MeV (about 2.2x heavier)"),
  ("Electric charge", "+2/3", "-1/3"),
  ("Spin", "1/2", "1/2  (same)"),
  ("Colour (strong charge)", "red, green or blue", "red, green or blue  (same)"),
  ("Weak role (left side)", "upper member of the pair", "lower member of the pair")],
 [("", "Up quark (u)", "Down quark (d)"),
  ("Hypercharge (left side / right side)", "+1/6  /  +2/3", "+1/6  /  -1/3"),
  ("Higgs link (strength of the hold)", "1.2e-5", "2.7e-5"),
  ("Physical size", "point-like: under 1e-19 m", "point-like: under 1e-19 m  (same)"),
  ("Quantum wave size (h-bar / mc)", "91 fm (bigger wave)", "42 fm"),
  ("Made of", "nothing known inside", "nothing known inside  (same)")],
 [("", "Up quark (u)", "Down quark (d)"),
  ("Stability inside a nucleon", "stays (proton > 2.4e34 years)", "turns into up (free neutron 878 s)"),
  ("Family members", "charm 1.27 GeV, top 172.6 GeV", "strange 93.5 MeV, bottom 4.18 GeV"),
  ("Ladder span (heaviest / lightest)", "about 80,000", "about 890"),
  ("Right-side partner (TRPG 78.2)", "a PAIR of colour directions", "a SINGLE direction"),
  ("Far side (TRPG 68)", "no heavy partner: no press", "heavy partner D: the fold presses it")],
]
rows = [r for blk in blocks for r in blk[1:]]
rows.insert(9, ("Z-boson link (vector part)", "+0.19", "-0.35"))
rows.insert(10, ("Electric pull felt (charge squared)", "4/9  (4x stronger)", "1/9"))
rows.insert(11, ("Mass when dressed inside a proton", "about 336 MeV", "about 340 MeV  (same)"))
fig = plt.figure(figsize=(13, 11.5)); ax = fig.add_axes([0.02, 0.1, 0.96, 0.84]); ax.axis("off")
fig.suptitle("Up quark vs down quark: what is the same and what changes", fontsize=17, weight="bold", y=0.985)
N = len(rows); ax.set_xlim(0, 1); ax.set_ylim(0, N + 1)
for j, h in enumerate(("Characteristic", "Up quark (u)", "Down quark (d)")):
    ax.text((0.01, 0.40, 0.70)[j], N + 0.5, h, fontsize=13, weight="bold", va="center", color=("k", "#b02e3c", "#1f6f94")[j])
for i, (lab, u, d) in enumerate(rows):
    yy = N - i - 0.5; same = "(same)" in d
    ax.add_patch(plt.Rectangle((0, yy - 0.47), 1, 0.94, fc="#e6f2e6" if same else "#fff1c9", ec="white"))
    ax.text(0.01, yy, lab, va="center", fontsize=11, weight="bold")
    ax.text(0.40, yy, u, va="center", fontsize=11, color="#b02e3c")
    ax.text(0.70, yy, d.replace("  (same)", "").replace(" (same)", ""), va="center", fontsize=11, color="#1f6f94")
fig.text(0.02, 0.07, "Green: the same for both.  Yellow: what changes. Only four things really change: charge, weak direction, Higgs hold (mass),\nand in TRPG the right-side partner's shape (pair vs single) and the far-side press. Everything else follows from those.", fontsize=11)
fig.text(0.02, 0.02, "First family is reversed (up lighter than down); in the heavier families up-type is heavier (charm > strange, top > bottom).\nThe reversal makes the neutron heavier than the proton, so protons last and hydrogen exists. Sources: PDG 2024; TRPG Guide 68, 78, 79.", fontsize=11, style="italic")
out = "TRPG_Up_vs_Down.png"; plt.savefig(out, dpi=150); print("saved", out)
