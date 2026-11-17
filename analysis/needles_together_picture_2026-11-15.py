"""Picture 1: the non-object drawn as an object. Quark pair and lepton pair, heavy ends up, light ends down, fat (log-drawn)
needles, the turn between partners, the twist ribbon, the far-side presses (down quarks; neutrinos' own hidden partner), the stretch before
and after, and the six faces of the quark turn at their true shapes. Honest key."""
import numpy as np, math, matplotlib, itertools
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
def semi(m): m = np.array(m, float); return 1 + np.log10(m/m[0])
uu, vv = np.meshgrid(np.linspace(0, 2*np.pi, 50), np.linspace(0, np.pi, 26))
sph = np.stack([np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv)])
def ell(a, R=np.eye(3), S=np.eye(3), off=np.zeros(3)):
    X = np.einsum("ij,jkl->ikl", S @ R @ np.diag(a), sph); return X + off[:, None, None]
def Rx(t): c, s = math.cos(t), math.sin(t); return np.array([[1,0,0],[0,c,-s],[0,s,c]])
def Ry(t): c, s = math.cos(t), math.sin(t); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rz(t): c, s = math.cos(t), math.sin(t); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
u = np.array([-1, -1, 2])/math.sqrt(6)
fig = plt.figure(figsize=(18, 12.5))
fig.suptitle("The hidden geometry drawn as if it were an object: two pairs of needles, their turn, twist, presses and faces", fontsize=16, weight="bold", y=0.985)
def pair(ax, aA, aB, colA, colB, R, title, press=None, press_label="", twist=True):
    ax.plot_wireframe(*ell(aA), color=colA, alpha=0.3, rstride=3, cstride=4)
    if press is not None:
        S = np.eye(3) - press*np.outer(u, u)
        ax.plot_wireframe(*ell(aB, R), color="gray", alpha=0.25, rstride=3, cstride=4)
        ax.plot_surface(*ell(aB, R, S), color=colB, alpha=0.5, linewidth=0)
        for sgn in (1, -1):
            p0 = sgn*7.3*u; d = -sgn*2.0*u; ax.quiver(*p0, *d, color="#bc4749", lw=3, arrow_length_ratio=0.35)
    else:
        ax.plot_surface(*ell(aB, R), color=colB, alpha=0.5, linewidth=0)
    for i in range(3):
        e = np.eye(3)[i]*7.0; ax.plot(*np.stack([-e, e]).T, color=colA, lw=1.8)
        f = R @ (np.eye(3)[i]*7.0); ax.plot(*np.stack([-f, f]).T, color=colB, lw=1.8, ls="--")
    if twist:
        t = np.linspace(-1, 1, 300); ang = 2*np.pi*2*t; r = 1.15*np.sqrt(np.clip(1 - t**2, 0.05, 1))
        H = R @ np.stack([r*np.cos(ang), r*np.sin(ang), aB[2]*0.95*t]); ax.plot(*H, color="#7b2cbf", lw=2.3)
    ax.text(0, 0, 7.6, "HEAVY end (third family)", ha="center", fontsize=10, weight="bold")
    ax.text(0, 0, -8.4, "light end (first family)", ha="center", fontsize=10)
    ax.set_xlim(-7.5, 7.5); ax.set_ylim(-7.5, 7.5); ax.set_zlim(-8.5, 8.5); ax.set_box_aspect((1, 1, 1.15), zoom=1.3)
    ax.set_axis_off(); ax.view_init(elev=12, azim=-58); ax.set_title(title, fontsize=12.5, pad=0)
aU, aD = semi((2.16, 1273, 172570)), semi((4.70, 93.5, 4183))
aE, aN = semi((0.511, 105.66, 1776.9)), semi((0.00036, 0.00866, 0.0503))
Rq = Rx(math.radians(2.4*5)) @ Ry(math.radians(0.21*20)) @ Rz(math.radians(13.0))
Rl = Rx(math.radians(49.0)) @ Ry(math.radians(8.6)) @ Rz(math.radians(33.4))
ax1 = fig.add_axes([0.0, 0.30, 0.36, 0.64], projection="3d")
pair(ax1, aU, aD, "#d1495b", "#2e86ab", Rq, "Quarks: up (red mesh) and down (blue)", press=0.25, press_label="far-side press on the\ndown needle (gray = before,\nblue = after: the stretch)")
fig.text(0.18, 0.375, "Red arrows: the far-side press on the down needle. Gray ghost = before the press,\nsolid blue = after: the stretch. Dashed axes = the down needle's turned directions.", ha="center", color="#bc4749", fontsize=10)
ax2 = fig.add_axes([0.33, 0.30, 0.36, 0.64], projection="3d")
pair(ax2, aE, aN, "#3b8b3b", "#8e6cc2", Rl, "Leptons: electron family (green mesh) and neutrinos (purple)", press=0.25, press_label="the neutrinos' own hidden\npartner presses them\n(why they are so light)")
fig.text(0.51, 0.375, "Red arrows: the neutrinos' own hidden partner presses them (why they are so light).\nGray ghost = before, solid purple = after. Dashed axes = the neutrinos' turned directions.", ha="center", color="#bc4749", fontsize=10)
# six faces of the quark turn, true shapes
lam, A, rb, eb = 0.22501, 0.826, 0.1591, 0.3523
s12 = lam; s23 = A*lam**2; z = A*lam**3*(rb + 1j*eb)*math.sqrt(1 - A**2*lam**4)/(math.sqrt(1 - lam**2)*(1 - A**2*lam**4*(rb + 1j*eb)))
s13 = abs(z); d = -np.angle(z); c12, c23, c13 = (math.sqrt(1 - x*x) for x in (s12, s23, s13)); e = np.exp(1j*d)
V = np.array([[c12*c13, s12*c13, s13/e], [-s12*c23 - c12*s23*s13*e, c12*c23 - s12*s23*s13*e, s23*c13],
              [s12*s23 - c12*c23*s13*e, -c12*s23 - s12*c23*s13*e, c23*c13]])
axT = fig.add_axes([0.02, 0.03, 0.66, 0.24]); axT.set_xlim(0, 6.6); axT.set_ylim(-0.45, 1.25); axT.axis("off"); axT.set_aspect("equal")
axT.text(0, 1.2, "The six faces of the quark turn, at their TRUE shapes (each scaled so its longest side is the same length). Their true areas are all equal (J/2 = 1.6e-5).",
         fontsize=10.5, va="bottom")
names = []
for idx, (kind, (i, j)) in enumerate(itertools.product(("columns", "rows"), itertools.combinations(range(3), 2))):
    sides = [V[k, i]*np.conj(V[k, j]) for k in range(3)] if kind == "columns" else [V[i, k]*np.conj(V[j, k]) for k in range(3)]
    P = np.array([0, sides[0], sides[0] + sides[1]]); L = max(abs(s) for s in sides)
    order = np.argsort([abs(s) for s in sides]); P = P/L
    # rotate so the longest side lies flat
    longest = int(np.argmax([abs(s) for s in sides])); a0 = P[longest]; a1 = P[(longest + 1) % 3]
    P = (P - a0)*np.exp(-1j*np.angle(a1 - a0)); P = P - P.real.min()
    if P.imag.mean() < 0: P = np.conj(P)
    x0 = 0.1 + idx*1.08; poly = np.array([[p.real*0.95 + x0, p.imag*0.95] for p in P])
    axT.fill(poly[:, 0], poly[:, 1], color="#F5C518", alpha=0.6, ec="#715800", lw=1.5)
    axT.text(x0 + 0.47, -0.17, f"{kind} {i+1}-{j+1}", ha="center", fontsize=9.5)
axT.text(0.1, -0.42, "Four of the six are so thin they look like lines: that IS their true shape (one side tiny, two nearly equal). The two fat ones (1-3) are the ones experiments measure best.", fontsize=9.5, style="italic")
kx = fig.add_axes([0.70, 0.04, 0.29, 0.9]); kx.axis("off")
key = ("HOW IT FITS TOGETHER\n\n"
       "• Each kind of particle is one needle. Its three\n"
       "  lengths are its three masses: heavy end up,\n"
       "  light end down.\n"
       "• Partner needles share one centre and are\n"
       "  turned against each other. The turn is how a\n"
       "  particle changes kind (the weak force).\n"
       "• The twist (purple coil) is a hidden timing on\n"
       "  the turn. Its size shows up as the area of the\n"
       "  six faces below, which are all equal.\n"
       "• The far side presses on one needle of each\n"
       "  pair at a slant (red arrows). The press bends\n"
       "  the needle's axes a little: that bend sets the\n"
       "  tilts. For neutrinos the press is why they\n"
       "  are so light.\n"
       "• All four needles live in the same three family\n"
       "  directions: the one relationship matrix.\n\n"
       "HONEST SCALE\n"
       "To scale (log rule): order and relative lengths\n"
       "  of each needle's axes; the 13°, 33° and 49°\n"
       "  turns; the faces' shapes.\n"
       "Exaggerated: needle widths (true ones are\n"
       "  lines), the quark turns of 2.4° (×5) and 0.21°\n"
       "  (×20), the press size (drawn 25%; about 12%).\n"
       "Cannot be drawn: family directions are not\n"
       "  space; the twist is a hidden angle; the far\n"
       "  side is a fourth, very heavy direction; the\n"
       "  faces live in the complex plane of the turn.\n\n"
       "Measured: masses, turns, faces (PDG 2024,\n"
       "NuFIT). TRPG: neutrino masses, the press.\n"
       "earlier tests")
kx.text(0.0, 1.0, key, va="top", ha="left", fontsize=10.5, bbox=dict(boxstyle="round,pad=0.6", fc="#fffbe6", ec="#F5C518", lw=2))
out = "TRPG_Needles_Together.png"; plt.savefig(out, dpi=160); print("saved", out)
