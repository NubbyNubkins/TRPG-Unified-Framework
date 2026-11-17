"""Picture: one view with everything described: the up and down needles (heavy axis up), the turn, the twist,
handedness as the side you observe, the Higgs vibration, a particle's swing between views, and the far-side press. Honest key."""
import numpy as np, math, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
def semi(m): m = np.array(m, float); return 1 + np.log10(m/m[0])
uu, vv = np.meshgrid(np.linspace(0, 2*np.pi, 60), np.linspace(0, np.pi, 30))
sph = np.stack([np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv)])
def ell(a, R=np.eye(3), S=np.eye(3), f=1.0): return np.einsum("ij,jkl->ikl", S @ R @ np.diag(np.array(a)*f), sph)
def Rx(t): c, s = math.cos(t), math.sin(t); return np.array([[1,0,0],[0,c,-s],[0,s,c]])
def Ry(t): c, s = math.cos(t), math.sin(t); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rz(t): c, s = math.cos(t), math.sin(t); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
aU = semi((2.16, 1273, 172570)); aD = semi((4.70, 93.5, 4183))
R = Rx(math.radians(2.4*5)) @ Ry(math.radians(0.21*20)) @ Rz(math.radians(13.0))
u = np.array([-1, -1, 2])/math.sqrt(6); S = np.eye(3) - 0.12*np.outer(u, u)
fig = plt.figure(figsize=(17, 11))
fig.suptitle("The whole hidden geometry in one view (up and down quarks shown; leptons work the same way)", fontsize=16, weight="bold", y=0.975)
ax = fig.add_axes([-0.02, 0.0, 0.66, 0.95], projection="3d")
# Higgs vibration: the up needle breathing in proportion
for f in (0.9, 1.1): ax.plot_wireframe(*ell(aU, f=f), color="#d1495b", alpha=0.10, rstride=4, cstride=6, linestyle="--")
ax.plot_wireframe(*ell(aU), color="#d1495b", alpha=0.35, rstride=3, cstride=4)
ax.plot_surface(*ell(aD, R, S), color="#2e86ab", alpha=0.45, linewidth=0)
for i in range(3):
    e = np.eye(3)[i]*7.2; ax.plot(*np.stack([-e, e]).T, color="#d1495b", lw=2)
    f_ = R @ (np.eye(3)[i]*7.2); ax.plot(*np.stack([-f_, f_]).T, color="#2e86ab", lw=2, ls="--")
ax.text(0, 0, 7.6, "third family (heavy)", color="#333", fontsize=10, ha="center")
ax.text(7.6, 0, 0, "first family\n(light)", color="#333", fontsize=9.5)
ax.text(-1.5, 8.0, -0.8, "second\nfamily", color="#333", fontsize=9.5)
# twist: helix around the down needle's long axis
t = np.linspace(-1, 1, 400); ang = 2*np.pi*2.2*t
r = 1.25*np.sqrt(np.clip(1 - t**2, 0.05, 1))
H = np.stack([r*np.cos(ang), r*np.sin(ang), aD[2]*t*0.95]); H = R @ H
ax.plot(*H, color="#7b2cbf", lw=2.5)
ax.text(*(R @ np.array([1.6, 0.4, 2.6])), "twist (phase)", color="#7b2cbf", fontsize=10.5, weight="bold")
# handedness ring with viewers
th = np.linspace(0, 2*np.pi*0.9, 200); z0 = 0.0; rr = 3.4
ax.plot(rr*np.cos(th), rr*np.sin(th), z0*np.ones_like(th), color="#e09f3e", lw=2.5)
ax.quiver(rr*np.cos(th[-1]), rr*np.sin(th[-1]), z0, -np.sin(th[-1])*1.2, np.cos(th[-1])*1.2, 0, color="#e09f3e", lw=2.5, arrow_length_ratio=0.5)
ax.text(0, 0, 9.0, "seen from above: counterclockwise", color="#a5650f", fontsize=10, ha="center")
ax.text(0, 0, -9.0, "seen from below: clockwise", color="#a5650f", fontsize=10, ha="center")
# a particle's swing between views
ph = np.linspace(0, np.pi/2, 50); sw = np.stack([4.6*np.cos(ph), 4.6*np.sin(ph), 1.2*np.sin(2*ph)])
ax.plot(*sw, color="#2a9d8f", lw=2.5); ax.quiver(*sw[:, -2], *(sw[:, -1] - sw[:, -2])*6, color="#2a9d8f", lw=2.5, arrow_length_ratio=0.8)
ax.text(5.2, 1.0, 3.0, "a particle swings\nbetween views", color="#2a9d8f", fontsize=10)
# far-side press
for sgn in (1, -1):
    p0 = sgn*7.4*u; d = -sgn*2.0*u; ax.quiver(*p0, *d, color="#bc4749", lw=3, arrow_length_ratio=0.35)
ax.text(*(7.8*u + np.array([0.3, 0, 0.4])), "far-side press\n(down only)", color="#bc4749", fontsize=10)
lim = 8; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim - 1.5, lim + 1.5); ax.set_box_aspect((1, 1, 1.15), zoom=1.45)
ax.set_axis_off(); ax.view_init(elev=14, azim=-55)
kx = fig.add_axes([0.63, 0.05, 0.36, 0.80]); kx.axis("off")
key = ("WHAT EACH PART IS\n\n"
       "Red mesh / red lines: up-quark needle and its axes\n"
       "Blue / blue dashed: down-quark needle, turned against it\n"
       "  (13.0° true; 2.4° ×5 and 0.21° ×20 so they show)\n"
       "Faint dashed shells: the Higgs vibration. The Higgs boson\n"
       "  makes every axis grow and shrink together, in proportion\n"
       "  (each particle's link to it is proportional to its mass)\n"
       "Purple helix: the twist (the turn's phase, about 66°)\n"
       "Orange ring: one turning seen from above and below. Same\n"
       "  motion, opposite names: handedness is the side you see\n"
       "Green arrow: a particle's state swinging between views\n"
       "  (meson and neutrino oscillation); the needle stays fixed\n"
       "Dark red arrows: the far-side press, on the down needle only\n\n"
       "HONEST SCALE\n"
       "To scale (log rule): order and relative axis lengths,\n"
       "  direction of the press, the 13° turn\n"
       "Exaggerated: needle widths (true up needle is 130,800 times\n"
       "  longer than wide), small turns, vibration size, press size\n"
       "Cannot be drawn: axes are families, not space; the twist is a\n"
       "  hidden angle, not a real coil; the far side is a fourth,\n"
       "  very heavy direction; the ring is a picture of the idea\n\n"
       "Measured: masses, turns, twist (PDG 2024), oscillation rates,\n"
       "Higgs links (ATLAS & CMS 2022). TRPG: far-side press, one i\n"
       "throughout (rule 4).")
kx.text(0.0, 1.0, key, va="top", ha="left", fontsize=10.6, bbox=dict(boxstyle="round,pad=0.6", fc="#fffbe6", ec="#F5C518", lw=2))
out = "TRPG_Needle_Full_Picture.png"
plt.savefig(out, dpi=170); print("saved", out)
