"""Picture: the hidden 'needles' (mass matrices as ellipsoids in family space), drawn with honest labels."""
import numpy as np, math, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
def semi(m): m = np.array(m, float); return 1 + np.log10(m/m[0])   # log-scaled semi-axes (shortest = 1)
kinds = [("Up quarks", (1.242, 627.8, 162500), "#d1495b"), ("Down quarks", (2.678, 53.36, 2873), "#2e86ab"),
         ("Charged leptons", (0.511, 105.66, 1776.9), "#3b8b3b"), ("Neutrinos (TRPG prediction)", (0.00036, 0.00866, 0.05034), "#8e6cc2")]
uu, vv = np.meshgrid(np.linspace(0, 2*np.pi, 60), np.linspace(0, np.pi, 30))
sph = np.stack([np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv)])
def ell(a, R=np.eye(3), S=np.eye(3)):
    X = np.einsum("ij,jkl->ikl", S @ R @ np.diag(a), sph); return X
def Rx(t): c, s = math.cos(t), math.sin(t); return np.array([[1,0,0],[0,c,-s],[0,s,c]])
def Ry(t): c, s = math.cos(t), math.sin(t); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rz(t): c, s = math.cos(t), math.sin(t); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def setup(ax, title, lim=6.5):
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim); ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("family 1", labelpad=-12); ax.set_ylabel("family 2", labelpad=-12); ax.set_zlabel("family 3", labelpad=-12); ax.set_title(title, fontsize=12, pad=4)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.view_init(elev=18, azim=-58)
fig = plt.figure(figsize=(16, 11.5))
fig.suptitle("The hidden geometry: each kind of particle is a 'needle' in family space", fontsize=17, weight="bold", y=0.985)
# Panel A: four needles
for idx, (name, m, col) in enumerate(kinds):
    ax = fig.add_subplot(2, 4, idx+1, projection="3d"); a = semi(m)
    X = ell(a); ax.plot_surface(*X, color=col, alpha=0.55, linewidth=0, shade=True)
    setup(ax, name)
    r = np.array(m)/m[0]
    ax.text2D(0.5, 0.93, f"three masses 1 : {r[1]:,.0f} : {r[2]:,.0f}", transform=ax.transAxes, fontsize=10.5, ha="center")
# Panel B: up and down turned against each other
ax = fig.add_subplot(2, 3, 4, projection="3d")
R = Rx(math.radians(2.40*5)) @ Ry(math.radians(0.214*20)) @ Rz(math.radians(13.0))
ax.plot_wireframe(*ell(semi(kinds[0][1])), color=kinds[0][2], alpha=0.35, rstride=3, cstride=4)
ax.plot_surface(*ell(semi(kinds[1][1]), R), color=kinds[1][2], alpha=0.45, linewidth=0)
for i in range(3):
    e = np.eye(3)[i]*7.0; ax.plot(*np.stack([-e, e]).T, color=kinds[0][2], lw=2.2)
    f = R @ (np.eye(3)[i]*7.0); ax.plot(*np.stack([-f, f]).T, color=kinds[1][2], lw=2.2, ls="--")
setup(ax, "Up (red) and down (blue) turned against each other")
ax.text2D(0.0, -0.09, "Lines = each needle's three axes (red solid = up, blue dashed = down).\nTurn: 13.0° true; 2.4° drawn ×5; 0.21° drawn ×20.\nThe turn is how particles change kind.", transform=ax.transAxes, fontsize=10)
# Panel C: far-side squeeze on down
ax = fig.add_subplot(2, 3, 5, projection="3d")
u = np.array([-1, -1, 2])/math.sqrt(6); S = np.eye(3) - 0.12*np.outer(u, u)
ad = semi(kinds[1][1])
ax.plot_wireframe(*ell(ad), color="gray", alpha=0.35, rstride=3, cstride=4)
ax.plot_surface(*ell(ad, S=S), color=kinds[1][2], alpha=0.6, linewidth=0)
for sgn in (1, -1):
    p0 = sgn*6.0*u; d = -sgn*1.8*u
    ax.quiver(*p0, *d, color="#e09f3e", linewidth=3, arrow_length_ratio=0.35)
setup(ax, "Far side presses on the down needle")
ax.text2D(0.0, -0.06, "Press along (−1, −1, 2), about 12%, plus a timing (phase).\nGray mesh = before, blue = after. Shown on the log-drawn shape.", transform=ax.transAxes, fontsize=10)
# Panel D: honest key
ax = fig.add_subplot(2, 3, 6); ax.axis("off")
key = ("HOW TO READ THIS PICTURE\n\n"
       "To scale (by the log rule)\n"
       " • The order and relative lengths of each needle's three axes:\n"
       "   drawn length = 1 + log10(mass ÷ lightest mass).\n"
       " • The direction of the far-side press.\n\n"
       "Exaggerated so it can be seen\n"
       " • The real needles are far thinner: the up quarks' are\n"
       "   130,800 times longer than wide. On a true scale they\n"
       "   would be invisible lines.\n"
       " • The small turns (2.4° ×5, 0.21° ×20).\n\n"
       "Cannot be drawn\n"
       " • The axes are families (first, second, third),\n"
       "   not directions in ordinary space.\n"
       " • The numbers are complex, so every axis carries a hidden\n"
       "   twist (its area is the six equal triangles, J/2).\n"
       " • Signs of masses, and the far side itself (a fourth,\n"
       "   very heavy direction) do not appear.\n\n"
       "Established: matrices as ellipsoids (Cauchy, Sylvester);\n"
       "turns (Cabibbo; Kobayashi & Maskawa). TRPG: neutrino\n"
       "masses, far-side press.")
ax.text(0.0, 1.0, key, va="top", ha="left", fontsize=10.5, family="DejaVu Sans",
        bbox=dict(boxstyle="round,pad=0.6", fc="#fffbe6", ec="#F5C518", lw=2))
plt.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.04, wspace=0.05, hspace=0.12)
out = "TRPG_Needle_Geometry.png"
plt.savefig(out, dpi=170); print("saved", out)
