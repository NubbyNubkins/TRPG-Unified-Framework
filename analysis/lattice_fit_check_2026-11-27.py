"""Does TRPG's main-line particle content fit the lattice construction of
N. Craig, "A gauge-invariant measure for lattice chiral fermions",
arXiv:2610.09093 (2026)?  Checks only what the construction asks of the
matter content.  All fields are written as left-handed Weyl fermions.
"""
from fractions import Fraction as F

# (name, SU3 dim, SU2 dim, hypercharge Y, SU(5) home)
SM_GEN = [
    ("Q",   3, 2, F(1, 6),  "10"), ("u^c", 3, 1, F(-2, 3), "10"),
    ("e^c", 1, 1, F(1),     "10"), ("d^c", 3, 1, F(1, 3),  "5bar"),
    ("L",   1, 2, F(-1, 2), "5bar"),
]
N_PARTNER = [("N^c", 1, 1, F(0), "1")]                     # hidden neutrino partner
FAR_SIDE = [                                               # E6 27 minus the 16
    ("D",    3, 1, F(-1, 3), "5"), ("Lfar",    1, 2, F(1, 2),  "5"),
    ("Dbar", 3, 1, F(1, 3),  "5bar"), ("Lfarbar", 1, 2, F(-1, 2), "5bar"),
    ("S", 1, 1, F(0), "1"),                                # the extra neutral state
]
T = {1: F(0), 2: F(1, 2), 3: F(1, 2)}                      # Dynkin index of fundamentals

def anomaly_sums(gen):
    # SU(3)^3: +1 per triplet, -1 per antitriplet (u^c, d^c, Dbar are antitriplets)
    anti3 = {"u^c", "d^c", "Dbar"}
    su3 = sum((-1 if f[0] in anti3 else 1) * f[2] for f in gen if f[1] == 3)
    su3sq_y = sum(T[3] * f[2] * f[3] for f in gen if f[1] == 3)
    su2sq_y = sum(T[2] * f[1] * f[3] for f in gen if f[2] == 2)
    y3 = sum(f[1] * f[2] * f[3] ** 3 for f in gen)
    grav = sum(f[1] * f[2] * f[3] for f in gen)
    doublets = sum(f[1] for f in gen if f[2] == 2)          # Witten SU(2) count
    return dict(SU3cubed=su3, SU3sq_Y=su3sq_y, SU2sq_Y=su2sq_y, Y_cubed=y3,
                grav_Y=grav, SU2_doublets=doublets)

def bl_anomalies(gen):
    # B-L charges (left-handed): Q 1/3, u^c,d^c -1/3, L -1, e^c,N^c +1,
    # vector-like far side: D 1/3, Dbar -1/3, Lfar +1, Lfarbar -1, S 0
    q = {"Q": F(1, 3), "u^c": F(-1, 3), "d^c": F(-1, 3), "L": F(-1), "e^c": F(1),
         "N^c": F(1), "D": F(1, 3), "Dbar": F(-1, 3), "Lfar": F(1), "Lfarbar": F(-1), "S": F(0)}
    bl3 = sum(f[1] * f[2] * q[f[0]] ** 3 for f in gen)
    blg = sum(f[1] * f[2] * q[f[0]] for f in gen)
    su2sq = sum(T[2] * f[1] * q[f[0]] for f in gen if f[2] == 2)
    su3sq = sum(T[3] * f[2] * q[f[0]] for f in gen if f[1] == 3)
    y_bl2 = sum(f[1] * f[2] * f[3] * q[f[0]] ** 2 for f in gen)
    bl_y2 = sum(f[1] * f[2] * q[f[0]] * f[3] ** 2 for f in gen)
    return dict(BL_cubed=bl3, BL_grav=blg, SU2sq_BL=su2sq, SU3sq_BL=su3sq,
                Y_BLsq=y_bl2, BL_Ysq=bl_y2)

# Craig's global condition: under the SU(2) acting on two fundamental colours of
# SU(5), every half-integer isospin must appear an even number of times.
SU5_DOUBLETS = {"10": 3, "5": 1, "5bar": 1, "1": 0}
SU5_CUBIC = {"10": 1, "5": 1, "5bar": -1, "1": 0}

def su5_content(gen):
    # one SU(5) irrep per group (Q,u^c,e^c)->10 ; (d^c,L)->5bar ; (D,Lfar)->5 ; (Dbar,Lfarbar)->5bar
    irreps = []
    if {"Q", "u^c", "e^c"} <= {f[0] for f in gen}: irreps.append("10")
    if {"d^c", "L"} <= {f[0] for f in gen}: irreps.append("5bar")
    if {"D", "Lfar"} <= {f[0] for f in gen}: irreps.append("5")
    if {"Dbar", "Lfarbar"} <= {f[0] for f in gen}: irreps.append("5bar")
    irreps += ["1"] * sum(1 for f in gen if f[4] == "1")
    return irreps

def beta_coeffs(n_gen, far_side, n_higgs=1):
    # one-loop b (positive = asymptotically free), b = 11/3 C2 - 2/3 sum_Weyl T - 1/3 sum_cplx-scalar T
    weyl3 = 4 * n_gen + (2 * n_gen if far_side else 0)       # Weyl triplets
    weyl2 = 4 * n_gen + (2 * n_gen if far_side else 0)       # Weyl doublets
    b3 = F(11) - F(2, 3) * weyl3 * F(1, 2)
    b2 = F(22, 3) - F(2, 3) * weyl2 * F(1, 2) - F(1, 3) * n_higgs * F(1, 2)
    return b3, b2

out = []
p = out.append
for label, gen in [("SM generation + hidden partner (16)", SM_GEN + N_PARTNER),
                   ("TRPG main-line generation (27)", SM_GEN + N_PARTNER + FAR_SIDE)]:
    p(f"== {label}: {sum(f[1]*f[2] for f in gen)} Weyl states")
    for k, v in anomaly_sums(gen).items():
        p(f"   {k:14s} {v}")
    for k, v in bl_anomalies(gen).items():
        p(f"   {k:14s} {v}")
    irreps = su5_content(gen)
    p(f"   SU(5) homes      {' + '.join(irreps)}")
    p(f"   SU(5) cubic sum  {sum(SU5_CUBIC[r] for r in irreps)}")
    nd = sum(SU5_DOUBLETS[r] for r in irreps)
    p(f"   doublets under SU(2) on two SU(5) colours: {nd} -> {'even, passes' if nd % 2 == 0 else 'ODD, fails'}")
p("")
p("One-loop coefficients (positive = asymptotically free), three generations:")
for lab, fs in [("SM + 3 hidden partners", False), ("TRPG with far side", True)]:
    b3, b2 = beta_coeffs(3, fs)
    p(f"   {lab:24s} b3 = {b3}, b2 = {b2}")
p("   (U(1)_Y is not asymptotically free in either; same status as the SM.)")
print("\n".join(out))
