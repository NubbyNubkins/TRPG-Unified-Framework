"""Independent TOV + tidal-deformability solver for colour-flavour-locked quark stars (geometric units, km).
Equation of state: Alford, Braby, Paris and Reddy (2005), Omega = -3 a4 mu^4/(4 pi^2) + 3 (m_s^2 - 4 Delta^2) mu^2/(4 pi^2) + B.
Tidal Love number: Hinderer (2008) with the surface-density correction of Postnikov, Prakash and Lattimer (2010)."""
import math, sys
import numpy as np
from scipy.integrate import solve_ivp

HC = 197.3269804                 # MeV fm
K = 1.3234e-6                    # 1 MeV/fm^3 in km^-2 (G/c^4)
MSUN = 1.476625                  # km

class EOS:
    def __init__(s, B, ms, D, a4):
        s.B4 = B * HC ** 3; s.c4 = 3 * a4 / (4 * math.pi ** 2); s.c2 = 3 * (ms ** 2 - 4 * D ** 2) / (4 * math.pi ** 2); s.B = B
    def x(s, P):                 # mu^2 in MeV^2, P in MeV/fm^3
        return (s.c2 + math.sqrt(s.c2 ** 2 + 4 * s.c4 * (P * HC ** 3 + s.B4))) / (2 * s.c4)
    def eps(s, P):               # MeV/fm^3
        return 3 * P + 4 * s.B + 2 * s.c2 * s.x(P) / HC ** 3
    def deps(s, P):
        return 3 + 2 * s.c2 / math.sqrt(s.c2 ** 2 + 4 * s.c4 * (P * HC ** 3 + s.B4))
    def EA(s):
        return 3 * math.sqrt(s.x(0.0))

def star(eos, Pc):
    def rhs(r, u):
        m, P, y = u
        P = max(P, 0.0)
        e, p = eos.eps(P) * K, P * K
        f = 1 - 2 * m / r
        dm = 4 * math.pi * r * r * e
        dP = -(e + p) * (m + 4 * math.pi * r ** 3 * p) / (r * r * f)
        nu = 2 * (m + 4 * math.pi * r ** 3 * p) / (r * r * f)
        F = (1 - 4 * math.pi * r * r * (e - p)) / f
        Q = 4 * math.pi * (5 * e + 9 * p + (e + p) * eos.deps(P)) / f - 6 / (r * r * f) - nu * nu
        dy = -(y * y + y * F + r * r * Q) / r
        return [dm, dP / K, dy]
    r0 = 1e-4
    m0 = 4 / 3 * math.pi * r0 ** 3 * eos.eps(Pc) * K
    ev = lambda r, u: u[1] - 1e-10; ev.terminal = True; ev.direction = -1
    sol = solve_ivp(rhs, [r0, 40], [m0, Pc, 2.0], events=ev, rtol=1e-9, atol=1e-12, method="LSODA")
    R = sol.t_events[0][0]; M, _, yR = sol.y_events[0][0]
    yR -= 4 * math.pi * R ** 3 * eos.eps(0.0) * K / M      # Postnikov surface jump
    C = M / R
    k2 = (8 * C ** 5 / 5) * (1 - 2 * C) ** 2 * (2 + 2 * C * (yR - 1) - yR) / (
        2 * C * (6 - 3 * yR + 3 * C * (5 * yR - 8)) + 4 * C ** 3 * (13 - 11 * yR + C * (3 * yR - 2) + 2 * C ** 2 * (1 + yR))
        + 3 * (1 - 2 * C) ** 2 * (2 - yR + 2 * C * (yR - 1)) * math.log(1 - 2 * C))
    return M / MSUN, R, 2 / 3 * k2 / C ** 5

def branch(eos):
    rows = [star(eos, Pc) for Pc in np.geomspace(5, 3000, 70)]
    M = np.array([r[0] for r in rows]); R = np.array([r[1] for r in rows]); L = np.array([r[2] for r in rows])
    i = int(np.argmax(M)); M, R, L = M[: i + 1], R[: i + 1], L[: i + 1]
    at = lambda arr, m: float(np.interp(m, M, arr)) if M[-1] >= m else float("nan")
    return M[-1], at(R, 2.08), at(L, 1.4), at(L, 2.0)

def ud_EA(B, a4, D=0.0):
    """Two-flavour neutral quark matter at zero pressure; D > 0 adds 2SC pairing, +mubar^2 D^2/pi^2 to the pressure
    (Alford and Rajagopal 2002), with mu_d = 2^(1/3) mu_u kept from neutrality (an approximation)."""
    B4 = B * HC ** 3; r = 2 ** (1 / 3)
    lo, hi = 50.0, 800.0
    for _ in range(100):
        mu = (lo + hi) / 2
        P = a4 * (mu ** 4 + (r * mu) ** 4) / (4 * math.pi ** 2) + ((1 + r) * mu / 2) ** 2 * D ** 2 / math.pi ** 2 - B4
        lo, hi = (lo, mu) if P > 0 else (mu, hi)
    return mu * (1 + 2 * r)

if __name__ == "__main__":
    print("CFL QUARK STARS, B = 60 MeV/fm^3 unless stated (independent solver)")
    print(" a4    m_s  Delta |  M_max   R_2.08   Lam_1.4  Lam_2.0   E/A(CFL)")
    cases = [(a4, 95, D) for a4 in (0.70, 0.85, 1.00) for D in (0, 70, 90, 105, 125)] + [(1.0, 100, 82), (1.0, 100, 97)]
    for a4, ms, D in cases:
        e = EOS(60, ms, D, a4)
        Mx, R208, L14, L20 = branch(e)
        print(f" {a4:.2f}  {ms:3d}  {D:4d}  | {Mx:6.3f}  {R208:7.2f}  {L14:8.0f}  {L20:7.0f}   {e.EA():7.1f}", flush=True)
    print("\nTWO-FLAVOUR (ud) MATTER E/A at zero pressure (must stay above 930 MeV for ordinary nuclei to be stable)")
    for B in (57, 60, 64):
        for a4 in (0.70, 1.00):
            print(f"  B={B} a4={a4:.2f}: unpaired {ud_EA(B, a4):6.1f} MeV | 2SC-paired, gap 90: {ud_EA(B, a4, 90):6.1f} | gap 125: {ud_EA(B, a4, 125):6.1f}")
