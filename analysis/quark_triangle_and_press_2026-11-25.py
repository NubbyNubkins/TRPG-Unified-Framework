"""Independent recheck of the quark-mixing triangle (lambda = 1/sqrt20, A = sqrt(2/3), alpha = 90, beta = 22.5),
the 3x3 down-quark 'press' matrix, the lepton twist, the generation-mass shares and the dark-energy formula.
Measured values: PDG 2026 where the repository already uses it (V_us, gamma, alpha, first row, V_ud),
PDG 2024 otherwise (V_cb, V_ub, sin 2beta, |V_td/V_ts|). One table, used for every test."""
import math, cmath
import numpy as np
from scipy.stats import chi2 as chi2dist

MEAS = {  # value, error (asymmetric errors: the side facing the prediction is used)
    "Vus": (0.22431, 0.00085),          # PDG 2026 (repo), error PDG 2024
    "Vcb": (0.0411, 0.0012),            # PDG 2024
    "Vub": (0.00382, 0.00020),          # PDG 2024
    "sin2b": (0.709, 0.011),            # PDG 2024
    "alpha": (84.1, 3.7),               # PDG 2026 (repo), +3.7 side
    "gamma": (66.4, 2.8),               # PDG 2026 (repo)
    "Vtd/Vts": (0.207, 0.0032),         # PDG 2024 (0.001 stat + 0.003 lattice combined)
}
VUD = (0.97367, 0.00032)                # superallowed beta decays, PDG 2024

def angles(V):
    z1 = V[0, 0] * np.conj(V[0, 2]); z2 = V[1, 0] * np.conj(V[1, 2]); z3 = V[2, 0] * np.conj(V[2, 2])
    beta = math.degrees(cmath.phase(-z2 / z3)); gamma = math.degrees(cmath.phase(-z1 / z2))
    return abs(beta), abs(gamma), 180 - abs(beta) - abs(gamma)

def observables(V):
    b, g, a = angles(V)
    return {"Vus": abs(V[0, 1]), "Vcb": abs(V[1, 2]), "Vub": abs(V[0, 2]), "sin2b": math.sin(math.radians(2 * b)),
            "alpha": a, "gamma": g, "Vtd/Vts": abs(V[2, 0]) / abs(V[2, 1])}

def report(label, V):
    o = observables(V); c2 = 0
    print(f"  {label}")
    for k, (m, e) in MEAS.items():
        p = (o[k] - m) / e; c2 += p * p
        print(f"     {k:8s} pred {o[k]:9.5f}  meas {m:9.5f} +- {e:.5f}  pull {p:+5.2f}")
    print(f"     chi2 = {c2:.2f} for 7 observables (p = {1 - chi2dist.cdf(c2, 7):.3f} if nothing had been chosen)")
    vud = abs(V[0, 0]); print(f"     |V_ud| = {vud:.5f} vs {VUD[0]} +- {VUD[1]}: {(vud - VUD[0]) / VUD[1]:+.2f} sigma;"
                              f" first row sum {abs(V[0,0])**2+abs(V[0,1])**2+abs(V[0,2])**2:.6f} vs 0.9984 +- 0.0007 (PDG 2026)")
    return c2

def standard(s12, s23, s13, d):
    c12, c23, c13 = (math.sqrt(1 - x * x) for x in (s12, s23, s13)); e = cmath.exp(1j * d)
    return np.array([[c12 * c13, s12 * c13, s13 / e],
                     [-s12 * c23 - c12 * s23 * s13 * e, c12 * c23 - s12 * s23 * s13 * e, s23 * c13],
                     [s12 * s23 - c12 * c23 * s13 * e, -c12 * s23 - s12 * c23 * s13 * e, c23 * c13]])

print("1. RIGHT UNITARITY TRIANGLE: lambda = 1/sqrt(20), A = sqrt(2/3), beta = 22.5, gamma = 67.5 (alpha = 90)")
lam, A = 1 / math.sqrt(20), math.sqrt(2 / 3)
s12, s23 = lam, A * lam ** 2
Ru = math.sin(math.radians(22.5))                    # right angle at alpha: R_u = sin beta
zb = Ru * cmath.exp(1j * math.radians(67.5))         # rhobar + i etabar
s13e = zb * s12 * s23 * math.sqrt(1 - s23 ** 2) / (math.sqrt(1 - s12 ** 2) * (1 - zb * s23 ** 2))  # exact (Charles et al. 2005)
V7 = standard(s12, s23, abs(s13e), cmath.phase(s13e))
print(f"     unitarity error {np.max(np.abs(V7.conj().T @ V7 - np.eye(3))):.1e}; J = {np.imag(V7[0,1]*V7[1,2]*np.conj(V7[0,2])*np.conj(V7[1,1])):.3e}")
report("triangle with PDG 2026/2024 direct measurements", V7)

print("\n2. 3x3 DOWN-QUARK MATRIX WITH SET WEIGHTS (M11 = 2 r_d, M22 = r_s + 4 r_s^2, |M12| = sqrt(r_d r_s),"
      " |M23| = 2 r_s, M13 = sqrt(2/3) sqrt(r_d r_s) e^{-i 67.5})")
for lab, rd, rs in (("at 162.5 GeV (r_d = 0.001, r_s = 0.020)", 0.00275 / 2.75, 0.055 / 2.75),
                    ("at M_Z     (r_d, r_s from 2.67, 53.16, 2855 MeV)", 2.67 / 2855, 53.16 / 2855)):
    m12, m23 = math.sqrt(rd * rs), 2 * rs
    m13 = math.sqrt(2 / 3) * m12 * cmath.exp(-1j * math.radians(67.5))
    M = np.array([[2 * rd, m12, m13], [m12, rs + m23 ** 2, m23], [np.conj(m13), m23, 1]], dtype=complex)
    w, U = np.linalg.eigh(M); i = np.argsort(np.abs(w)); w, U = w[i], U[:, i]
    # rephase to the standard convention (V_ud, V_us, V_cb, V_tb real)
    U = U @ np.diag([np.exp(-1j * np.angle(U[0, 0])), 1, 1])
    U[:, 1] *= np.exp(-1j * np.angle(U[0, 1])); U[:, 2] *= np.exp(-1j * np.angle(U[1, 2]))
    print(f"     eigenvalues/m_b: {w[0]:.6f} (target {rd:.6f}), {w[1]:.5f} (target {rs:.5f})")
    report(lab, U)

print("\n3. LEPTON TWIST (Brannen 2006 form sqrt(m_k) = a (1 + sqrt2 cos(delta + 2 pi k / 3)))")
me, mmu, mtau = 0.51099895069, 105.6583755, (1776.93, 0.09)
def masses(a, d): return [(a * (1 + math.sqrt(2) * math.cos(d + 2 * math.pi * k / 3))) ** 2 for k in range(3)]
d = 2 / 9
ms = masses(1, d); order = sorted(range(3), key=lambda k: ms[k])
a = math.sqrt(mmu / ms[order[1]])
me_p, mmu_p, mtau_p = (masses(a, d)[k] for k in order)
print(f"     delta = 2/9 anchored on m_mu: m_e = {me_p:.8f} MeV ({(me_p/me-1)*1e6:+.1f} ppm), m_tau = {mtau_p:.4f} MeV"
      f" ({(mtau_p-mtau[0])/mtau[0]*100:+.4f}%, {(mtau_p-mtau[0])/mtau[0]*0+ (mtau_p-mtau[0])/mtau[1]:+.1f} sigma)")
Q = (me + mmu + mtau[0]) / (math.sqrt(me) + math.sqrt(mmu) + math.sqrt(mtau[0])) ** 2
print(f"     Koide Q from measured masses = {Q:.7f} (2/3 = 0.6666667); this 45-degree rule is Koide (1982)")

print("\n4. GEN-3 SHARE OF EACH TRACE (running masses at 162.5 GeV)")
for lab, m in (("up", (0.00127, 0.62, 162.5)), ("down", (0.00275, 0.055, 2.75)), ("lepton", (me, mmu, mtau[0]))):
    print(f"     {lab:6s}: heaviest/trace = {m[2]/sum(m)*100:.2f}%")

print("\n5. DARK ENERGY FORMULA rho = B^(3/2) / [3 (16 pi^2)^2 M_Pl,red^2]")
hc3 = 0.1973269804 ** 3  # GeV^3 fm^3
MPlr, H0, OmL = 2.435e18, 1.437642e-42, 0.6847
obs = (3 * MPlr ** 2 * H0 ** 2 * OmL) ** 0.25
for B in (57, 60, 64):
    Bg = B * 1e-3 * hc3   # MeV/fm^3 -> GeV^4
    print(f"     B = {B}: rho^(1/4) = {(Bg**1.5/(3*(16*math.pi**2)**2*MPlr**2))**0.25*1e12:.3f} meV  (observed {obs*1e12:.3f} meV)")
print("     This was the closest of 10 coefficient/power forms tried; with 10 tries")
print("     spread over ~3 decades, landing within 5% by chance is expected about once.")
