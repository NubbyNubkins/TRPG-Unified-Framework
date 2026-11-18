"""Check (18 Nov 2026): heavy partner masses from the FULL 3x3 seesaw with the literal shared-matrix rule.
Rule: Dirac link = D = diag(m_e, m_mu, m_tau) in the charged-lepton basis. m_nu = U diag(m) U^T (U = PMNS, TM1).
M_R = -D m_nu^-1 D ; heavy masses = Takagi (singular) values of M_R.
Earlier versions (up to v1.2.2) quoted M_k = m_l,k^2 / m_nu,k, which is exact only for U = 1 (no mixing).
Also: m_bb for the natural reading (Majorana phases zero) across the TM1 delta range.
Credit: Minkowski 1977; Yanagida 1979; Gell-Mann, Ramond & Slansky 1979; Casas & Ibarra 2001; Xing & Zhou 2007; Albright & Rodejohann 2009; NuFIT 6.0."""
import math, numpy as np
D = np.diag([0.51099895e-3, 0.1056583755, 1.77693])
MNU = np.array([0.00036, 0.0087, 0.0501])*1e-9
s13 = 0.02195
def pmns(s12, s13, s23, d, a21=0, a31=0):
    c12, c13, c23 = math.sqrt(1-s12), math.sqrt(1-s13), math.sqrt(1-s23); s12, s13, s23 = math.sqrt(s12), math.sqrt(s13), math.sqrt(s23)
    e = np.exp(1j*d)
    U = np.array([[c12*c13, s12*c13, s13/e], [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, s23*c13],
                  [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13]])
    return U @ np.diag([1, np.exp(1j*a21/2), np.exp(1j*a31/2)])
def tm1(s23):
    s12 = 1 - 2/(3*(1-s13))           # TM1: |U_e1|^2 = 2/3
    t13 = math.sqrt(s13); c13 = math.sqrt(1-s13); t23 = math.asin(math.sqrt(s23))
    cosd = -(1-5*s13)*(1/math.tan(2*t23))/(2*math.sqrt(2)*t13*math.sqrt(1-3*s13))
    cosd = max(-1, min(1, cosd)); return s12, math.acos(cosd)
def heavy(U):
    mnu = U @ np.diag(MNU) @ U.T; MR = -D @ np.linalg.inv(mnu) @ D
    return np.sort(np.linalg.svd(MR, compute_uv=False))
print("diagonal (earlier) formula m_l^2/m_nu:", ", ".join(f"{x:.2e}" for x in np.diag(D)**2/MNU), "GeV")
print("\nFull matrix, TM1, Majorana phases 0, delta in the upper half (~270) and lower half:")
for s23 in [0.47, 0.50, 0.52, 0.55, 0.57]:
    s12, d = tm1(s23)
    for dd in (2*math.pi-d, d):
        M = heavy(pmns(s12, s13, s23, dd))
        print(f"  s23={s23:.2f} s12={s12:.3f} delta={math.degrees(dd):5.1f}: heavy = {M[0]:.2e}, {M[1]:.2e}, {M[2]:.2e} GeV")
rng = np.random.default_rng(18); lo = np.full(3, 1e99); hi = np.zeros(3)
for _ in range(20000):
    s23 = rng.uniform(0.435, 0.585); s12, d = tm1(s23); dd = 2*math.pi - d
    M = heavy(pmns(s12, s13, s23, dd, *rng.uniform(0, 2*math.pi, 2))); lo = np.minimum(lo, M); hi = np.maximum(hi, M)
print("\nRange over theta23 (3 sigma) and all Majorana phases, TM1 delta near 270:")
for k in range(3): print(f"  partner {k+1}: {lo[k]:.2e} to {hi[k]:.2e} GeV")
print("\nm_bb, natural reading (Majorana phases zero), TM1:")
for s23 in [0.45, 0.50, 0.55, 0.58]:
    s12, d = tm1(s23); dd = 2*math.pi - d; U = pmns(s12, s13, s23, dd)
    print(f"  s23={s23:.2f} delta={math.degrees(dd):.0f}: m_bb = {abs((U[0]**2*MNU).sum())*1e12:.2f} meV")
vals = []
for _ in range(20000):
    s23 = rng.uniform(0.435, 0.585); s12, d = tm1(s23); U = pmns(s12, s13, s23, 2*math.pi-d, *rng.uniform(0, 2*math.pi, 2)); vals.append(abs((U[0]**2*MNU).sum())*1e12)
print(f"  all phases free: {min(vals):.2f} to {max(vals):.2f} meV")
