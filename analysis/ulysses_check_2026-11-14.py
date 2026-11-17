"""Independent check of the an earlier test-31 matter-excess cases with the public leptogenesis code ULYSSES (14 Nov 2026).
ULYSSES: Granelli, Moffat, Perez-Gonzalez, Schulz & Turner, Comput. Phys. Commun. 262 (2021) 107813; v2 (2023). Version used: 3.0.2.
Solver: EtaB_3DME = density-matrix equations for three decaying heavy neutrinos (Blanchet, Di Bari, Jones & Marzola 2013).

Install (the plain 'pip install ulysses' fails on its old 'progressbar' helper; this works):
    pip install numpy scipy numba odeintw matplotlib progressbar2
    pip install --no-deps ulysses

Convention mapping. Ours: h = sqrt(M) R sqrt(m) U^dagger / v (rows = partners). ULYSSES: h = (i/v) U sqrt(m) R^T sqrt(M) (rows = flavours).
These are the same matrix (transposed, times i) when ULYSSES is given the SAME R angles and the three PMNS phases with the sign flipped.
ULYSSES integrates in z = M1/T, so z must start below M1/M3 (about 1e-5 here); zmin = 1e-7 is used (default 1e-3 would miss partner 3)."""
import json, math, warnings
warnings.filterwarnings("ignore")
import numpy as np, ulysses
cases = json.load(open("best_leptogenesis_settings_2026-11-14.json"))
deg = 180/math.pi; OBS = 6.1e-10
def pdict(c):
    m1, m2, m3 = c["m_light_eV"]; M = c["M_heavy_GeV"]
    return dict(m=math.log10(m1), M1=math.log10(M[0]), M2=math.log10(M[1]), M3=math.log10(M[2]),
                delta=-c["delta_rad"]*deg, a21=-c["a21_rad"]*deg, a31=-c["a31_rad"]*deg,
                x1=c["x_rad"][0]*deg, x2=c["x_rad"][1]*deg, x3=c["x_rad"][2]*deg,
                y1=c["y_rad"][0]*deg, y2=c["y_rad"][1]*deg, y3=c["y_rad"][2]*deg,
                t12=math.degrees(math.asin(math.sqrt(c["sin2_12"]))), t13=math.degrees(math.asin(math.sqrt(c["sin2_13"]))),
                t23=math.degrees(math.asin(math.sqrt(c["sin2_23"]))), m2solar=(m2**2-m1**2)*1e-18, m2atm=(m3**2-m1**2)*1e-18)
def run(model, c, init, zmin=1e-7, zsteps=3000):
    L = getattr(ulysses, model)(zmin=zmin, zmax=1000, zsteps=zsteps, ordering=0, initial_abundance=init)
    return L(pdict(c)), L
print("ULYSSES version 3.0.2; observed eta_B = 6.1e-10. |eta_B| is quoted (the sign flips if every phase is mirrored).\n")
print(f"{'case':>4} {'ours':>18} {'ULYSSES 3DME':>13} {'with scatterings':>17} {'zero start':>11} {'ULYSSES/observed':>17}")
for c in cases:
    e, L = run("EtaB_3DME", c, 1); es, _ = run("EtaB_3DMEsct", c, 1); e0, _ = run("EtaB_3DME", c, 0)
    print(f"{c['case']:>4} {c['our_eta_B_full_equations']:>18.2e} {abs(e):>13.2e} {abs(es):>17.2e} {abs(e0):>11.2e} {abs(e)/OBS:>17.2f}")
c = cases[0]
print("\nNumerical stability (case 1):", ", ".join(f"zmin {zm:g}/steps {zs}: {abs(run('EtaB_3DME', c, 1, zm, zs)[0]):.4e}" for zm, zs in [(1e-7, 3000), (1e-8, 3000), (1e-7, 10000), (1e-6, 1000)]))
print("\nRule 5 check: sizes of the neutrino's Dirac links in each case (should be the lepton masses 0.000511, 0.10566, 1.77693 GeV)")
for c in cases:
    e, L = run("EtaB_3DME", c, 1); sv = np.sort(np.linalg.svd(L.h*L.v, compute_uv=False))
    print(f"  case {c['case']}: {sv[0]:.5f}, {sv[1]:.5f}, {sv[2]:.5f} GeV   (first two are off by x{sv[0]/0.000511:.1f} and x{sv[1]/0.10566:.2f})")
print("  -> the earlier 'best' cases used a quarter-turn (x3 near 90 degrees) that swaps which neutrino pairs with which partner;")
print("     that changes the first two links, so those cases did not hold rule 5. See matter_excess_rule5_exact_2026-11-14.py.")
