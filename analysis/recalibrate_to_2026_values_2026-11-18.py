"""Recalibration of the headline numbers to the newest measured values (18 November 2026).

Every number below is recomputed from one table of inputs, so that no headline rests on an older value.
Inputs: PDG 2026 (Review of Particle Physics, summary tables and CKM review of 1 June 2026),
JUNO 2026 (arXiv:2609.37895, 207.2 live days), NuFIT 6.0 for the global fit, and the PDG 2026
table of physical constants for alpha_s. Nothing is fitted here; each line says what changed.
The Higgs-mass rule is rerun with the functions of higgs_at_the_grain_2026-11-17.py.
"""
import math, os, sys, io, contextlib, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- inputs
PDG24 = dict(me=0.51099895, mmu=105.6583755, mtau=(1776.93, 0.09), ms=(93.5, 0.8), md=(4.70, 0.07),
             mu=(2.16, 0.07), mt=(172.57, 0.29), mH=(125.20, 0.11), mW=80.3692, Vus=0.2243)
PDG26 = dict(me=0.51099895069, mmu=105.6583755, mtau=(1776.93, 0.09), ms=(92.9, 0.7), md=(4.70, 0.07),
             mu=(2.16, 0.07), mt=(172.60, 0.27), mH=(125.13, 0.11), mW=80.3625, Vus=0.22431,
             mu_md=(0.462, 0.020), first_row=(0.9984, 0.0007), gamma=(66.4, 2.8), alpha=(84.1, 3.7, 3.0),
             tau_tau=(290.3, 0.5), s13sq=0.0217, alpha_s=(0.1180, 0.0009))
JUNO26 = dict(dm21=(7.388e-5, 0.078e-5), dm31=2.509e-3, s12sq=(0.3036, 0.0064))
NUFIT6 = dict(dm21=7.49e-5, dm31_noSK=2.513e-3, dm31_SK=2.534e-3, s13sq=0.02195)
LATTICE_AS = (0.11873, 0.00056)   # Dalla Brida and others 2025
Q_ETA = (22.1, 0.7)               # Colangelo, Lanz, Leutwyler, Passemar 2018 (PDG 2026 quark-mass review)

def z(pred, meas, err): return (pred - meas) / err
def line(s): print(s, flush=True)

line("RECALIBRATION OF THE HEADLINE NUMBERS TO PDG 2026 AND JUNO 2026")
line("=" * 100)

# 1. tau mass from the 45 degree rule
a, b = math.sqrt(PDG26["me"]), math.sqrt(PDG26["mmu"])
A, B, C = 1.0, -4 * (a + b), 3 * (PDG26["me"] + PDG26["mmu"]) - 2 * (a + b) ** 2
mtau = ((-B + math.sqrt(B * B - 4 * A * C)) / 2) ** 2
line(f"1. tau mass (45 deg rule): {mtau:.3f} MeV vs PDG 2026 {PDG26['mtau'][0]} +- {PDG26['mtau'][1]}: "
     f"{z(mtau, *PDG26['mtau']):+.1f} sigma   [unchanged from v1.3.2]")

# 2. down quark from the Cabibbo relation sin(theta_C) = sqrt(m_d/m_s)
ratio = 1 / PDG26["Vus"] ** 2
md_pred = PDG26["ms"][0] / ratio
rat_meas = PDG26["ms"][0] / PDG26["md"][0]
rat_err = rat_meas * math.hypot(PDG26["ms"][1] / PDG26["ms"][0], PDG26["md"][1] / PDG26["md"][0])
line(f"2. m_s/m_d from the Cabibbo angle: {ratio:.2f}. PDG 2026 masses give {rat_meas:.2f} +- {rat_err:.2f}: "
     f"{z(ratio, rat_meas, rat_err):+.1f} sigma; m_d = {md_pred:.2f} MeV vs {PDG26['md'][0]} +- {PDG26['md'][1]}: "
     f"{z(md_pred, *PDG26['md']):+.1f} sigma")
line("   v1.3.2 compared with '20.4 +- 0.3' (-1.7 sigma). That number has no source in PDG 2024 or 2026, and the PDG masses the")
line("   same script used (93.5 / 4.70 = 19.9) disagree with it. The check passes either way; the comparison value is corrected.")
line("   The sine form was chosen over the tangent form (18.9, -2.7 sigma here) after both were tried: a fitted choice of form.")

# 3. up quark through eta -> 3 pi (Leutwyler's ellipse)
r_ud = math.sqrt(1 - (ratio / Q_ETA[0]) ** 2)
mu_pred = r_ud * md_pred
lo = math.sqrt(1 - (ratio / (Q_ETA[0] - Q_ETA[1])) ** 2) * md_pred
hi = math.sqrt(1 - (ratio / (Q_ETA[0] + Q_ETA[1])) ** 2) * md_pred
line(f"3. up quark: m_u/m_d = {r_ud:.3f} (PDG 2026 0.462 +- 0.020: {z(r_ud, *PDG26['mu_md']):+.1f} sigma); "
     f"m_u = {mu_pred:.2f} MeV (range {lo:.2f}-{hi:.2f} from Q) vs 2.16 +- 0.07: {z(mu_pred, *PDG26['mu']):+.1f} sigma")

# 4. neutrino masses from the 45 degree rule with roots (-,+,+)
def nu_masses(d21, d31):
    def Q(m1):
        m2, m3 = math.sqrt(m1 * m1 + d21), math.sqrt(m1 * m1 + d31)
        return (m1 + m2 + m3) / (-math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)) ** 2 - 2 / 3
    lo, hi = 1e-7, 2e-3
    for _ in range(200):
        mid = (lo + hi) / 2
        if Q(lo) * Q(mid) <= 0: hi = mid
        else: lo = mid
    m1 = (lo + hi) / 2
    return m1, math.sqrt(m1 * m1 + d21), math.sqrt(m1 * m1 + d31)
line("4. neutrino masses (45 deg rule, roots -,+,+):")
for lab, d21, d31 in (("JUNO 2026", JUNO26["dm21"][0], JUNO26["dm31"]),
                      ("NuFIT 6.0 without SK", NUFIT6["dm21"], NUFIT6["dm31_noSK"]),
                      ("NuFIT 6.0 with SK", NUFIT6["dm21"], NUFIT6["dm31_SK"])):
    m = nu_masses(d21, d31)
    line(f"   {lab:22s}: m1 {m[0]*1e3:.3f}, m2 {m[1]*1e3:.2f}, m3 {m[2]*1e3:.2f} meV, sum {sum(m)*1e3:.2f} meV")
mJ = nu_masses(JUNO26["dm21"][0], JUNO26["dm31"])

# 5. the small-splitting hint m2/m3 = 3 - 2 sqrt2
r = 3 - 2 * math.sqrt(2)
def hint(d31, m1=0.36e-3): return r * r * (m1 * m1 + d31) - m1 * m1
line("5. small splitting from m2/m3 = 3 - 2 sqrt2 (a lead):")
for lab, d31 in (("NuFIT 6.0-era 2.519e-3 (v1.2.2)", 2.519e-3), ("JUNO 2026 2.509e-3 (v1.3.2)", JUNO26["dm31"])):
    p = hint(d31)
    line(f"   {lab:32s}: {p*1e5:.3f}e-5 eV^2 vs JUNO 7.388 +- 0.078: {z(p, *JUNO26['dm21']):+.2f} sigma")
p0 = hint(JUNO26["dm31"], m1=0.0); pJ = hint(JUNO26["dm31"], m1=mJ[0])
line(f"   dependence on m1 (which the 45 deg rule took from the measured small splitting): m1 = 0 gives {p0*1e5:.3f}, "
     f"m1 = {mJ[0]*1e3:.3f} meV gives {pJ*1e5:.3f}; so the reuse moves the hint by {abs(pJ-p0)*1e5:.3f}e-5, 0.3% (not circular in practice)")
derr = r * r * 0.027e-3
line(f"   error from the large splitting alone (+-0.026e-3): +-{derr*1e5:.3f}e-5; fail band at JUNO's present error: "
     f"{(hint(JUNO26['dm31'])-0.078e-5)*1e5:.2f} to {(hint(JUNO26['dm31'])+0.078e-5)*1e5:.2f}e-5 (v1.3.2 printed 7.32 to 7.48)")
line("   how the ratio was found: 3 - 2 sqrt2 = 0.1716 (the midpoint of a side of the mass triangle) was proposed when the")
line("   measured m2/m3 (0.172-0.173) was already known, so it is a hint found after the fact, not a prediction made in advance.")

# 6. TM1
for lab, s13 in (("NuFIT 6.0 s13^2 0.02195", NUFIT6["s13sq"]), ("PDG 2026 s13^2 0.0217", PDG26["s13sq"]), ("v1.3.2 used 0.0220", 0.0220)):
    t12 = 1 - 2 / (3 * (1 - s13))
    line(f"6. TM1 sin^2 theta12 ({lab}): {t12:.4f} vs JUNO 2026 0.3036 +- 0.0064: {z(t12, *JUNO26['s12sq']):+.2f} sigma")

# 7. 0vbb effective mass with TM1, JUNO masses, phases free
s12 = 1 - 2 / (3 * (1 - NUFIT6["s13sq"])); c12 = 1 - s12; s13 = NUFIT6["s13sq"]; c13 = 1 - s13
ph = np.linspace(0, 2 * np.pi, 721)
P1, P2 = np.meshgrid(ph, ph)
mbb = np.abs(c12 * c13 * mJ[0] + s12 * c13 * mJ[1] * np.exp(1j * P1) + s13 * mJ[2] * np.exp(1j * P2))
z0 = abs(c12 * c13 * mJ[0] + s12 * c13 * mJ[1] + s13 * mJ[2] * np.exp(-2j * math.radians(270)))
line(f"7. neutrinoless double beta (TM1, JUNO masses, Majorana phases free): {mbb.min()*1e3:.2f} to {mbb.max()*1e3:.2f} meV "
     f"(v1.3.2: 1.4 to 4.1 with m2 = 8.7); Majorana phases zero and delta = 270 deg: {z0*1e3:.2f} meV")

# 8. strong coupling
pa = 0.1178
line(f"8. alpha_s 0.1178 vs PDG 2026 world average {PDG26['alpha_s'][0]} +- {PDG26['alpha_s'][1]}: {z(PDG26['alpha_s'][0], pa, PDG26['alpha_s'][1]):+.1f} sigma; "
     f"vs lattice {LATTICE_AS[0]} +- {LATTICE_AS[1]}: {z(LATTICE_AS[0], pa, LATTICE_AS[1]):+.1f} sigma   [world average unchanged]")
line("   The inputs of the strength calculation (sin^2 0.23122, 1/alpha(MZ) 127.951, MZ) are unchanged in PDG 2026 at the")
line("   precision used (MZ 91.1876 -> 91.1879 moves alpha_s by far less than 0.0001), so 0.1178 stands.")

# 9. Higgs-mass rule rerun with PDG 2026 top, W and Higgs masses
spec = importlib.util.spec_from_file_location("hg", os.path.join(HERE, "higgs_at_the_grain_2026-11-17.py"))
hg = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(hg)
mt, emt = PDG26["mt"]; mw = PDG26["mW"]; mH, emH = PDG26["mH"]
for far in (False, True):
    c = hg.mh_needed(hg.GRAIN, far, mt=mt, mw=mw)
    dt = (hg.mh_needed(hg.GRAIN, far, mt=mt + emt, mw=mw) - hg.mh_needed(hg.GRAIN, far, mt=mt - emt, mw=mw)) / 2
    da = (hg.mh_needed(hg.GRAIN, far, mt=mt, a3=0.1189, mw=mw) - hg.mh_needed(hg.GRAIN, far, mt=mt, a3=0.1171, mw=mw)) / 2
    dpole = abs(dt) / emt * 0.5
    s = math.sqrt(dt ** 2 + da ** 2 + dpole ** 2 + 0.3 ** 2 + emH ** 2)
    s_nopole = math.sqrt(dt ** 2 + da ** 2 + 0.3 ** 2 + emH ** 2)
    line(f"9. Higgs rule, {'with far side' if far else 'SM only     '}: m_H = {c:.2f} +- {s:.2f} GeV vs PDG 2026 {mH} +- {emH}: "
         f"{z(c, mH, s):+.1f} sigma ({z(c, mH, s_nopole):+.1f} without the 0.5 GeV pole allowance)")
need = __import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda m: hg.mh_needed(hg.GRAIN, True, mt=m, mw=mw) - mH, 165, 180, xtol=1e-4)
line(f"   top pole mass the rule needs for m_H = {mH}: {need:.2f} GeV vs PDG 2026 {mt} +- {emt}: {z(mt, need, emt):.1f} sigma on the top mass alone "
     f"(v1.3.2: 171.63 vs 172.57 +- 0.29, 3.2 sigma)")

# 10. quark-mixing checks
fr, efr = PDG26["first_row"]
line(f"10. first-row unitarity: PDG 2026 {fr} +- {efr} ({(1-fr)/efr:.1f} sigma shortfall). TRPG: exactly 1. Still open, not failed.")
g, eg = PDG26["gamma"]
line(f"    gamma: PDG 2026 {g} +2.7 -2.8 deg; lead 64.1-67.6 deg -> inside")
al = PDG26["alpha"]
line(f"    alpha (direct, PDG 2026): {al[0]} +{al[1]} -{al[2]} deg; the lead needs alpha = 90 deg (fits 89-93): "
     f"90 is {(90-al[0])/al[1]:.1f} sigma above the direct value. Not stated in v1.3.2.")

# 11. tau lifetime
line(f"11. tau lifetime 289.7 fs vs PDG 2026 {PDG26['tau_tau'][0]} +- {PDG26['tau_tau'][1]}: {z(289.7, *PDG26['tau_tau']):+.1f} sigma")
line("=" * 100)
line("Summary: the Higgs rule is under more strain (the top mass it needs is now 3.7 sigma below PDG 2026; the Higgs")
line("mass itself is 1.4 sigma off, 2.1 without the pole allowance). With JUNO's splittings, used throughout from v1.3.2, m2")
line("is 8.6 meV and the sum 0.0591 eV, so the quoted ranges become 8.6-8.7 meV and 0.0591-0.0594 eV, and 0vbb 1.3-4.1 meV.")
line("The down-quark comparison value is corrected (the pass improves). The direct alpha angle puts the alpha = 90 deg lead")
line("1.6 sigma off. Everything else moves in the last digit or not at all.")
