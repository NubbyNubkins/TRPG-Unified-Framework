"""Three-loop check of the Higgs-mass rule (18 November 2026), fixed before running.
Question: does replacing two-loop running (plus a 0.3 GeV allowance for the missing terms) by the three-loop result move the
prediction, and in which direction? The answer is accepted whichever way it goes.
Method: Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio & Strumia (2013, JHEP 12, 089) give the Standard Model's lambda at the
Planck mass at three loops (NNLO matching) as a linear fit:
    lambda(M_Pl) = -0.0143 - 0.0066 (M_t - 173.34) + 0.0018 (alpha_s - 0.1184)/0.0007 + 0.0029 (M_h - 125.15).
First the two-loop code of higgs_at_the_grain_2026-11-17.py is checked against the same fit (value and slopes), then the
difference between three and two loops at the Planck mass is added to the far-side result at the framework's grain.
Inputs: PDG 2026 (M_t 172.60 +- 0.27, M_h 125.13 +- 0.11, M_W 80.3625, alpha_s 0.1180 +- 0.0009).
"""
import os, io, contextlib, importlib.util
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("hg", os.path.join(HERE, "higgs_at_the_grain_2026-11-17.py"))
hg = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(hg)

def lam3(mt, a3, mh): return -0.0143 - 0.0066 * (mt - 173.34) + 0.0018 * (a3 - 0.1184) / 0.0007 + 0.0029 * (mh - 125.15)
ref = dict(mt=173.34, a3=0.1184, mw=80.384)
l2 = hg.lam_at(hg.MPL, 125.15, False, **ref)[4]
d_mt = (hg.lam_at(hg.MPL, 125.15, False, mt=173.84, a3=0.1184, mw=80.384)[4] - hg.lam_at(hg.MPL, 125.15, False, mt=172.84, a3=0.1184, mw=80.384)[4])
d_mh = (hg.lam_at(hg.MPL, 125.65, False, **ref)[4] - hg.lam_at(hg.MPL, 124.65, False, **ref)[4])
d_as = (hg.lam_at(hg.MPL, 125.15, False, mt=173.34, a3=0.11875, mw=80.384)[4] - hg.lam_at(hg.MPL, 125.15, False, mt=173.34, a3=0.11805, mw=80.384)[4])
print("1. Two-loop code against the three-loop fit at the Buttazzo point (M_Planck):")
print(f"   lambda: two-loop {l2:+.4f}, three-loop fit -0.0143")
print(f"   slopes per GeV of M_t: {d_mt:+.4f} (fit -0.0066); per GeV of M_h: {d_mh:+.4f} (fit +0.0029); per 0.0007 of alpha_s: {d_as:+.4f} (fit +0.0018)")

mt, mh, a3, mw = 172.60, 125.13, 0.1180, 80.3625
need3 = brentq(lambda m: lam3(mt, a3, m), 100, 160)
need2 = hg.mh_needed(hg.MPL, False, mt=mt, a3=a3, mw=mw)
shift = need3 - need2
print("\n2. Standard Model only, lambda = 0 at the Planck mass, PDG 2026 inputs:")
print(f"   two-loop code {need2:.2f} GeV; three-loop fit {need3:.2f} GeV; three minus two = {shift:+.2f} GeV")
far2 = hg.mh_needed(hg.GRAIN, True, mt=mt, a3=a3, mw=mw)
far3 = far2 + shift
dt = (hg.mh_needed(hg.GRAIN, True, mt=mt + 0.27, mw=mw) - hg.mh_needed(hg.GRAIN, True, mt=mt - 0.27, mw=mw)) / 2
da = (hg.mh_needed(hg.GRAIN, True, mt=mt, a3=0.1189, mw=mw) - hg.mh_needed(hg.GRAIN, True, mt=mt, a3=0.1171, mw=mw)) / 2
dpole = abs(dt) / 0.27 * 0.5
loops = 0.1   # what remains after the three-loop step (Buttazzo et al. quote about 0.1 GeV from the four-loop QCD term)
s = (dt**2 + da**2 + dpole**2 + loops**2 + 0.11**2) ** 0.5
s0 = (dt**2 + da**2 + loops**2 + 0.11**2) ** 0.5
print("\n3. With the far side, at the framework's grain (Planck/sqrt2), the three-loop step added:")
print(f"   two-loop {far2:.2f} GeV -> three-loop {far3:.2f} +- {s:.2f} GeV (pole allowance included) vs 125.13 +- 0.11: {(far3-mh)/s:+.1f} sigma;"
      f" without the allowance +- {s0:.2f}: {(far3-mh)/s0:+.1f} sigma")
need_t = brentq(lambda m: hg.mh_needed(hg.GRAIN, True, mt=m, mw=mw) + shift - mh, 165, 180, xtol=1e-4)
print(f"   top pole mass the rule needs: {need_t:.2f} GeV vs 172.60 +- 0.27: {(172.60-need_t)/0.27:.1f} sigma on the top mass alone")
print("\nVerdict line fixed before running: within 2 sigma PASS, beyond 3 sigma FAIL, between PARTIAL.")
v = abs(far3 - mh) / s; v0 = abs(far3 - mh) / s0
lab = lambda z: "PASS" if z < 2 else ("FAIL" if z > 3 else "PARTIAL")
print(f"   with the allowance: {lab(v)}; without it: {lab(v0)}")
