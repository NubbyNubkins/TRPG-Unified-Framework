"""Check (18 Nov 2026): is the 4.6e13 GeV meeting a unification, and could the strengths be equal at the grain instead?
Uses the two-loop runner of force_strengths_final_2026-11-15.py (Machacek & Vaughn 1984; Luo, Wang & Xiao 2003), far-side weak pairs
(spin 1/2, Dirac, 3 copies) at M_L and D partners (3 copies) at M_D. Grain = Planck energy / sqrt2 = 8.63e18 GeV.
(a) the three strengths from M_Z to the grain with M_L = 2.4e8, M_D = 4.6e13 and alpha_s = 0.1178;
(b) can any (M_L, M_D) make all three equal AT THE GRAIN?  (two conditions, two unknowns)
(c) Sakharov-style count of states with the far side."""
import math, importlib.util, os, numpy as np
from scipy.optimize import minimize
up = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fs", os.path.join(up, "force_strengths_final_2026-11-15.py")); fs = importlib.util.module_from_spec(spec)
src = open(os.path.join(up, "force_strengths_final_2026-11-15.py")).read().split("\nhead(")[0]   # definitions only
exec(compile(src, "fs", "exec"), fs.__dict__)
EG = fs.MPL/math.sqrt(2); mp = 0.938272
def at(mu, ML, MD, a3=0.1178):
    _, g = fs.run(a3=a3, thr=fs.thrLD(ML, MD), stop=mu); return [fs.ai(x) for x in g[:3]]
print("(a) 1/alpha1, 1/alpha2, 1/alpha3 (GUT-normalised), alpha_s(MZ) = 0.1178, M_L = 2.4e8, M_D = 4.6e13 GeV")
for mu in (1e3, 2.4e8, 1e12, 4.6e13, 1e15, 1e17, EG):
    a = at(mu, 2.4e8, 4.6e13); print(f"   mu = {mu:9.2e} GeV : {a[0]:6.2f} {a[1]:6.2f} {a[2]:6.2f}   spread {max(a)-min(a):5.2f}")
print("   -> they cross near 4.6e13 GeV and separate again above it; at the grain they are unequal.")
print("\n(b) equality at the grain: minimise the spread at E_grain over log10 M_L, log10 M_D (both between 1e3 and 1e18 GeV)")
def spread(p):
    ML, MD = 10**p[0], 10**p[1]
    if not (3 <= p[0] <= 18 and 3 <= p[1] <= 18): return 1e3
    a = at(EG, ML, MD); return max(a)-min(a)
best = min((minimize(spread, x0, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-4)) for x0 in ([5, 10], [8, 13], [12, 6], [15, 15], [4, 17])), key=lambda r: r.fun)
print(f"   smallest spread at the grain: {best.fun:.2f} (in 1/alpha) at M_L = 10^{best.x[0]:.2f}, M_D = 10^{best.x[1]:.2f} GeV -> {'FAILS' if best.fun > 0.5 else 'possible'}")
a = at(EG, 10**best.x[0], 10**best.x[1]); print(f"   there: {a[0]:.2f} {a[1]:.2f} {a[2]:.2f}")
print("\n(c) states for a Sakharov-style count: SM with 3 hidden partners = 96 fermion + 24 carrier + 4 Higgs = 124;")
print("   far side adds 3 x (10 + 1) Weyl = 33 Weyl = 66 -> 190; with the majoron 191.")
print(f"\n(d) ln(E_grain/m_p) with E_grain = E_Pl/sqrt2: {math.log(EG/mp):.2f}  (44.0 is the unreduced Planck energy, {math.log(fs.MPL/mp):.2f})")
a3g = at(EG, 2.4e8, 4.6e13)[2]; print(f"   simple count 2 pi/(7 alpha_s) at the grain, with the far side: {2*math.pi*a3g/7:.1f}  (SM states alone, earlier: 47.5)")
