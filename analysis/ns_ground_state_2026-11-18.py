"""Check (18 Nov 2026): what kind of star the neutron-star calculation is, and the Bodmer-Witten check.
The TOV run (neutron_stars_no_knobs_2026-11-14.py) applies the quark-matter pressure up to a bare surface with a density jump
(neutron_star_tidal_fix_2026-11-13.py, self-bound surface correction): it is a self-bound strange quark star with no hadronic crust.
That is only consistent if three-flavour (strange) quark matter is the true ground state (E/A < 930 MeV, Bodmer 1971; Witten 1984)
while two-flavour (u,d) quark matter is NOT (E/A > ~934 MeV, else ordinary nuclei would convert; Farhi & Jaffe 1984).
Here: E/A at zero pressure for the same EoS (B, a4 from running alpha_s at X=3, m_s = 100 MeV, gap) and for u,d matter
(two flavours, no strange, paired 2SC-like with the same gap as a rough upper bound on binding; Alford, Rajagopal, Reddy & Wilczek 2001)."""
import math, numpy as np, importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__)); up = here
spec = importlib.util.spec_from_file_location("ns", os.path.join(up, "neutron_stars_no_knobs_2026-11-14.py")); ns = importlib.util.module_from_spec(spec); spec.loader.exec_module(ns)
hc = 197.3269
def EA_cfl(B, gap, X=3.0, ms=100.0):
    mu = np.linspace(200, 600, 20001); a4 = np.array([1-2*ns.alpha_s(X*m/1000)/math.pi for m in mu])
    p = (3/(4*math.pi**2))*(a4*mu**4 - (ms**2-4*gap**2)*mu**2)/hc**3 - B
    i = np.argmax(p > 0); return 3*mu[i]          # at p = 0, E/A = 3 mu (baryon = 3 quarks)
def EA_ud(B, gap, X=3.0):
    # two massless flavours, 2SC pairing of 2 of 3 colours: Omega = -(2/(4pi^2)) a4 mu^4 ... use p = (1/(2pi^2)) a4 mu^4 *(2/... )
    # free u,d with charge neutrality ignored (symmetric): p = (Nc Nf/(12 pi^2)) a4 mu^4 = (6/(12pi^2)) a4 mu^4 ; 2SC adds (Delta^2 mu^2)/pi^2
    mu = np.linspace(200, 600, 20001); a4 = np.array([1-2*ns.alpha_s(X*m/1000)/math.pi for m in mu])
    p = ((6/(12*math.pi**2))*a4*mu**4 + gap**2*mu**2/math.pi**2)/hc**3 - B
    i = np.argmax(p > 0); return 3*mu[i]
print(" B (MeV/fm3)  gap (MeV)  E/A strange matter  E/A u,d matter   strange matter stable?  ud matter unstable?")
for B in (57.1, 60.1, 63.7):
    for gap in (0.0, 82.0, 90.0, 97.0):
        s = EA_cfl(B, gap); u = EA_ud(B, gap)
        print(f"  {B:5.1f}      {gap:5.1f}      {s:7.1f} MeV          {u:7.1f} MeV        {'yes' if s < 930 else 'NO'}                    {'yes' if u > 934 else 'NO'}")
