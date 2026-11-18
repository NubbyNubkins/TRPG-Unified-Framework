"""Check (18 Nov 2026): running strange-quark mass.
m_s(2 GeV) = 93.5 MeV (PDG 2024) runs up to about 120-140 MeV at mu ~ 0.6-0.8 GeV. In the pressure only (m_s^2 - 4 Delta^2) enters,
so the stars fix 4 Delta^2 - m_s^2; the gap window moves with m_s. Re-check the Bodmer-Witten margins at the shifted gaps."""
import math, importlib.util, os
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("g", os.path.join(here, "ns_ground_state_2026-11-18.py"))
src = open(os.path.join(here, "ns_ground_state_2026-11-18.py")).read().split('print(" B (MeV/fm3)')[0]
g = {"__file__": os.path.join(here, "x.py")}; exec(src, g)
lo, hi = 4*82**2 - 100**2, 4*97**2 - 100**2
print(f"stars fix 4*Delta^2 - m_s^2 = {lo} to {hi} MeV^2 (at B = 60.1)")
for ms in (100, 120, 130, 140):
    d1, d2 = math.sqrt((lo+ms*ms)/4), math.sqrt((hi+ms*ms)/4)
    print(f"  m_s = {ms} MeV -> gap window {d1:.1f} to {d2:.1f} MeV")
print("\nBodmer-Witten at the m_s = 130 window (gap 92-106 MeV):")
for B in (57.1, 60.1, 63.7):
    for gap in (92.0, 106.0):
        s = g["EA_cfl"](B, gap, ms=130.0); u = g["EA_ud"](B, gap)
        print(f"  B = {B}: gap {gap}: strange {s:.1f} MeV ({'stable' if s < 930 else 'NOT'}), up-down {u:.1f} MeV ({'nuclei safe' if u > 934 else 'NUCLEI UNSAFE'})")
