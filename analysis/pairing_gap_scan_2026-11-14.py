"""Neutron-star maximum mass, radius and tidal deformability against the bag
constant B and the quark pairing gap, using the model in
code/neutron_stars_no_knobs_2026-11-14.py. Run from anywhere; the code folder
may be given as the first argument."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "code")
sys.path.insert(0, CODE); os.chdir(CODE)
import importlib, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    m = importlib.import_module("neutron_stars_no_knobs_2026-11-14")
GRID = [(B, g) for B in (57.1, 60.1, 63.7) for g in (70.0, 80.0, 90.0)]
GRID += [(57.1, 75.0), (57.1, 85.0), (60.1, 85.0), (60.1, 95.0), (60.1, 100.0),
         (63.7, 95.0), (63.7, 105.0), (63.7, 110.0)]
for B, gap in GRID:
    M, L, R, a4, es = m.sequence(B, 3.0, gap)
    print(f"B={B} gap={gap} Mmax={M:.2f} R14={R:.2f} L14={L:.0f}", flush=True)
