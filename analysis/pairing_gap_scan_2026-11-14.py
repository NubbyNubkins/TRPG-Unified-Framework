import sys,os
sys.path.insert(0, sys.argv[1]); os.chdir(sys.argv[1])
import importlib; m=importlib.import_module("neutron_stars_no_knobs_2026-11-14")
for B in (57.1,60.1,63.7):
  for gap in (70.0,80.0,90.0):
    M,L,R,a4,es=m.sequence(B,3.0,gap)
    print(f"B={B} gap={gap} Mmax={M:.2f} R14={R:.2f} L14={L:.0f}",flush=True)
