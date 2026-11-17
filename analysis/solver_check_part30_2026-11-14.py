import math, numpy as np
from scipy.integrate import solve_ivp
__file__ = "leptogenesis_boltzmann_2026-11-14.py"
src = open(__file__).read()
head, tail = src.split("# find candidates")
exec(head)
cand_code = tail.split("\n",1)[1].split("cands.sort(key=lambda t: -t[0])")[0]
exec(cand_code)
cands.sort(key=lambda t: -t[0])
for idx in range(4):
    val, pv, s23 = cands[idx]
    out = []
    for meth, rt in (("Radau",1e-6),("Radau",1e-8),("BDF",1e-8),("LSODA",1e-9)):
        import re
        g = {"__file__": __file__}
        code = head.replace('method="Radau", rtol=1e-6', f'method="{meth}", rtol={rt}')
        exec(code, g)
        r = g["solve_BE"](pv, s23, limit=60.0)[0]
        out.append(f"{meth} {rt:g}: {abs(r):.2e}")
    print(idx, s23, f"approx {val:.2e} |", " ; ".join(out), flush=True)
