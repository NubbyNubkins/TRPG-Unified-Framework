"""Export the best leptogenesis settings from an earlier test (same seed, same search) for an independent ULYSSES run."""
import math, numpy as np, json
__file__ = "leptogenesis_boltzmann_2026-11-14.py"
src = open(__file__).read(); head, tail = src.split("# find candidates")
exec(head)
exec(tail.split("\n",1)[1].split("cands.sort(key=lambda t: -t[0])")[0])
cands.sort(key=lambda t: -t[0])
out = []
for k, (val, pv, s23) in enumerate(cands[:4]):
    be = solve_BE(pv, s23, limit=60)[0]
    d = {"case": k+1, "convention": "R = R23(z1) R13(z2) R12(z3), z = x + i y in radians, R_ij(z): R[i,i]=R[j,j]=cos z, R[i,j]=sin z, R[j,i]=-sin z; "
         "h = sqrt(M) R sqrt(m) U^dagger / v, v = 174 GeV; U = standard PDG PMNS times diag(1, e^{i a21/2}, e^{i a31/2})",
         "m_light_eV": [0.00036, 0.0087, 0.0501], "M_heavy_GeV": list(map(float, M)),
         "x_rad": list(map(float, pv[0:3])), "y_rad": list(map(float, pv[3:6])),
         "delta_rad": float(pv[6]), "a21_rad": float(pv[7]), "a31_rad": float(pv[8]),
         "sin2_12": 0.307, "sin2_13": 0.02195, "sin2_23": s23, "our_eta_B_full_equations": float(abs(be))}
    out.append(d); print(json.dumps(d))
json.dump(out, open("best_leptogenesis_settings_2026-11-14.json", "w"), indent=1)
