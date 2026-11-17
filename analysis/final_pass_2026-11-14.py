"""Final pass on the open items (14 Nov 2026).
C. Is the matter-excess misalignment the same size as the neutrino-mixing misalignment (TM1, 14.9 deg = 0.26 rad; an earlier test)?
F. Which delta_CP values allow the observed matter excess? (approximate flavoured formula, then full rate equations)
B. A principle for the lepton 45 deg rule: equal chance of 'self' vs 'relationship' (maximum uncertainty of a yes/no split).
Uses the functions of leptogenesis_boltzmann_2026-11-14.py."""
import math, numpy as np
__file__ = "leptogenesis_boltzmann_2026-11-14.py"
exec(open(__file__).read().split("# find candidates")[0])
rng = np.random.default_rng(31)
obs = 6.1e-10

print("="*78); print("C. One misalignment for both? (|angle| = 0.26 rad, the TM1 value)"); print("="*78)
def eta_single(plane, phi, d, a21, a31, s23, mag=0.26):
    z = mag*np.exp(1j*phi)
    R = rot(*plane, z); U = pmns(0.307, 0.02195, s23, d, a21, a31)
    return abs(eta(R, U)[0]), R, U
for plane in [(1, 2), (0, 2), (0, 1)]:
    best = (0, None)
    for _ in range(20000):
        args = (rng.uniform(0, 2*math.pi), rng.uniform(0, 2*math.pi), rng.uniform(0, 2*math.pi), rng.uniform(0, 2*math.pi))
        s23 = 0.47 if rng.random() < 0.5 else 0.561
        val = eta_single(plane, *args, s23)[0]
        if val > best[0]: best = (val, (args, s23))
    (phi, d, a21, a31), s23 = best[1]
    # full equations at that point: express as pv for solve_BE is not possible (single plane) -> direct BE with custom R
    print(f"  plane {plane}: best approximate eta_B = {best[0]:.2e} (phase {math.degrees(phi)%360:.0f} deg, delta {math.degrees(d)%360:.0f} deg)")

print("\n"+"="*78); print("F. Which delta_CP allows the observed excess?"); print("="*78)
def pv_eta(pv, s23): return params_eta(pv, s23)
rows = []
for dbin in range(0, 360, 30):
    best = (0, None, None)
    for s23 in (0.47, 0.561):
        pts = []
        for _ in range(1500):
            pv = np.concatenate([rng.uniform(0, 2*math.pi, 3), rng.uniform(-0.5, 0.5, 3),
                                 [math.radians(rng.uniform(dbin, dbin+30))], rng.uniform(0, 2*math.pi, 2)])
            pts.append((pv_eta(pv, s23), pv))
        pts.sort(key=lambda t: -t[0]); val, pv = pts[0]; step = 0.3
        for it in range(400):
            tr = pv + rng.normal(scale=step, size=9); tr[3:6] = np.clip(tr[3:6], -1, 1)
            tr[6] = math.radians(dbin + (math.degrees(tr[6]) - dbin) % 30)   # keep delta inside the bin
            tv = pv_eta(tr, s23)
            if tv > val: val, pv = tv, tr
            if it % 100 == 99: step *= 0.5
        if val > best[0]: best = (val, pv, s23)
    be, K = solve_BE(best[1], best[2])
    rows.append((dbin, best[0], abs(be)))
    print(f"  delta {dbin:3d}-{dbin+30:3d} deg: best approximate {best[0]:.1e}, full equations {abs(be):.1e} "
          f"-> {'ALLOWS observed' if abs(be) >= obs else ('within 2x' if abs(be) >= obs/2 else 'too small')}")

print("\n"+"="*78); print("B. 45 deg as 'equal chance of self vs relationship'"); print("="*78)
for name, lean in [("leptons", 44.99991), ("down quarks (own scale)", 47.535), ("up quarks (own scale)", 51.195),
                   ("down quarks (common scale, P2)", 48.06), ("up quarks (common scale, P2)", 52.15)]:
    p = math.cos(math.radians(lean))**2
    H = -(p*math.log2(p) + (1-p)*math.log2(1-p))
    print(f"  {name:32s} chance of the 'fully related' (democratic) view = {p:.5f}; uncertainty = {H:.5f} bits (max 1)")
print("  Leptons sit at the maximum (1.00000 bits). Quarks are below it; colour (the strong force's extra relationships)")
print("  would be the thing pulling them off. Principle, not derivation.")
