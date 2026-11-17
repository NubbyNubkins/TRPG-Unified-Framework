"""Test: add a value for the third state and perspective across the three scales: quarks with gluons; protons with electrons; the top.
Each scale gets a push/pull (a state that holds partners) and a perspective (a turn of identity). The third scale is the top's.
Credit: Nambu 1989; Miransky, Tanabashi & Yamawaki 1989; Bardeen, Hill & Lindner 1990 (top pairs make the Higgs); Pendleton & Ross 1981, Hill 1981 (fixed point);
Weinberg 1979 & Susskind 1979 (a scaled-up strong force); Strassler & Peskin 1991 (Higgs between tops); Gatto, Sartori & Tonin 1968; Cabibbo 1963;
Kobayashi & Maskawa 1973; Buttazzo et al. 2013; CMS & ATLAS 2025; PDG 2024."""
import math, numpy as np
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
v = 246.21965; mt_pole = 172.57; mH = 125.20; me = 0.51099895e-3; Mp = 0.93827209
head("1. The three scales, each with a push/pull and a perspective (energy-matrix shapes from an earlier test)")
pr("  scale                      push/pull (product-shaped, holds partners)          perspective (rotation, turns identity)")
pr("  1 quarks with gluons       colour-blind total: makes the proton (938 MeV from    colour turn: a gluon turns red to blue")
pr("                             ~9 MeV of quark mass)")
pr("  2 protons with electrons   electric: makes the atom (13.6 eV)                    W turn: neutron to proton, muon to electron")
pr("  3 the top                  THIRD STATE: the Higgs link, charge = mass  THIRD PERSPECTIVE: the family turn, the tilt")
pr("                                                                                   (the W's family-crossing entries, an earlier test)")
pr("  The pattern says: at each scale the push/pull should make something. Scale 1 makes the proton, scale 2 the atom. What does scale 3 make?")

head("2. THIRD STATE, value A: what the Higgs link adds between a top and an antitop (toponium, seen as a threshold excess in 2025)")
aH = mt_pole**2/(4*math.pi*v**2)
for a_s in (0.13, 0.15, 0.17):
    a = 2/(mt_pole*(4/3)*a_s); Bc = ((4/3)*a_s)**2*(mt_pole/2)/2
    shift = aH*4/(a*(2 + mH*a)**2)
    pr(f"  alpha_s {a_s}: colour binding {Bc*1e3:5.0f} MeV; the Higgs link adds {shift*1e3:4.0f} MeV ({100*shift/Bc:.1f}%)   [first order; Strassler & Peskin 1991]")
pr(f"  Higgs-link strength between tops: {aH:.4f} (electric 0.0073, strong ~0.15). Can it bind tops alone? needs strength x (top/2) / Higgs mass > 0.84:")
pr(f"    {aH*mt_pole/2/mH:.3f}. NO. The third state binds nothing on its own (as an earlier test found), but it is real and measurable at the top: ~2% of toponium.")

head("3. THIRD STATE, value B: suppose scale 3 MAKES something, like scales 1 and 2: the Higgs itself as a bound top pair")
pr("  This is a known idea (Nambu; Miransky, Tanabashi & Yamawaki 1989; Bardeen, Hill & Lindner 1990): top pairs bind into the Higgs above some energy L.")
pr("  Condition: the top's Higgs link and the Higgs's own strength both become unlimited at L. Running down gives the top and Higgs masses.")
def rge(t, y):
    g1, g2, g3, yt, la = y; k = 1/(16*math.pi**2)
    return [k*(41/6)*g1**3, k*(-19/6)*g2**3, k*(-7)*g3**3, k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - (17/12)*g1**2),
            k*(24*la**2 + 12*la*yt**2 - 6*yt**4 - 3*la*(3*g2**2 + g1**2) + (3/8)*(2*g2**4 + (g2**2 + g1**2)**2))]
mt_ms = 0.93690*v/math.sqrt(2)
gsol = solve_ivp(lambda t, y: rge(t, list(y) + [0, 0])[:3], [0, math.log(1.22e19/mt_pole)], [0.35830, 0.64779, 1.1666], dense_output=True, rtol=1e-10)
pr(f"  {'L (GeV)':>10s} {'top (running mass)':>20s} {'measured':>9s} {'Higgs':>8s} {'measured':>9s}")
for L in (1e4, 1e6, 1e10, 1e15, 1e19):
    tL = math.log(L/mt_pole); g = gsol.sol(tL); res = []
    for yL, lL in ((6.0, 9.0), (12.0, 36.0)):
        s = solve_ivp(rge, [tL, 0], [*g, yL, lL], rtol=1e-10, atol=1e-12); yt, la = s.y[3, -1], s.y[4, -1]; res.append((yt*v/math.sqrt(2), math.sqrt(max(2*la, 0))*v))
    pr(f"  {L:10.0e} {res[0][0]:12.0f} ({res[1][0]:.0f}) {mt_ms:9.0f} {res[0][1]:8.0f} ({res[1][1]:.0f}) {mH:9.1f}")
pr("  (two starting values shown; the answer forgets them, which is the point of the condition)")
pr("  FAIL: even with L at gravity's grain the top comes out ~35% heavy and the Higgs nearly twice its mass. With L near the top scale itself it is far worse.")
pr("  The scaled-up-strong-force version (Weinberg 1979, Susskind 1979) also fails: it predicts no light Higgs and new heavy mesons near 1-3 TeV; none")
pr("  are seen, and the Higgs found in 2012 behaves like a single elementary state (LHC couplings within ~10%). So scale 3 does NOT make the Higgs.")
pr("  What survives: the Higgs value IS the third scale's 'made thing' in a weaker sense: the push/pull of scale 3 sets every mass, as the colour")
pr("  push/pull sets the proton and the electric push/pull sets the atom. Its size (246 GeV) stays a measured property of space.")

head("4. THIRD PERSPECTIVE: the family turn. A potential value from the masses?")
md, ms, mb, mu_, mc, mtq = 2.678e-3, 53.36e-3, 2.873, 1.242e-3, 0.6278, 162.5    # all at one energy (the top's)
pr(f"  family 1 <-> 2: sqrt(down/strange) = {math.sqrt(md/ms):.4f}; measured turn 0.2250 ({100*(math.sqrt(md/ms)/0.2250 - 1):+.1f}%)   HOLDS (Gatto, Sartori & Tonin 1968)")
pr(f"  family 2 <-> 3: sqrt(strange/bottom) = {math.sqrt(ms/mb):.4f}; measured 0.0418 ({math.sqrt(ms/mb)/0.04183:.1f}x)   FAIL")
pr(f"                  strange/bottom = {ms/mb:.4f} ({ms/mb/0.04183:.2f}x)   FAIL;  sqrt(charm/top) = {math.sqrt(mc/mtq):.4f} ({math.sqrt(mc/mtq)/0.04183:.1f}x)   FAIL")
pr("  So the turn between the first two families follows from the masses; the turn to the third (the top's family) does not, by any single-ladder")
pr("  rule. That is the same finding as an earlier test: the third perspective needs the up and down ladders together. It is the joint rule (next test).")

head("5. Connecting the three scales (the tertiary results): is there a pattern in the scales themselves?")
q3 = Mp/3; pr(f"  electron {me*1e3:.3f} MeV, one third of the proton {q3*1e3:.2f} MeV, top {mt_pole:.2f} GeV")
pr(f"  steps: proton/3 / electron = {q3/me:.0f}; top / (proton/3) = {mt_pole/q3:.0f}. Equal steps would need middle^2 = electron x top:")
pr(f"    (proton/3)^2 / (electron x top) = {q3**2/(me*mt_pole):.3f}  (with the meson home mass 309.6 MeV: {0.3096**2/(me*mt_pole):.3f})")
pr(f"    geometric middle sqrt(electron x top) = {math.sqrt(me*mt_pole)*1e3:.1f} MeV vs 312.8 (-{100*(1 - math.sqrt(me*mt_pole)/q3):.0f}%)")
pr("  NEAR MISS, FAIL as a law (11% on the square). Honest count: this was 1 of about 4 scale patterns I looked at (also proton/electron vs top/proton,")
pr("  1836 vs 184, a factor 10 to 0.2%, but 10 is only special in base ten, so no meaning; and Higgs value / proton/3 = 787, matching nothing). Recorded, not used.")

head("6. Summary")
pr("  HOLDS as a classification: three scales, each with a push/pull and a perspective; the third state is the Higgs link, the third perspective is the family turn.")
pr("  Third state's value at the top: adds ~30 MeV (~2%) to toponium's binding. It cannot bind anything alone.")
pr("  FAILS: scale 3 making the Higgs out of top pairs (top ~35% heavy, Higgs nearly double, at best); the family turn to the third family from one ladder;")
pr("  equal steps between the three scales (11% off). The third perspective needs the joint rule. Scorecard unchanged.")
