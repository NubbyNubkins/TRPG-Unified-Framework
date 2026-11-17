"""Ideas 4, 1 and 3 (14 Nov 2026). Credits in ideas_for_last_four_2026-11-14.py."""
import math, itertools, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
f = 6.3e10; Mpl = 1.22e19; Mred = 2.435e18
head("IDEA 4. The majoron's mass from gravity leaking into the kept spot (Planck-suppressed breaking of matter number)")
pr("  A term of power n breaks the spot's symmetry: V = c (f/sqrt2)^n / M^(n-4) cos(n J/f), so m^2 = c n^2 (f/sqrt2)^n / (M^(n-4) f^2).")
pr("  Target for all of the dark matter: 1.4 keV / theta^4; theta from 0.5 to 2 gives 0.09 to 23 keV.")
lo, hi = 1.42e3/16, 1.42e3*16
hits = []
for Mname, M in (("Planck mass", Mpl), ("reduced Planck mass", Mred)):
    for n in range(5, 13):
        m = math.sqrt(n**2*(f/math.sqrt(2))**n/(M**(n-4)*f**2))*1e9
        tag = "  <== lands" if lo <= m <= hi else ""
        if tag: hits.append((Mname, n, m))
        pr(f"   {Mname:20s} n = {n:2d}: m = {m:9.2e} eV{tag}")
pr(f"  Lands for: {', '.join(f'{a}, n = {b} ({c:.2g} eV)' for a, b, c in hits) or 'none'}.")
step = math.log10(math.sqrt(Mpl/(f/math.sqrt(2))))
pr(f"  Fair-test check: each step in n changes the mass by a factor {10**step:.1e}; the target window is a factor {hi/lo:.0f} wide.")
pr(f"  So for any settling value, about {100*math.log10(hi/lo)/step:.0f}% of the time SOME power lands, per choice of mass scale (two choices here).")
pr("  And every lower power (n = 5, 6, 7) must be switched off by hand; nothing in the matrix does that yet. Landing is not evidence.")

head("IDEA 1. The far-side mass as a balance point between the theory's own scales")
pr("  Target: 2.6e7 GeV (one loop); two-loop terms move it by roughly a factor 3, so the window is 0.9e7 to 8e7 GeV.")
S = {"electron": 0.000511, "muon": 0.10566, "tau": 1.77686, "proton": 0.93827, "top": 172.57, "Higgs value": 246.22, "Higgs mass": 125.25,
     "Z": 91.19, "settling value": f, "Planck": Mpl, "reduced Planck": Mred}
lo1, hi1 = 0.9e7, 8e7; hits = []; vals = []
names = list(S)
for a, b in itertools.combinations_with_replacement(names, 2):
    x = math.sqrt(S[a]*S[b]); vals.append(x)
    if lo1 <= x <= hi1: hits.append(f"sqrt({a} x {b}) = {x:.2e}")
for a, b in itertools.permutations(names, 2):
    x = S[a]**2/S[b]; vals.append(x)
    if lo1 <= x <= hi1: hits.append(f"{a}^2/{b} = {x:.2e}")
pr(f"  {len(vals)} simple balance points (geometric means and seesaws) from {len(S)} scales. Inside the window: {len(hits)}")
for h in hits: pr("   " + h)
rng = np.random.default_rng(3); arr = np.log10(np.array(vals)); cnt = 0
for i in range(2000):
    c = rng.uniform(2, 12); cnt += np.any((arr > c) & (arr < c + math.log10(hi1/lo1)))
pr(f"  Control: a window this wide placed anywhere from 1e2 to 1e14 GeV catches at least one such combination {100*cnt/2000:.0f}% of the time.")

head("IDEA 3. The pop is the settling itself: the settling piece drives the start and hands its energy to the hidden partners")
pr("  The settling piece (mass m_s ~ settling value) decays to hidden partners with the coupling y = M_N/f that gives them their mass.")
pr("  For no hot majorons the space it heats must stay well below the partners' mass: T_reheat < M_N/10.")
pr("  Decay rate >= y^3 f / (4 pi) (it must be heavy enough to make the pair), so T_reheat / M_N ~ 0.4 sqrt(y M_red / (4 pi f)):")
for y in (1e-1, 1e-2, 1e-4, 1e-6, 1e-8):
    MN = y*f; G = y**3*f/(4*math.pi); T = (90/(math.pi**2*106.75))**0.25*math.sqrt(G*Mred)
    eps = 3*MN*0.0503e-9/(16*math.pi*174.1**2)
    pr(f"   y = {y:7.0e}: partner mass {MN:8.1e} GeV, T_reheat/M_N = {T/MN:9.2e}, largest CP asymmetry (Davidson-Ibarra) {eps:8.1e}")
pr("  Cool enough needs y below about 1e-8, i.e. a partner near 1e3 GeV, whose CP asymmetry (about 1e-13) is far too small: the matter excess")
pr("  needs (T/m) x asymmetry of about 1.6e-10, with T/m below 1. Either the start is too hot (hot majorons) or too weak (no excess).")
