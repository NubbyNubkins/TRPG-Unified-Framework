"""The cyclic reading (simple to complex; singularity takes in too much and releases; black holes consume all; rebirth).
Credit: Clausius 1865; Boltzmann 1877; Schroedinger 1944; Prigogine 1977; Bekenstein 1973; Hawking 1974/75; Tolman 1934; Penrose 2010 (CCC);
Steinhardt & Turok 2002; Smolin 1992; Poplawski 2010; Egan & Lineweaver 2010; Riess et al. 1998; Perlmutter et al. 1999; DESI 2024/2025."""
import math
def pr(*a): print(*a, flush=True)
G, c, hbar, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.381e-23; Msun = 1.989e30; yr = 3.156e7; Gly = 9.461e24
pr("1. Where the universe's entropy is (Egan & Lineweaver 2010, in units of k):")
for n, s in (("stars", "1e81"), ("cosmic background light", "2e89"), ("relic neutrinos", "2e89"), ("stellar black holes", "~1e97"), ("supermassive black holes", "3e104")):
    pr(f"  {n:26s} {s}")
SBH = lambda M: 4*math.pi*G*M**2/(hbar*c)   # S/k
pr(f"  One Sun-mass black hole: S/k = {SBH(Msun):.1e}; the Sun as a star ~1e58. Gobbling by black holes IS the growth of disorder.")
pr("\n2. Bigger black holes leak less (Hawking): leak temperature and lifetime")
for M in (1, 4.3e6, 1e11):
    T = hbar*c**3/(8*math.pi*G*M*Msun*kB); t = 5120*math.pi*G**2*(M*Msun)**3/(hbar*c**4)/yr
    pr(f"  {M:9.2e} Suns: temperature {T:.1e} K, lifetime {t:.1e} years")
pr("  So taking in more makes a black hole MORE contained, not less; it becomes uncontainable only at the small end (final burst).")
H0 = 67.4e3/3.086e22; RH = c/H0; M = c**3/(2*G*H0)
pr(f"\n3. Universe as a black hole: mass inside the Hubble radius at critical density {M:.1e} kg has a horizon {2*G*M/c**2/Gly:.1f} Gly = Hubble radius {RH/Gly:.1f} Gly")
pr("  This is an exact identity for a flat universe (r_s = c/H), so it is why the idea recurs, not evidence for it (Pathria 1972; Poplawski 2010).")
pr("\n4. Can everything recombine? With constant dark energy (TRPG an earlier test: a kept amount) the expansion accelerates forever.")
pr("  Galaxies beyond ~16 Gly (today's cosmic event horizon) can never reach us; only the Local Group merges. Black holes then evaporate")
pr("  (lifetimes above). End state: cold, dilute light, not one singularity. A recollapse needs dark energy to change sign; TRPG's own rule says it does not.")
pr("  The closest working cycle is Penrose's conformal cyclic cosmology: after black holes evaporate only massless light is left, nothing")
pr("  measures size, and the cold vast end looks like a new hot beginning. No recollapse needed; its claimed evidence ('Hawking points',")
pr("  An, Meissner & Penrose 2018) is disputed.")
