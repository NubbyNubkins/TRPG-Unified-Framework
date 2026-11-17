"""Question: since quarks and leptons affect each other and their math is close but not precise, do they push and pull
on each other's values, even slightly?
Part A: the pushes that known physics already contains, with their sizes. Part B: is any of them the right size and direction to close a near miss?
Credit: Gell-Mann & Low 1954 (values change with energy); Cheng, Eichten & Li 1974; Machacek & Vaughn 1984 (how each Higgs link pushes the others);
Arason et al. 1992; Xing, Zhang & Zhou 2008 (masses at common energies); Koide 1982; Sumino 2009 (why is the lepton rule exact for masses as weighed?);
Li & Ma 2006; Rodejohann & Zhang 2011; Jegerlehner 2017 (quark loops in the photon); FLAG 2024; PDG 2024."""
import math, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
L16 = 16*math.pi**2; al = 1/137.036; alq = 1/133.0
def cone(m):
    x = np.sqrt(np.sort(np.asarray(m, float))); a = x.mean(); c = x/a - 1; r2 = (2/3)*(c*c).sum(); return a*a, r2/2, math.degrees(math.atan(math.sqrt(r2/2)))

head("A. The pushes that are already in known physics, and how big each is")
yt, yb, ytau, yc, ys, ymu = 0.937, 2.873/174.10, 1.77693/174.10*0.983, 0.628/174.10, 0.0534/174.10, 0.10566/174.10
pr("  How to read the sizes: 'per step' = per factor of 2.7 in the energy at which the particle is looked at. There are about 32 such steps between")
pr("  the top's energy and the energy where the forces meet, and about 12 between the electron's mass and the top's.")
pr("  A1. THROUGH THE HIGGS, shared by all. Every Higgs link is pulled by the sum of the squares of ALL the links, quarks counted three times (colour).")
pr(f"      Quarks on every lepton link: {3*(yt**2 + yb**2 + yc**2)/L16*100:.2f}% per step (almost all of it the top).")
pr(f"      Leptons on every quark link: {(ytau**2 + ymu**2)/L16*100:.1e}% per step (almost all of it the tau).")
pr("      The SAME for all three rungs of a ladder: it changes a ladder's SIZE and leaves its SHAPE (opening, twist) alone.")
pr(f"      It is lopsided by a factor of about {round(3*yt**2/(ytau**2), -3):,.0f}: the quarks push the leptons; the leptons barely push back.")
pr("  A2. INSIDE A WEAK PAIR, rung by rung. The two members of a pair push on each other's links with OPPOSITE signs: a true push/pull.")
pr(f"      top on bottom: {-1.5*yt**2/L16*100:+.2f}% per step;  bottom on top: {-1.5*yb**2/L16*100:+.1e}%;  charm on strange: {-1.5*yc**2/L16*100:+.1e}%;  tau on its neutrino: {-1.5*ytau**2/L16*100:+.1e}%.")
pr("      This one DOES change shapes: it is why the down ladder's shape and the tilt between the ladders drift above the top's energy.")
pr("      It acts only between partners of one pair. Nothing of this kind joins a quark to a lepton.")
pr("  A3. THROUGH THE ELECTRIC FORCE, each on itself. A charged particle's own electric cloud shifts its mass by an amount set by its charge squared")
pr("      and by its own mass, so it changes shapes. Between the lightest and heaviest rung of each ladder (mass as weighed against one common energy):")
for nme, q2, lo, hi in (("leptons", 1.0, 0.000511, 1.777), ("down ladder", 1/9, 0.00268, 2.873), ("up ladder", 4/9, 0.00124, 162.5)):
    pr(f"         {nme:<12} {3*al*q2/(4*math.pi)*math.log((hi/lo)**2)*100:.2f}%")
pr("  A4. THROUGH THE PHOTON, quark on lepton and lepton on quark. Quark pairs flicker inside the photon that makes a lepton's electric cloud, and")
pr(f"      lepton pairs inside a quark's. Size: about (strength/pi)^2 = {(al/math.pi)**2:.1e} times a few, so 1e-5 to 1e-4 of the mass, differing rung to rung.")
pr("  -> So YES: quarks and leptons do push on each other's values, and it is small. Directly (A4): about a part in 100,000. Through the Higgs (A1):")
pr("     large but shape-blind. The only strong rung-by-rung push/pull (A2) is between the two members of one pair, quark with quark.")

head("B. Is any push the right size, and the right direction, to close a near miss?")
Mp = np.array([0.51099895, 105.6583755, 1776.93]); mprot = 938.27208
def lep_at(mu): return Mp*(1 - (al/math.pi)*(1 + 0.75*np.log(mu*mu/Mp**2)))          # one-loop electric cloud, mass as weighed -> mass at energy mu (MeV)
pr("  B1. The lepton cone, and its size against one third of the proton (313.85 against 312.76 MeV, 0.35% apart):")
for nme, m in (("masses as weighed", Mp), ("all at the proton's own energy", lep_at(mprot)), ("all at the Z's energy", lep_at(91187.6))):
    a2, rs, op = cone(m); pr(f"         {nme:<32} size {a2:.2f} MeV ({100*(a2/(mprot/3) - 1):+.2f}% from proton/3)   opening {op:.3f} deg   relationship/self {rs:.5f}")
mu = np.exp(np.linspace(math.log(200), math.log(5000), 4000)); sz = np.array([cone(lep_at(x))[0] for x in mu]); hit = mu[np.argmin(abs(sz - mprot/3))]
pr(f"      The size equals proton/3 exactly at a common energy of about {hit:.0f} MeV. But there the opening is {cone(lep_at(hit))[2]:.3f} deg, no longer 45.")
pr("      So the electric push is the RIGHT SIZE to move the 0.35%, but it cannot make both numbers exact at once: the cone is exact only for masses as")
pr("      weighed, the size match only near 1.3 GeV, and nothing in TRPG picks 1.3 GeV. FAIL as a fix. The 0.35% stays a near miss.")
pr("      The real puzzle it exposes (Sumino 2009): the electric push should spoil the balance (relationship/self) by 0.35%, yet for masses as weighed it holds to 0.001%.")
D = np.array([2.678, 53.36, 2873.0]); Uq = np.array([1.242, 627.8, 162500.0]); mu0 = 162500.0
pr("  B2. The quark ladders against 5/4 and 5/3. Put the quarks on the leptons' footing: add each quark's own electric cloud down to its own")
pr("      mass (or, since a light quark is never free, only down to 300 or 1000 MeV):")
for cut in (0.0, 300.0, 1000.0):
    fD = 1 + (3*alq/9/(4*math.pi))*np.log(mu0**2/np.maximum(D, cut)**2); fU = 1 + (3*alq*4/9/(4*math.pi))*np.log(mu0**2/np.maximum(Uq, cut)**2)
    pr(f"         cloud counted down to {'its own mass' if cut == 0 else f'{cut:.0f} MeV':<13}: down ladder {cone(D*fD)[1]:.4f} (was {cone(D)[1]:.4f}; 5/4 = 1.2500)    up ladder {cone(Uq*fU)[1]:.4f} (was {cone(Uq)[1]:.4f}; 5/3 = 1.6667)")
pr("      Both move AWAY from the simple values, and by less than a fifth of the gap. Too small AND the wrong direction. FAIL.")
pr("  B3. The larger misses: charm in the colour rule (+6.8%), one strength for all forces (10.6%), the top against the full link (0.9%).")
pr(f"      The quark-lepton pushes are at most {3*al/(4*math.pi)*math.log((1.777/0.000511)**2)*100:.1f}% (A3) or shape-blind (A1). They cannot reach the first two.")
pr("      The top: its link is itself pushed about 5% per step by its own size and the strong force; 'equal to 1' would hold at one energy only. Already FAIL (34.4).")

head("C. Summary")
pr("  YES  : quarks and leptons push on each other's values. Three real routes, sizes above. Standard physics; nothing added.")
pr("  SHAPE: the big push (through the Higgs) cannot change a ladder's shape; the pushes that can are small, or act only inside one weak pair.")
pr("  FAIL : as the reason the math is 'close but not precise'. No miss is closed: the quark ladders move the wrong way and too little, and the")
pr("         lepton size can be matched to the proton only by spoiling the 45 degrees.")
pr("  KEPT : where a rule holds exactly matters. The lepton cone is exact for masses as weighed and for no common energy. Any quark rule must say")
pr("         which footing it uses BEFORE it is tested; the electric push alone moves the answer by tenths of a per cent.")
