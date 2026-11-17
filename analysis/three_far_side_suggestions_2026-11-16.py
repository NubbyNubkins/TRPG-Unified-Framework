# Tests of the three suggestions for the far-side rules
import numpy as np, math
print("=== 1. MMG partner ===")
# Far side = the 10 of the 27 after fixing the spot a: 27 -> 1(a) + 10(b,c,x) + 16(y,z). det term a(bc-|x|^2):
# every one of the 10 gets weight <a>: MMG (b,c,x0,xu: 4 colourless) and D (6 coloured parts of x) EQUAL.
print("Minimal matrix (determinant cubic, Jordan/Albert): 10 parts all weigh <a>; MMG and D tied -> missing partner FAILS (= an earlier test A-i).")
# Dimopoulos-Wilczek form: far-side weight proportional to each part's matter number (B-L)
parts={"D (colour triplet, charge -1/3)":-2/3,"D-bar":+2/3,"MMG doublet":0.0,"MMG doublet-bar":0.0}
for k,v in parts.items(): print(f"  {k:32s} matter number {v:+.3f} -> direct weight {'F x '+format(abs(v),'.3f') if v else 'ZERO: only through the grain, F^2/grain'}")
F=4.56e13; G=8.63e18; print(f"  MMG = F^2/grain = {F*F/G:.3e} GeV (needed 2.2e8 two-loop strengths) -> automatic, no tuning")
print("  Label: LEAD. Structure follows exactly (MMG has matter number 0); assumption = the far side's weight points along matter number.")
print("\n=== 2. D link at a fixed point ===")
a=1/38.6; g3=math.sqrt(4*math.pi*a)
for c in (3.0,4.5,6.0):
    lam=math.sqrt(8*g3**2/c); print(f"  self-coefficient {c}: fixed-point link coupling {lam:.3f}; own coupling needed for 34.5: {lam/34.5:.4f} (alpha at meeting {a:.4f}, ratio {lam/34.5/a:.2f})")
print("  A fixed point sets ONE coupling (the link, ~0.6-0.9), not the ratio. Own ~ alpha_meeting is an earlier test B-iii: chi2 31.8 FAIL.")
print("\n=== 3. Angles = 45 + push of each family's strongest link ===")
def op(m):
    m=np.array(m,float); s=np.sqrt(m); return math.degrees(math.acos(s.sum()/math.sqrt(3*m.sum())))
up=[1.23,620,168300]; dn=[2.67,53.2,2839]
t=math.log(8.63e18/91.19)
for yt2 in (0.3,0.4,0.5):
    f=math.exp(1.5*yt2*t/(16*math.pi**2))
    print(f"  top-link running MZ->grain (avg yt^2={yt2}): light/top ratio x{f:.3f}; up opening {op(up):.2f} -> {op([up[0]*f,up[1]*f,up[2]]):.2f} at the grain")
# what ratio change would bring up to 45?
for f in (2,5,10,20,50):
    print(f"   ratio x{f}: up opening {op([up[0]*f,up[1]*f,up[2]]):.2f}")
print("  Running from the top's link moves the up angle by ~0.3 deg, not 7: FAIL as the cause. Ordering of the four families by strongest link: pattern only (post hoc, ~4% by chance).")
