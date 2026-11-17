"""Idea tested: the triangle may be one face of a many-sided pyramid; from the side you see one face; a pyramid is a
containment shape that nothing escapes. Tests on the quark turn (CKM): the six unitarity triangles; their common area; the volume that joins
the ladders and the turn; and whether the figure closes (nothing escapes). Credit: Jarlskog 1985; Bjorken & Dunietz 1987; Wolfenstein 1983;
PDG 2024; Belfatto, Beradze & Berezhiani 2020 and Branco et al. 2021 (a heavy down-type partner and the first-row shortfall)."""
import numpy as np, math, itertools
def pr(*a): print(*a, flush=True)
lam, A, rb, eb = 0.22501, 0.826, 0.1591, 0.3523
s12 = lam; s23 = A*lam**2; s13e = A*lam**3*(rb + 1j*eb)*math.sqrt(1 - A**2*lam**4)/(math.sqrt(1 - lam**2)*(1 - A**2*lam**4*(rb + 1j*eb)))
s13 = abs(s13e); d = -np.angle(s13e)
c12, c23, c13 = (math.sqrt(1 - x**2) for x in (s12, s23, s13)); e = np.exp(1j*d)
V = np.array([[c12*c13, s12*c13, s13/e], [-s12*c23 - c12*s23*s13*e, c12*c23 - s12*s23*s13*e, s23*c13],
              [s12*s23 - c12*c23*s13*e, -c12*s23 - s12*c23*s13*e, c23*c13]])
pr("Six faces: each pair of rows or columns closes a triangle in the complex plane; area of each:")
for kind in ("columns", "rows"):
    for i, j in itertools.combinations(range(3), 2):
        sides = [V[k, i]*np.conj(V[k, j]) for k in range(3)] if kind == "columns" else [V[i, k]*np.conj(V[j, k]) for k in range(3)]
        a, b = sides[0], sides[1]; area = abs((a.conjugate()*b).imag)/2
        pr(f"  {kind} {i+1}-{j+1}: sides {', '.join(f'{abs(s):.2e}' for s in sides)}; area {area:.3e}")
J = (V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1])).imag
pr(f"All six areas equal J/2 = {abs(J)/2:.3e} (Jarlskog 1985). Some faces are long and thin, one is nearly equilateral (the one measured at B factories).")
mu2 = np.array([1.242, 627.8, 162500.0])**2; md2 = np.array([2.678, 53.36, 2873.0])**2
Mu = np.diag(mu2); Md = V @ np.diag(md2) @ V.conj().T
C = Mu @ Md - Md @ Mu; det = np.linalg.det(C)
prod = (mu2[2]-mu2[1])*(mu2[2]-mu2[0])*(mu2[1]-mu2[0])*(md2[2]-md2[1])*(md2[2]-md2[0])*(md2[1]-md2[0])
pr(f"The joining volume: det[up ladder, down ladder] = {det.imag:.3e} i (MeV^12); formula 2 J x (mass-squared gaps) = {2*J*prod:.3e}: equal.")
pr("  So the ladders and the turn form one perspective-free object: a determinant (a volume) whose faces are the six triangles.")
row1 = abs(V[0,0])**2 + abs(V[0,1])**2 + abs(V[0,2])**2
pr(f"Closure (nothing escapes) holds by construction in the 3-family fit ({row1:.6f}). Measured directly (PDG 2024): |Vud|^2+|Vus|^2+|Vub|^2 =")
vud, vus, vub = 0.97367, 0.22431, 0.00382; s = vud**2 + vus**2 + vub**2; err = math.sqrt((2*vud*0.00032)**2 + (2*vus*0.00085)**2)
pr(f"  {s:.5f} +- {err:.5f}: short of 1 by {1-s:.5f} ({(1-s)/err:.1f} sigma), the 'Cabibbo angle anomaly'. A heavy down-type partner mixing with")
pr(f"  the down quark at an angle of about {math.sqrt(1-s):.3f} would account for it (Belfatto, Beradze & Berezhiani 2020; Branco et al. 2021).")

pr("\nCan the six triangles close into a solid? (diamond = triangular bipyramid: 6 faces, 9 edges, each edge shared by two faces)")
tris = {}
for kind in ("col", "row"):
    for i, j in itertools.combinations(range(3), 2):
        tris[f"{kind}{i+1}{j+1}"] = sorted(abs(V[k, i]*V[k, j]) if kind == "col" else abs(V[i, k]*V[j, k]) for k in range(3))
L = sorted((l, n) for n, t in tris.items() for l in t)
ratios = [L[2*k+1][0]/L[2*k][0] for k in range(9)]
pr("  18 sides sorted and paired as closely as possible (the best possible pairing for gluing):")
for k in range(9): pr(f"    {L[2*k][0]:.3e} ({L[2*k][1]}) with {L[2*k+1][0]:.3e} ({L[2*k+1][1]}): mismatch x{ratios[k]:.2f}")
pr(f"  Worst mismatch x{max(ratios):.2f}. Gluing needs every pair equal; sides differing up to x{max(ratios):.1f} cannot be glued. No closed solid.")
pr("  Also: a bipyramid's faces meet three or four at each corner, so its faces must share edge lengths in a fixed pattern; sizes here span")
pr(f"  x{L[-1][0]/L[0][0]:.0f} from shortest to longest side, so even faces scaled freely could not keep equal areas and close.")
pr("\nWhat the six faces DO show: they come in three near-mirror pairs, row view vs column view (up-type side vs down-type side):")
for a, b in (("col12", "row12"), ("col13", "row13"), ("col23", "row23")):
    pr(f"  {a} {', '.join(f'{x:.3e}' for x in tris[a])}  vs  {b} {', '.join(f'{x:.3e}' for x in tris[b])}")
pr(f"  Pairs match because |Vus|~|Vcd| ({abs(V[0,1]):.4f} vs {abs(V[1,0]):.4f}) and |Vcb|~|Vts| ({abs(V[1,2]):.4f} vs {abs(V[2,1]):.4f});")
pr(f"  they differ in the shortest side because |Vub| != |Vtd| ({abs(V[0,2]):.5f} vs {abs(V[2,0]):.5f}, x{abs(V[2,0])/abs(V[0,2]):.2f}). Two copies glued")
pr("  edge to edge would make a flat 'pillow' with no inside (zero volume), not a containing solid.")

pr("\nOne partner for both the fold and the leak? (rough estimate, simplest form)")
pr("  A heavy partner D coupled with strength g along u gives a fold k*p0 = g^2/M_D and a mixing angle theta = g/M_D, so M_D = k*p0/theta^2.")
p0 = 2873.0/3; th2 = (1 - s)*6   # the down quark carries 1/6 of u = (-1,-1,2)/sqrt6
for name, k in (("V3 (down squeezed, up stretched)", 0.0660), ("V1 (down only)", 0.1162)):
    pr(f"  {name}: fold {k:.4f}, leak 0.00164 -> M_D ~ {k*p0/th2/1000:.1f} GeV")
pr("  The LHC excludes heavy down-type partners below about 1.5 TeV (ATLAS and CMS searches for vector-like quarks). A partner at 6-11 GeV")
pr("  is ruled out, so the fold amount and the first-row leak cannot both come from one partner in this simplest form.")
