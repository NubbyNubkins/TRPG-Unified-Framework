"""Lemma B check: a global charge X is an exact symmetry of the lattice fermion
measure in every topological sector iff its mixed anomalies with the gauge
groups vanish.  Reason: rotating each left-handed field f by exp(i x_f a)
multiplies the measure by exp(-i a sum_f x_f n_f), n_f = index of the overlap
operator in rep f, and by the lattice index theorem (Hasenfratz, Laliena and
Niedermayer 1998; Luscher 1998) sum_f x_f n_f = n3 A(SU3^2 X) + n2 A(SU2^2 X)
+ n1 A(Y^2 X) in each sector.  Part 1 checks the mechanism numerically on a
2D U(1) lattice (index of charge q in flux m is q m); part 2 evaluates the
three sums for TRPG's 27 with X = B - L and, as a control, X = L alone.
"""
import numpy as np, itertools
from fractions import Fraction as F
s1 = np.array([[0,1],[1,0]],complex); s2 = np.array([[0,-1j],[1j,0]]); s3 = np.diag([1,-1]).astype(complex)

def index_2d(q, m, L=8):
    sites = list(itertools.product(range(L), repeat=2)); idx = {s:i for i,s in enumerate(sites)}
    phi = 2*np.pi*m/L**2
    def link(s, mu):                       # constant flux m, Luscher-style gauge
        x, y = s
        if mu == 0: return np.exp(-1j*phi*y)
        return np.exp(1j*phi*L*x) if y == L-1 else 1.0
    N = 2*len(sites); D = np.zeros((N,N), complex); g = [s1, s2]
    for s in sites:
        i = idx[s]; D[2*i:2*i+2, 2*i:2*i+2] += 2*np.eye(2)
        for mu in range(2):
            sp = list(s); sp[mu] = (sp[mu]+1) % L; j = idx[tuple(sp)]
            u = link(s, mu)**q
            D[2*i:2*i+2, 2*j:2*j+2] += -0.5*(np.eye(2) - g[mu])*u
            D[2*j:2*j+2, 2*i:2*i+2] += -0.5*(np.eye(2) + g[mu])*np.conj(u)
    G5 = np.kron(np.eye(len(sites)), s3)
    H = G5 @ (D - np.eye(N)); H = (H + H.conj().T)/2
    w = np.linalg.eigvalsh(H)
    dim_field = int(round(np.sum(w > 0)))       # rank of P- = (1+sign H)/2
    return dim_field - N//2                    # minus the antifield dimension

out = []; p = out.append
p("Part 1: 2D U(1), 8x8 lattice, index n(q, m) of the left-handed overlap field")
tab = {}
for q in (1, 2, 3):
    row = [index_2d(q, m) for m in (0, 1, 2)]; tab[q] = row
    p(f"   charge {q}: flux 0,1,2 -> index {row}")
lin = all(tab[q][m] == tab[1][1]*q*m for q in (1,2,3) for m in (0,1,2))
p(f"   index = (+/-) q*m in every case: {lin}")
for x in [(2, -1, 0), (3, 0, -1), (1, 1, -1), (1, 1, 1)]:
    phases = [sum(x[k]*tab[k+1][m] for k in range(3)) for m in (0,1,2)]
    mixed = sum(x[k]*(k+1) for k in range(3))
    p(f"   global charges x={x} on fields of charge (1,2,3): mixed anomaly sum x.q = {mixed:+d}; "
      f"measure phase factor per unit angle in flux 0,1,2 = {phases} -> {'exact symmetry' if all(v==0 for v in phases) else 'broken by topology'}")

p("")
p("Part 2: TRPG's 27 per generation (left-handed Weyl), mixed sums with the gauge groups")
#            name     d3 d2  Y        B-L       L
fields = [("Q",3,2,F(1,6),F(1,3),F(0)),("u^c",3,1,F(-2,3),F(-1,3),F(0)),("d^c",3,1,F(1,3),F(-1,3),F(0)),
          ("L",1,2,F(-1,2),F(-1),F(1)),("e^c",1,1,F(1),F(1),F(-1)),("N^c",1,1,F(0),F(1),F(-1)),
          ("D",3,1,F(-1,3),F(1,3),F(0)),("Dbar",3,1,F(1,3),F(-1,3),F(0)),
          ("Lfar",1,2,F(1,2),F(1),F(-1)),("Lfarbar",1,2,F(-1,2),F(-1),F(1)),("S",1,1,F(0),F(0),F(0))]
for lab, k in (("B - L (majoron's symmetry)", 4), ("L alone (control)", 5)):
    a3 = sum(F(1,2)*f[2]*f[k] for f in fields if f[1] == 3)
    a2 = sum(F(1,2)*f[1]*f[k] for f in fields if f[2] == 2)
    a1 = sum(f[1]*f[2]*f[3]**2*f[k] for f in fields)
    p(f"   X = {lab:28s} SU3^2 X = {a3}, SU2^2 X = {a2}, Y^2 X = {a1}  -> "
      f"{'exact on the lattice in every sector' if a3 == a2 == a1 == 0 else 'broken by gauge topology (sphalerons), as in nature'}")
print("\n".join(out))
