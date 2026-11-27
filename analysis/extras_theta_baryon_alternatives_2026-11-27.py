"""Extras for Appendix B: (1) baryon and lepton number on the lattice for
TRPG's 27; (2) whether the recorded alternatives fit Craig's covered families."""
from fractions import Fraction as F
#            name     d3 d2  Y        B        L
fields = [("Q",3,2,F(1,6),F(1,3),0),("u^c",3,1,F(-2,3),F(-1,3),0),("d^c",3,1,F(1,3),F(-1,3),0),
          ("L",1,2,F(-1,2),0,1),("e^c",1,1,F(1),0,-1),("N^c",1,1,F(0),0,-1),
          ("D",3,1,F(-1,3),F(1,3),0),("Dbar",3,1,F(1,3),F(-1,3),0),
          ("Lfar",1,2,F(1,2),0,-1),("Lfarbar",1,2,F(-1,2),0,1),("S",1,1,F(0),0,0)]
G = 3
def mixed(k):
    a3 = G*sum(F(1,2)*f[2]*f[k] for f in fields if f[1] == 3)
    a2 = G*sum(F(1,2)*f[1]*f[k] for f in fields if f[2] == 2)
    return a3, a2
for lab, k in (("B", 4), ("L", 5)):
    a3, a2 = mixed(k)
    # change of charge in one unit of SU(2) topology = 2 * A(SU2^2 X) (index of a doublet is 1 per unit)
    p = f"{lab}: SU3^2 = {a3}, SU2^2 = {a2} -> one unit of weak topology changes {lab} by {2*a2}"
    print(p)
print("=> only Delta B = Delta L = 3 (one unit per generation, B - L kept): no proton decay from the lattice measure,")
print("   matching the framework's 'no proton decay' statement at the renormalisable level.")
print()
# Craig's covered families (Table 1): SU(N) with antisym2 + (N-4) antifund (N>=5); SU(N) with sym2 + (N+4) antifund (N>=3); SM generations
print("Recorded alternatives vs Craig's covered families:")
print(" main line   : SM gauge group, 3 x (10 + 5bar) + vector-like pairs + singlets -> covered (Craig + Lemma A)")
print(" gauged B-L  : extra U(1)_{B-L}; the 16 per generation is the Spin(10) spinor -> not in Craig's list;")
print("               Kikukawa's Spin(10) 16 measure: smoothness and locality open (as cited by Craig)")
print(" 4x4 (Pati-Salam SU(4)xSU(2)xSU(2)): (4,2,1)+(4bar,1,2) is not of Craig's form -> not covered;")
print("               it also embeds in Spin(10), same open status")
