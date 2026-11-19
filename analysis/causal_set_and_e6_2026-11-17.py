"""Borrowed tests 3 and 4.
3. Causal sets, 'everpresent Lambda' (Sorkin 1987/1991; Ahmed, Dodelson, Greene & Sorkin 2004):
   dark energy ~ +- 1/sqrt(V), V = four-volume of our past (number of space-time points N = V / l_P^4, fluctuation sqrt N).
   TEST FIXED BEFORE RUNNING: compute 1/sqrt(V) for our past light cone in the measured expansion (Planck 2018: H0 67.4,
   Omega_m 0.315, Omega_r 9.1e-5) and compare with the measured Lambda = 3 Omega_L H0^2 / c^2.
   It REPLACES TRPG's dark-energy input only if the coefficient is fixed by the theory (it is not: the published model carries
   a free order-one amplitude and a random sign). Otherwise it is at most an explanation of the SIZE.
4. E6 tables (Guersey, Ramond & Sikivie 1976; Slansky 1981): exact charges of the 27, anomaly sums, and whether any
   renormalisable interaction can break baryon number with TRPG's content (12 known carriers, one Higgs pair, the majoron
   singlet, no coloured scalars). Then the proton lifetime from the grain, compared with Super-Kamiokande (2.4e34 years).
"""
import math, numpy as np
from fractions import Fraction as Fr
from scipy.integrate import quad
print("3. CAUSAL-SET DARK ENERGY")
c=2.99792458e8; Mpc=3.0857e22; H0=67.4e3/Mpc; Om=0.315; Or=9.1e-5; OL=1-Om-Or
H=lambda a: H0*math.sqrt(Or/a**4+Om/a**3+OL)
chi=lambda a: quad(lambda x: c/(x*x*H(x)),a,1,limit=200)[0]          # comoving distance from scale a to today
V=quad(lambda a: (4*math.pi/3)*chi(a)**3*a**3/(a*H(a))*c,1e-8,1,limit=200)[0]   # four-volume of past cone, m^4 (c dt)
Lam_meas=3*OL*H0**2/c**2
print(f"  four-volume of our past: {V:.3e} m^4;  1/sqrt(V) = {1/math.sqrt(V):.3e} m^-2;  measured Lambda = {Lam_meas:.3e} m^-2")
print(f"  amplitude the model must be given: Lambda * sqrt(V) = {Lam_meas*math.sqrt(V):.2f}  (order one: the SIZE is explained)")
print("  the coefficient and the sign are not fixed by the theory -> input count unchanged.")
print("  the value wanders with time in this model; TRPG's present reading (the kept amount of the bag) is a constant.")
print()
print("4. E6 TABLES")
F=lambda x: Fr(x)
# left-handed Weyl fields of one 27: (name, colour dim, colour sign (+1 = 3, -1 = 3bar, 0 = singlet), weak dim, Y, baryon number)
st=[("Q",3,1,2,F('1/6'),F('1/3')),("u^c",3,-1,1,F('-2/3'),F('-1/3')),("d^c",3,-1,1,F('1/3'),F('-1/3')),
    ("L",1,0,2,F('-1/2'),0),("e^c",1,0,1,F(1),0),("nu^c",1,0,1,F(0),0),
    ("D",3,1,1,F('-1/3'),F('1/3')),("D^c",3,-1,1,F('1/3'),F('-1/3')),("L'",1,0,2,F('-1/2'),0),("L'^c",1,0,2,F('1/2'),0),("S",1,0,1,F(0),0)]
n=sum(cd*wd for _,cd,_,wd,_,_ in st); print(f"  states: {n} (16 + 10 + 1)")
for name,cd,cs,wd,Y,B in st:
    Qs=[Y+Fr(1,2),Y-Fr(1,2)] if wd==2 else [Y]
    print(f"   {name:5s} colour {cd} weak {wd} Y {str(Y):5s} charges {[str(q) for q in Qs]}")
T={1:0,3:Fr(1,2)}; T2={1:0,2:Fr(1,2)}
an={"SU(3)^2 U(1)":sum(T[cd]*wd*Y for _,cd,_,wd,Y,_ in st),
    "SU(2)^2 U(1)":sum(T2[wd]*cd*Y for _,cd,_,wd,Y,_ in st),
    "U(1)^3":sum(cd*wd*Y**3 for _,cd,_,wd,Y,_ in st),
    "gravity^2 U(1)":sum(cd*wd*Y for _,cd,_,wd,Y,_ in st),
    "SU(3)^3":sum(cs*wd for _,cd,cs,wd,_,_ in st if cd==3)}
for k,v in an.items(): print(f"  anomaly {k:16s} = {v}")
nd=sum(cd for _,cd,_,wd,_,_ in st if wd==2); print(f"  weak doublets per family: {nd} ({'even: Witten anomaly absent' if nd%2==0 else 'ODD'})")
# baryon number: every renormalisable term with fermions is (fermion)(fermion)(scalar) or a fermion mass term.
# With only colourless scalars, the two fermions must form a colour singlet: 3 x 3bar or 1 x 1. Check B of every such allowed pair.
viol=[]
fields=[(nm,cd,cs,B) for nm,cd,cs,wd,Y,B in st]
for i,a in enumerate(fields):
    for b in fields[i:]:
        singlet=(a[2]==0 and b[2]==0) or (a[2]*b[2]==-1)
        if singlet and a[3]+b[3]!=0: viol.append((a[0],b[0]))
print(f"  colour-singlet fermion pairs that change baryon number: {viol if viol else 'none'}")
print("  -> with no coloured scalars, no renormalisable term breaks baryon number (three quarks are needed, i.e. dimension 6).")
mp=0.938272; hbar=6.582e-25; yr=3.156e7
tau=lambda M,a: M**4/(a*a*mp**5)*hbar/yr
print(f"  proton lifetime if dimension-6 decay comes only from the grain (8.63e18 GeV, strength 1): {tau(8.63e18,1):.1e} years")
print(f"  for comparison, a true one-force merger at F = 4.6e13 GeV (alpha ~ 1/40): {tau(4.6e13,1/40):.1e} years  (Super-K floor 2.4e34: FAIL)")
print("  risk: if the Higgs had to sit in a whole 27, its coloured parts would be scalar D's that link quarks to leptons and quark pairs -> decay returns.")
