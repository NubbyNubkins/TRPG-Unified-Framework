# Le Verrier test - which single missing particle would close the force strengths (one loop, PDG 2024 inputs)
import math
MZ=91.1876; a_em=1/127.951; s2=0.23122; a3=0.1180
a1=(5/3)*a_em/(1-s2); a2=a_em/s2
A=[1/a1,1/a2,1/a3]; b=[41/10,-19/6,-7]
def run(mu,extra):
    t=math.log(mu/MZ); out=[]
    for i in range(3):
        v=A[i]-b[i]/(2*math.pi)*t
        for db,M in extra:
            if mu>M: v-=db[i]/(2*math.pi)*math.log(mu/M)
        out.append(v)
    return out
def mismatch(extra):
    # scale where 1 and 2 meet, then how far 3 is off
    lo,hi=math.log(1e3),math.log(1e19)
    for _ in range(200):
        m=(lo+hi)/2; v=run(math.exp(m),extra)
        lo,hi=(m,hi) if v[0]>v[1] else (lo,m)
    mu=math.exp(lo); v=run(mu,extra); return mu,(v[2]-v[0])/v[0]
mu,d=mismatch([]); print(f"SM alone: 1&2 meet at {mu:.2e} GeV, strong line off by {100*d:+.1f}% (one loop; two loop gives -8.9%)")
cands={ # (dB1,dB2,dB3) for the whole multiplet, GUT-normalised hypercharge
 "far-side weak pair (vector-like lepton doublet, Dirac)":(0.4,2/3,0),
 "weak triplet, no colour (Majorana fermion)":(0,4/3,0),
 "colour octet, no weak charge (Majorana fermion)":(0,0,2),
 "extra Higgs-like doublet (scalar)":(0.1,1/6,0),
 "coloured far-side part D (vector-like down quark)":(4/15,0,2/3),
 "D and weak pair together, same mass (whole 5+5bar)":(4/15+0.4,2/3,2/3),
}
print("\nOne missing particle at mass M: which M closes the meeting?")
for nm,db in cands.items():
    best=None
    for k in range(0,1700):
        M=10**(2+k*0.01); mu,d=mismatch([(db,M)])
        if best is None or abs(d)<abs(best[2]): best=(M,mu,d)
    M,mu,d=best
    ok=abs(d)<0.002
    print(f"  {nm:55s}: {'closes' if ok else 'cannot close'}  M={M:.2e} GeV, meeting {mu:.2e} GeV, left {100*d:+.2f}%")
print("\nTwo-loop shifts these masses upward roughly tenfold (weak pair: 2.6e7 one loop -> 2.2e8 two loop, an earlier test).")
print("Proton safety for a true meeting of forces needs the meeting above ~1e15.5 GeV (Super-K).")

print("\nLe Verrier step 2: the far-side weak pairs come in three placements (one per family)")
for n in (2,3):
    db=tuple(n*x for x in (0.4,2/3,0)); best=None
    for k in range(0,1700):
        M=10**(2+k*0.01); mu,d=mismatch([(db,M)])
        if best is None or abs(d)<abs(best[2]): best=(M,mu,d)
    print(f"  {n} weak pairs at one mass: M={best[0]:.2e} GeV, meeting {best[1]:.2e}, left {100*best[2]:+.2f}%")
print("\nLe Verrier step 3: what the leftover demands of ANY single missing body")
print("  closing needs the weak line bent more than hypercharge and colour: (dB2 - dB1) large, dB3 small")
for nm,db in cands.items():
    print(f"  {nm:55s}: weak-minus-hypercharge {db[1]-db[0]:+.2f}, colour-minus-weak {db[2]-db[1]:+.2f}")
print("  => the missing body for the strengths carries WEAK charge and NO colour")
print("  => the press needs COLOUR and a link to the down quarks (the D, already in the matrix)")
print("  => one body cannot do both jobs: the leftovers point at two different kinds")
