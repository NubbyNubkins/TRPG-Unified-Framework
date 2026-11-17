# The two partials on the geometry
import math
L10=math.log10
print("A) Strengths: the leak as a reflection on the centre line (log weight)")
grain=8.63e18; F=4.56e13
pair=F*F/grain
print(f"   grain {grain:.3g}, meeting F {F:.3g} -> reflection of grain through F = {pair:.3g} GeV (weak pairs 2.2-2.5e8 needed)")
print(f"   steps of ten: grain-F = {L10(grain/F):.2f}, F-pair = {L10(F/pair):.2f} (equal by construction)")
for nm,mD,mnu in (("nu3 / tau link",1.77686,50.2e-12),("nu2 / muon link",0.10566,8.66e-12)):
    M=mD*mD/mnu; print(f"   neutrino seesaw same shape: partner {M:.2g} GeV reflected through {mD} GeV -> {mnu*1e9*1e3:.3g} meV ; steps {L10(M/mD):.2f} = {L10(mD/mnu):.2f}")
print("   => leak and seesaw are one operation: a colourless thing with no partner on its own half is the mirror image of a heavy partner")
print("\nB) Quark masses: cube angles vs ladder openings (opening = angle of the ladder's cone)")
cube={"face diagonal to body diagonal  35.26":math.degrees(math.atan(1/math.sqrt(2))),
      "45":45.0,"edge to body diagonal  54.74":math.degrees(math.acos(1/math.sqrt(3))),"60":60.0}
op={"neutrinos (m1 0.36 meV)":37.13,"charged leptons":45.0,"down":48.15,"up":52.147}
for k,v in op.items():
    best=min(cube.items(),key=lambda kv:abs(kv[1]-v)); print(f"   {k:24s} {v:7.3f}  nearest cube angle {best[1]:.2f} ({v-best[1]:+.2f} deg)")
def kq(m):
    Q=sum(m)/sum(x**.5 for x in m)**2; e=math.sqrt(2*(3*Q-1)); return math.degrees(math.atan(e/math.sqrt(2)))
d21,d31=7.42e-5,2.511e-3
print("   neutrino opening vs lightest mass m1 (normal order):")
for m1 in (0,0.1e-3,0.36e-3,1e-3,3e-3,10e-3):
    print(f"     m1 {m1*1e3:5.2f} meV -> {kq([m1,(m1*m1+d21)**.5,(m1*m1+d31)**.5]):.2f} deg")
lo,hi=0.36e-3,3e-3
for _ in range(60):
    mid=(lo+hi)/2; f=kq([mid,(mid*mid+d21)**.5,(mid*mid+d31)**.5])
    lo,hi=(mid,hi) if f>35.2644 else (lo,mid)
print(f"   35.26 deg needs m1 = {lo*1e3:.2f} meV; rule 5 gives 0.36 meV -> cube angle contradicts rule 5: FAIL in-system")
print("   up 52.15 vs 54.74: -2.6 deg; up opening uncertainty ~0.3 deg from m_u, m_c -> FAIL")
print("\n   Press direction on the cube: 89% along the body diagonal (1,1,1) = the line through both overlap points;")
print("   10% along (1,1,-2) = toward/away from one family corner (octahedron direction); 1% (1,-1,0) quarter turn")
print("\nC) Rule 7's 45 deg shared between quark and lepton tilts (prior art: quark-lepton complementarity, Raidal 2004; Minakata & Smirnov 2004)")
thC=13.04; s12=0.307; e=0.0125; th12=math.degrees(math.asin(s12**.5)); dth=math.degrees(e/(2*math.sqrt(s12*(1-s12))))
S=thC+th12; print(f"   Cabibbo {thC} + solar {th12:.2f}+-{dth:.2f} = {S:.2f} deg ; 45 is {(S-45)/dth:.1f} sigma away -> near-miss FAIL")
