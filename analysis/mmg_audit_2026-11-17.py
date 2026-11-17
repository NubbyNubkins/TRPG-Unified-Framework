# Audit of the MMG. Is it needed? Is it a problem? Can it be measured? Does anything need another particle?
import math
M=2.41e8; grain=8.633e18; F=math.sqrt(M*grain); v=174.0; MPl=1.22e19; g=106.75
HoverT2=1.66*math.sqrt(g)/MPl
print("1. NEEDED? (numbers from earlier tests, two loops)")
print("   strengths without it: miss -8.9% +- 0.1% (60 sigma). With three far-side weak pairs at 2.23e8 GeV (1.8e8 to 2.7e8): meet at 4.6e13 GeV.")
print("   nothing else in the record uses the MMG: not the weights, the tilt, the twist, the matter excess or the dark matter.")
print("   it adds no state: it is the colourless part of the far side already in the matrix; without the split it would sit at the meeting with D.")
print(f"   MMG = F^2/grain with F = {F:.2e} GeV: the same statement as the gravity leak, not a second assumption.\n")
print("2. A PROBLEM? If the MMG has no link at all it is stable, and any made in the heat are still here.")
print("   made while the light was building up (Chung, Kolb & Riotto 1998; Giudice, Kolb & Riotto 2001): amount/dark matter ~ (M^2 sigma v / 0.12) x (2000 T_start / M)^7")
s=0.01
for T in (1.7e5,1e6,5.5e7,2.4e8):
    print(f"   start level {T:.1e} GeV: {s*(2000*T/M)**7/0.12:.1e} times the dark matter")
Tmax=M/2000*(0.12/s)**(1/7)
print(f"   a stable MMG needs the start level below {Tmax:.1e} GeV")
X=1.68e-3; M2=1.63e10
print(f"   but then the handover number is at most {Tmax/(2*M2):.1e} against {X:.2e} needed: short by {X/(Tmax/(2*M2)):.0f}x")
print("   => a stable MMG and the matter-excess lead cannot both stand. The MMG must have a link that lets it decay.")
lo_bbn=math.sqrt(6.6e-25*8*math.pi/M)
rho_s=3.6e-9*0.1*(M/1.1e3)**2        # GeV, if it reached thermal amounts (Higgsino-like scaling)
Tdom=4/3*rho_s; lo_dil=math.sqrt(HoverT2*Tdom**2*8*math.pi/M)
hi_wash=math.sqrt(HoverT2*M*8*math.pi)
print(f"   link window: above {lo_bbn:.1e} (decays within a second); above {lo_dil:.1e} (decays before it would dominate at {Tdom:.0f} GeV and thin out the matter excess);")
print(f"                below {hi_wash:.1e} (does not erase the matter excess); below 6e-4 (neutrino weights, an earlier test)")
for nm,e in (("MMG weight / grain",M/grain),("F / grain",F/grain),("Higgs value / grain",v/grain)):
    ok = lo_dil<e<hi_wash
    print(f"   natural size {nm:20s} = {e:.1e}: {'inside the window' if ok else 'OUTSIDE (too small)' if e<lo_dil else 'outside (too big)'};"
          f" lifetime {6.58e-25/(e*e*M/(8*math.pi)):.1e} s; neutrino weight it adds {(e*v)**2/M*1e9:.1e} eV")
print("\n3. MEASURABLE?")
th=0.0003   # remaining theory error after the three-loop running (leak_establish_three_loop_2026-11-16): top threshold
print("   direct: no (about 15,000 times beyond the LHC). Indirect: the strong coupling, three-loop running done: 0.1178 (grain = Planck/root 2) or 0.1196 (grain = Planck)")
for nm,a,e in (("ALPHA 2025",0.11873,0.00056),("PDG average",0.1180,0.0009)):
    print(f"   {nm}: {a} +- {e}: {(a-0.1178)/math.hypot(e,th):+.1f} sigma from 0.1178, {(a-0.1196)/math.hypot(e,th):+.1f} sigma from 0.1196")
print(f"   grain factor that would fit ALPHA exactly: about {1.414-(0.11873-0.1178)/(0.1196-0.1178)*0.414:.2f} (no reason behind it)")
print(f"   at 0.3% (+-0.00035): the two grain choices are {0.0018/math.hypot(0.00035,th):.1f} sigma apart")
print("   also: no proton decay at any rate (Hyper-Kamiokande seeing one would break it); no primordial gravity waves.")
