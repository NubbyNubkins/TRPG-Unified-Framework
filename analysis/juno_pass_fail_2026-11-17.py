# Pass/fail numbers for the neutrino midpoint hint m2/m3 = 3 - 2*sqrt(2)
import math
r=3-2*math.sqrt(2); m1=0.36e-3
print(f"rule: m2/m3 = {r:.5f};  lightest weight from rule 5: m1 = {m1*1e3} meV")
def ratio(d21,d31,m1): return math.sqrt((m1*m1+d21)/(m1*m1+d31))
def pred21(d31,m1): return r*r*(m1*m1+d31)-m1*m1
data={"JUNO first 59 days alone (Capozzi et al. 2025 reading)":(7.500e-5,0.118e-5),
      "global fit with JUNO (Capozzi et al. 2025)":(7.48e-5,0.10e-5),
      "global fit with JUNO (Esteban et al. 2026)":(7.530e-5,0.097e-5)}
for d31,e31,lab in ((2.519e-3,0.019e-3,"with SK atmospheric"),(2.529e-3,0.021e-3,"without SK atmospheric")):
    p=pred21(d31,m1); ep=r*r*e31
    print(f"\nlarge gap {d31*1e3:.3f} +- {e31*1e3:.3f} e-3 eV2 ({lab}):")
    print(f"  predicted small gap = {p*1e5:.3f} +- {ep*1e5:.3f} e-5 eV2   (with m1 = 0: {pred21(d31,0)*1e5:.3f})")
    for k,(d,e) in data.items():
        rr=ratio(d,d31,m1); er=rr*0.5*math.hypot(e/d,e31/d31)
        print(f"  {k}: small gap {d*1e5:.3f}+-{e*1e5:.3f}; m2/m3 = {rr:.4f} +- {er:.4f}; off by {(rr-r)/er:+.1f} sigma")
d31=2.519e-3; p=pred21(d31,m1)
print("\nJUNO design goal (JUNO 2022): small gap to 0.3%, large gap to 0.2% after 6 years")
e=math.hypot(0.003*p, r*r*0.002*d31)
print(f"  then the test has width {e*1e5:.3f} e-5: PASS band (2 sigma) {1e5*(p-2*e):.3f} to {1e5*(p+2*e):.3f}; FAIL (3 sigma) outside {1e5*(p-3*e):.3f} to {1e5*(p+3*e):.3f}")
for c in (7.48e-5,7.50e-5,7.53e-5):
    print(f"  if the small gap stays at {c*1e5:.2f}: {abs(c-p)/e:.1f} sigma away")
print(f"  the rule in one line: small gap / large gap = {p/d31:.5f}  (measured now {7.53e-5/2.519e-3:.5f} +- {7.53e-5/2.519e-3*math.hypot(0.097/7.53,0.019/2.519):.5f})")
