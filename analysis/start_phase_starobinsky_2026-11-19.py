"""Short start phase (space stiffness, Starobinsky) added to TRPG, with knock-on checks (19 Nov 2026). Exploration only.
Formulas: Starobinsky (1980); Mukhanov-Chibisov (1981); Gorbunov-Panin / Gorbunov-Tokareva (scalaron decay);
Davidson-Ibarra (2002); Espinosa et al. (2015) / Herranen et al. (2014) (Higgs during inflation)."""
import math
pi=math.pi
MP=2.435e18          # reduced Planck mass, GeV
As=2.10e-9           # Planck 2018 ln(1e10 As)=3.044
gstar=106.75
hbar=6.582e-25       # GeV s
print("="*72); print("1. THE START PHASE: one number in, everything else out")
for N in (53,54,55,56):
    M=math.sqrt(24*pi**2*As)/N*MP
    ns=1-2/N; r=12/N**2; run=-2/N**2
    print(f"  N={N}: scalaron mass {M:.2e} GeV | tilt n_s {ns:.4f} | wave ratio r {r:.4f} | running {run:.5f}")
N=54.5; M=math.sqrt(24*pi**2*As)/N*MP; r=12/N**2
H=pi*MP*math.sqrt(As*r/2)
print(f"  central (N~54-55 for this start; Bezrukov-Gorbunov 2012): M = {M:.2e} GeV, expansion rate during start H = {H:.2e} GeV")
print("  Tilt against published constraints (McDonough & Ferreira 2025, Table 1; ACT DR6 2025):")
for lab,c,e in [("Planck",0.9638,0.0040),("Planck + DESI BAO",0.9690,0.0035),("CMB-SPA (Planck+ACT+SPT)",0.9693,0.0029),
                ("ACT DR6 P-ACT-LB",0.9743,0.0034),("CMB-SPA + DESI BAO",0.9737,0.0025)]:
    print(f"     {lab:28s} {c:.4f} +/- {e:.4f}: Starobinsky {1-2/N:.4f} is {(1-2/N-c)/e:+.1f} sigma")
print("     -> passes CMB-only data within ~2 sigma; FAILS (>3 sigma) once DESI BAO is combined. The DESI-driven shift comes from")
print("        a ~2 sigma DESI-CMB disagreement (McDonough & Ferreira 2025). PENDING, UNDER STRAIN.")
print("  r: limit < 0.036 (BICEP/Keck 2021). Simons Observatory ~0.003, LiteBIRD ~0.001 (early 2030s) decide it.")
print("  Flatness: |Omega_k| ~ 1e-5 or less (measured |Omega_k| < 0.002). Same temperature everywhere: explained. No isocurvature, f_NL ~ 0.01.")

print("="*72); print("2. REHEATING: how hot the universe gets after the start")
def Trh(xi):
    G=4*M**3/(192*pi*MP**2)*(1-6*xi)**2   # scalaron -> 4 real Higgs scalars
    return (90/(pi**2*gstar))**0.25*math.sqrt(G*MP), G
for xi in (0.0,0.06,0.10,0.12,0.15):
    T,G=Trh(xi); print(f"  Higgs-curvature coupling xi = {xi:.2f}: T_reheat = {T:.1e} GeV")
T0,_=Trh(0.0)

print("="*72); print("3. KNOCK-ON CHECKS AGAINST TRPG'S OWN PIECES")
# MMG
Mmmg=2.4e8
print(f"  MMG (2.4e8 GeV) vs reheat {T0:.0e}: MADE in the heat. A stable MMG would outweigh dark matter (paper sec. 9) -> must decay.")
for y in (2.6e-16,8.8e-12,1e-10):
    G=y**2*Mmmg/(8*pi); print(f"     decay link {y:.1e}: lifetime {hbar/G:.1e} s")
print("     needs link > ~9e-12 (decays before it dominates near 23 GeV; TRPG mmg_audit). Was optional (cold start); now REQUIRED: +1 number.")
print(f"  Coloured far side D (4.6e13) vs scalaron {M:.1e}: never made. HOLDS (no D relics).")
# Higgs safety during the start
print(f"  Higgs during the start: H = {H:.1e} GeV. With the measured top mass the Higgs self-pull turns negative well below H")
print("     (this is why TRPG's grain rule gives 127.3, not 125.1). Safe only if xi >~ 0.06 (Espinosa et al. 2015; Herranen et al. 2014).")
# Majoron
for f in (1e11,1.5e12):
    print(f"  Majoron settling value f = {f:.1e}: H/2pi = {H/(2*pi):.1e} {'>' if H/(2*pi)>f else '<'} f -> {'start angle is averaged (post-start case): no DM isocurvature, angle input disappears' if H/(2*pi)>f else 'isocurvature cap applies'}")
print("     the majoron lead's r < 1e-11 is REPLACED by r ~ 0.004 (the two cannot both hold).")
print(f"     majoron coldness wanted reheat below ~1.6e9 GeV (paper sec. 9); xi = 0: {T0:.0e} (too hot); xi ~0.12-0.15: {Trh(0.12)[0]:.0e}-{Trh(0.15)[0]:.0e} (fits).")
# Matter excess
v=174.0; m3=0.0502e-9
for name,Mn in (("middle partner 4e8",4e8),("heaviest partner 1.5e12",1.5e12)):
    eps=3/(16*pi)*Mn*m3/v**2; print(f"  Davidson-Ibarra ceiling on asymmetry per decay, {name}: {eps:.1e}")
T12,_=Trh(0.12); br=(1.5e12/M)**2/(1-6*0.12)**2
YB=(28/79)*2*br*1.48e-4*0.75*(T12/M)
print(f"  Non-thermal route (scalaron -> heaviest partner, branching {br:.1e} at xi=0.12), at the ceiling: Y_B = {YB:.1e} vs needed 8.7e-11 -> {YB/8.7e-11:.2f} of it")
print("  Thermal route from the middle partner (ceiling 4e-8, needs about 7e-8 even at full efficiency): FAILS (standard Ibarra bound M > ~1e9 GeV).")
print("     caveats: the ceiling is a best case (strictly for the lightest partner); the middle partners, in the heat at this temperature,")
print("     would erase part of it (not included). So the start phase can at best JUST reach the excess. BORDERLINE, not a result.")
print("  -> matter excess stays a measured start number (it already is in the ledger). Count unchanged; the start phase supplies a candidate source.")

print("="*72); print("4. COUNT")
print("  Seeds of galaxies: 2 inputs (size, tilt) -> 1 input (scalaron mass from the size); tilt becomes a PREDICTION (pending, under strain).")
print("  New: wave ratio r = 0.0040 (+/-0.0003) PREDICTED; flatness and same-temperature EXPLAINED (were 'consistent, not derived').")
print("  Added: MMG decay link required (+1 number); xi >~ 0.06 for Higgs safety, ~0.12-0.15 if the majoron lead is kept (fitted choice 7).")
print("  Removed: majoron start angle becomes an average (no longer a start condition); majoron lead's r < 1e-11 withdrawn.")
