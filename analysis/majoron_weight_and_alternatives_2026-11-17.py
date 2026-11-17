# (A) dark matter with no particle at all: black holes from the start? (B) one setting for both far-side jobs; (C) a weight rule for the majoron
import math
grain=8.633e18; MMG=2.41e8; F=math.sqrt(MMG*grain); v=174.0
print("(A) Black holes made at the start as the dark matter (Hawking 1971; Carr & Hawking 1974): needs no particle")
Teq=0.8e-9  # GeV
for T in (1e6,1e7,1e8,1e9):
    M_g=2e38*(1e-3/T)**2      # horizon mass in grams: 1e5 suns at 1 MeV
    beta=Teq/T                # share that must collapse for all the dark matter
    s_rad=0.45/math.sqrt(2)/ (math.sqrt(math.log(1/beta))) # erfc ~ exp(-x^2): rough
    s_mat=(beta/0.056)**0.2   # matter-like era before the light (Khlopov & Polnarev 1980; Harada et al. 2016)
    print(f"   start level {T:.0e} GeV: hole mass {M_g:.0e} g; ripple needed {s_mat:.1e} (held-energy era) or {s_rad:.2f} (light era); measured on large sizes 4.6e-5: short by {s_mat/4.6e-5:.0f}x")
print("   open window for such holes as all the dark matter: 1e17 to 1e22 g (Carr, Kohri, Sendouda & Yokoyama 2021)")
print("\n(B) One setting for both far-side jobs (D's weight F and the hidden partners' weights)")
print(f"   F = {F:.2e} GeV; a link cannot exceed about 3.5 -> settling value f >= {F*math.sqrt(2)/3.5:.2e} GeV (the examples in an earlier test used 6.3e10)")
K=8.8e3*(6.3e10)**4  # eV GeV^4 : m_J * angle^4 * f^4
for nm,f in (("largest link 3.5",F*math.sqrt(2)/3.5),("link root 2 (f = F)",F),("link 1",F*math.sqrt(2))):
    print(f"   {nm:20s}: f = {f:.2e}; majoron weight x angle^4 = {K/f**4:.1e} eV; heaviest partner's link {6.3e10*math.sqrt(2)/f:.1e}; jitter limit on the start's stretching {1.8e6*f/6.3e10:.1e} GeV x angle")
print("   with angles 0.3 to 3: weight between about 4e-10 and 2e-4 eV (was: far below 1 eV up to 350 keV)")
print("\n(C) Weight rule tried: the tilt comes through the grain like the MMG's weight: majoron weight = (Higgs value)^2 / grain")
for nm,vv in (("174 GeV",174.0),("246 GeV",246.0)):
    m=vv*vv/grain*1e9
    for fn,f in (("largest link",F*math.sqrt(2)/3.5),("f = F",F)):
        print(f"   Higgs value {nm}: weight {m:.1e} eV; {fn}: starting angle needed {(K/f**4/m)**0.25:.2f} (random average 1.81; largest 3.14)")
print("   luck: the allowed weights span about 5 tens, and 8 simple balance points were already tried in an earlier test: a landing here is expected by chance. HINT only.")
print("\n(D) CORRECTION to (B): the settling piece carries matter number 2 and a D pair carries none, so 'link x f' is not allowed for D.")
print("    The allowed form is F = f^2 / (a larger energy). Largest energy available: the grain. Smallest: f itself.")
f_hi=math.sqrt(F*grain)
for nm,f in (("f = F (smallest)",F),("f = root(F x grain): everything through the grain",f_hi)):
    print(f"   {nm:48s}: f = {f:.2e} GeV; majoron weight x angle^4 = {K/f**4:.1e} eV")
print(f"   chain through the grain: f = {f_hi:.2e} -> F = f^2/grain = {f_hi**2/grain:.2e} -> MMG = F^2/grain = {F*F/grain:.2e} GeV")
for lim,nm in ((2e-20,"Rogers & Peiris 2021"),(3e-19,"Dalal & Kravtsov 2022, as recalled")):
    print(f"   wave-too-light limit {lim:.0e} eV ({nm}): starting angle must be below {(K/f_hi**4/lim)**0.25:.2f}")
print(f"   heaviest hidden partner's link: {6.3e10*math.sqrt(2)/f_hi:.1e}; jitter limit on the start's stretching: {1.8e6*f_hi/6.3e10:.1e} GeV x angle")
