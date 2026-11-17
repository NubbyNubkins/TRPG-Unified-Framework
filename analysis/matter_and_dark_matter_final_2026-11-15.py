"""Open items 4 and 5 together (15 Nov 2026): the matter excess and dark matter from ONE settling event.
The lepton-number-breaking piece sigma = (f + rho) e^{iJ/f} / sqrt2 settles at f. It gives the hidden partners their masses
(M_i = y_i f / sqrt2), leaves the majoron J (mass m_J from a tiny breaking), and its heavy radial part is rho ("the settling piece").
Asked: is there a region where eta_B = 6.1e-10 AND Omega_J h^2 = 0.120 together, how many free numbers, and does TRPG fix any?
Deterministic (fixed seeds). Runtime about 1 minute.  Usage: python3 matter_and_dark_matter_final_2026-11-15.py > ...txt
Credit: Chikashige, Mohapatra & Peccei 1981 (majoron); Rothstein, Babu & Seckel 1993; Akhmedov, Berezhiani, Mohapatra & Senjanovic 1993;
Berezinsky & Valle 1993; Lattanzi & Valle 2007; Garcia-Cely & Heeck 2017; Heeck & Patel 2019; Preskill, Wise & Wilczek 1983; Abbott & Sikivie 1983;
Dine & Fischler 1983; Turner 1986; Arias et al. 2012 (misalignment); Linde 1985, Lyth 1990, Planck 2018 X (isocurvature); Sikivie 1982 (domain walls);
Lazarides & Shafi 1991; Asaka, Hamaguchi, Kawasaki & Yanagida 1999; Fukugita & Yanagida 1986; Covi, Roulet & Vissani 1996; Davidson & Ibarra 2002;
Abada et al. 2006; Nardi et al. 2006; Buchmuller, Di Bari & Plumacher 2005; Harvey & Turner 1990; Hall, Jedamzik, March-Russell & West 2010;
Frigerio, Hambye & Masso 2011; Heeck & Teresi 2017; Kamada & Yanagi 2019; Boucenna, Morisi, Shafi & Valle 2014; Irsic et al. 2017, 2024;
Audren et al. 2014; Poulin, Serpico & Lesgourgues 2016; Nygaard, Tram & Hannestad 2021; Simon et al. 2022; Albright & Rodejohann 2009 (TM1);
NuFIT 6.0 (2024); DESI 2024, 2025; KamLAND-Zen 2024; Planck 2018; PDG 2024."""
import math, os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*110); pr(t); pr("="*110)

# ---------------------------------------------------------------- constants (GeV unless said)
MPL, MRED = 1.220890e19, 2.435323e18                 # Planck mass, reduced Planck mass
V174, V246, MH = 174.104, 246.22, 125.20            # Higgs value, Higgs mass (PDG 2024)
ME, MMU, MTAU, MP = 0.51099895e-3, 0.1056583755, 1.77693, 0.93827209
ALPHA = 1/137.035999
S0, RHOC = 2891.2, 1.053672e-5                       # entropy today (cm^-3), critical density / h^2 (GeV cm^-3)
OM_DM, OM_DM_ERR = 0.1200, 0.0012                    # Planck 2018
OM_B = 0.02237
ETA_B, ETA_B_ERR = 6.12e-10, 0.04e-10                # Planck 2018 (baryon-to-photon ratio)
YB = ETA_B/7.04                                      # baryons per entropy
RHO_OVER_S_DM = OM_DM*RHOC/S0                        # GeV: dark matter mass density per entropy (0.437 eV)
HBAR = 6.582119569e-25                               # GeV s
GYR = 3.15576e16
MNU = np.array([0.00036, 0.00866, 0.05034])*1e-9     # TRPG neutrino masses, sum 0.0594 eV
def gstar(T):                                        # rough table of the number of light species (g* ~ g*s above 1 MeV)
    tT = np.log10([1e-7, 1e-4, 1e-3, 0.1, 0.2, 0.4, 1, 5, 80, 300, 1e19]); tg = [3.36, 3.36, 10.75, 17.0, 30.0, 60.0, 73.0, 86.0, 100.0, 106.75, 106.75]
    return float(np.interp(math.log10(T), tT, tg))
def gstar_s(T): return 3.91 if T < 2e-4 else gstar(T)

# ================================================================= 1. misalignment
head("1. Dark matter as the majoron's slow settling wave (misalignment). The exact oscillator, not a rule of thumb")
# phi'' + 3H phi' + m^2 phi = 0 with H = 1/(2t). In units m = 1. Late-time rho a^3 compared with (1/2) m^2 phi0^2 (H/m)^(3/2).
def ode(t, y): return [y[1], -1.5/t*y[1] - y[0]]
sol = solve_ivp(ode, [1e-4, 400.0], [1.0, 0.0], rtol=1e-10, atol=1e-13, dense_output=True)
ts = np.linspace(300, 400, 4001); ph, dph = sol.sol(ts)
rho = 0.5*(dph**2 + ph**2); c32_num = float(np.mean(rho/(0.5*(1/(2*ts))**1.5)))
c32_exact = math.gamma(1.25)**2*2**1.5/math.pi*2**1.5
C_OSC = c32_exact**(2/3)
pr(f"  late density = (1/2) m^2 phi0^2 x K x (H/m)^(3/2):  K numeric {c32_num:.4f}, exact (Bessel) {c32_exact:.4f}")
pr(f"  -> the swing effectively starts when H = m/{C_OSC:.3f}.  The rule of thumb '3H = m' (used in an earlier test) overstates the amount by {(3/C_OSC)**1.5:.2f}x,")
pr(f"     so every mass quoted from it must go up by {(3/C_OSC)**3:.2f}x.")
def omega_mis(mJ_eV, f, theta=1.0, n=1):
    """Omega h^2 of a majoron wave started at angle theta (in units of its own period 2 pi f/n). Harmonic. T_R above the onset."""
    m = mJ_eV*1e-9; T = 1e3
    for _ in range(40): T = math.sqrt(m*MPL/(C_OSC*1.6602*math.sqrt(gstar(T))))
    s = 2*math.pi**2/45*gstar_s(T)*T**3
    return 0.5*m*m*(theta*f/n)**2/s*S0/RHOC, T
def mJ_for_dm(f, theta=1.0, n=1):
    lo, hi = -30.0, 12.0
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if omega_mis(10**mid, f, theta, n)[0] < OM_DM: lo = mid
        else: hi = mid
    return 10**(0.5*(lo+hi))
F0 = 6.3e10
m0 = mJ_for_dm(F0); o14, Tosc14 = omega_mis(1400.0, F0)
pr(f"  f = 6.3e10 GeV, starting angle 1: all of the dark matter needs m_J = {m0/1e3:.2f} keV (an earlier test said 1.4 keV; with 1.4 keV the wave is only {o14/OM_DM*100:.0f}% of it).")
pr(f"  The swing starts at T = {omega_mis(m0, F0)[1]:.2e} GeV. Scaling: m_J = {m0/1e3:.1f} keV x (6.3e10/f)^4 / angle^4 (one relation among three numbers).")
pr("   f (GeV)     m_J for angle 0.3      angle 1        angle 1.81 (random average)   angle 2.5")
for f in (2.5e10, 6.3e10, 1.0e11, 1.4e11, 1.0e12, 1.4e12, 3e13, 1e15):
    pr(f"   {f:8.1e}   " + "   ".join(f"{mJ_for_dm(f, th):12.3e} eV" for th in (0.3, 1.0, 1.81, 2.5)))

# ================================================================= 2. hot majorons
head("2. Hot majorons: when the heat itself makes them (this is what kills the hot start)")
Yeq = (1.2020569/math.pi**2)/(2*math.pi**2/45*106.75)
m_th = RHO_OVER_S_DM/Yeq*1e9
pr(f"  A majoron once in contact with the heat keeps {Yeq:.2e} per unit entropy. All of the dark matter then needs m_J = {m_th:.0f} eV,")
pr(f"  and it is a thermal relic of that mass. Lyman-alpha needs a thermal relic above 5.3 keV (Irsic et al. 2017, 95%) or 5.7 keV (Irsic et al. 2024):")
pr(f"  {m_th:.0f} eV is {5300/m_th:.0f}x too light. FAIL. At m_J = 8.6 keV thermal majorons would be {8600*1e-9*Yeq/RHO_OVER_S_DM:.0f}x too much dark matter. FAIL either way.")
pr(f"  With the earlier numbers (130 eV): hot share {130e-9*Yeq/RHO_OVER_S_DM*100:.0f}% (an earlier test quoted 73%; same statement).")
def Y_hot(f, Tmax, Ms, csig=1.0):
    """Majorons made by pairs of heavy partners in the heat (N N -> J J), cross-section ~ csig M^2/(16 pi f^4): order-of-magnitude freeze-in."""
    Y = 0.0
    for M in Ms:
        sv = csig*M*M/(16*math.pi*f**4); nN = 0.75*1.2020569/math.pi**2*2
        y = 2*nN**2*sv*MPL*Tmax/(2*math.pi**2/45*106.75*1.6602*math.sqrt(106.75))
        if Tmax < M: y *= math.exp(-2*M/Tmax)*(M/Tmax)**3
        Y += y
    return min(Y, Yeq)
def f_min_cold(Tmax, Ms, mJ_of_f, hot_frac=0.01):
    lo, hi = 9.0, 17.0
    for _ in range(60):
        mid = 0.5*(lo+hi); f = 10**mid
        if Y_hot(f, Tmax, Ms)*mJ_of_f(f)*1e-9 > hot_frac*RHO_OVER_S_DM or Tmax > f: lo = mid
        else: hi = mid
    return 10**hi
pr("  How large must f be so that a heat of temperature T_max makes under 1% hot majorons (the wave still being all of the dark matter, angle 1)?")
MS_CORNER = [3.59e4, 1.63e10, 1.00e11]
for Tm in (1e9, 1e10, 1e11, 3e11):
    fm = f_min_cold(Tm, MS_CORNER, lambda f: mJ_for_dm(f))
    pr(f"   T_max = {Tm:.0e} GeV (partners {MS_CORNER[1]:.1e}, {MS_CORNER[2]:.1e}): f > {fm:.1e} GeV, where m_J = {mJ_for_dm(fm):.1e} eV and the heaviest partner's link y_3 = {math.sqrt(2)*MS_CORNER[2]/fm:.1e}")
F_HOT = f_min_cold(1e11, MS_CORNER, lambda f: mJ_for_dm(f))
pr(f"  -> A hot start (T_max ~ 1e11 GeV, what thermal leptogenesis needs) is compatible with majoron dark matter only if f is above about {F_HOT:.0e} GeV,")
pr(f"     where the majoron is a light wave (m_J about {mJ_for_dm(F_HOT):.1f} eV or less). At f = 6.3e10 to 1e11 GeV it is not: hot share at f = 6.3e10, T_max = 1e11: "
   f"{min(Y_hot(F0, 1e11, MS_CORNER), Yeq)*mJ_for_dm(F0)*1e-9/RHO_OVER_S_DM:.0f}x the whole dark matter (must be under a few % of it).")
pr("     (Cross-section known to a factor of a few; f_min moves as its 1/7 power along the dark-matter line.)")
pr("\n  2b. Freeze-in instead of the wave (majorons leaking slowly out of the heat; Hall et al. 2010; Frigerio, Hambye & Masso 2011):")
pr(f"   They are born with momentum about 2.5 to 3.2 x temperature, so Lyman-alpha (thermal relic > 5.3 keV) needs m_J > {7.2*(5.3/3)**(4/3)*2.5/3.15:.0f} to {7.2*(5.3/3)**(4/3)*3.2/3.15:.0f} keV")
pr("   (same range as Heeck & Teresi 2017: 12 to 16 keV).")
mfi = 16e3; Yfi = RHO_OVER_S_DM/(mfi*1e-9)
lam_eff = math.sqrt(Yfi/(2*135/(8*math.pi**3*1.6602*100**1.5)*MPL/MH**2)/(V246**2/(32*math.pi*MH)))
pr(f"   Through the Higgs (h -> J J by mixing with the settling piece): needs an effective link {lam_eff:.1e}, i.e. (Higgs-to-settling link)/(settling self-link) = "
   f"{2*lam_eff*F0**2/MH**2:.0e} at f = 6.3e10.")
pr("   That takes two more numbers, both extremely small, and drags the Higgs mass up unless tuned. Not in the main line: NOT PURSUED (allowed only by adding numbers).")
ffi = (1.0*MS_CORNER[1]**2/(16*math.pi)*2*(0.75*1.2020569/math.pi**2*2)**2*MPL*MS_CORNER[1]/(2*math.pi**2/45*106.75*1.6602*math.sqrt(106.75))/Yfi)**0.25
pr(f"   Through the partners (N N -> J J in a heat at T ~ M_2 = 1.6e10): right amount for f about {ffi:.0e} GeV. But there the wave overfills unless the starting angle is under "
   f"{(mJ_for_dm(ffi)/mfi)**0.25:.3f}, and the matter excess is the thermal one (2.6 to 3x short). FAIL as a joint solution.")

# ================================================================= 3. decay and other limits
head("3. Limits on the majoron itself, each with its own error bar")
SUM_M2 = float(np.sum((MNU*1e9)**2))                                     # eV^2
def tau_nunu(mJ_eV, f): return HBAR/((mJ_eV*1e-9)*SUM_M2*1e-18/(16*math.pi*f*f))
def tau_gamgam(mJ_eV, f, trmD2=ME**2+MMU**2+MTAU**2):
    K = trmD2/(V174*f); x = (mJ_eV*1e-9)**2/(4*ME**2)
    amp = 0.5*x/3 if x < 0.5 else 0.5                                     # electron term of sum N_c T3 Q^2 g(x); the constant pieces cancel (no anomaly)
    return HBAR/(ALPHA**2*K*K*(mJ_eV*1e-9)**3/(4096*math.pi**7*V174**2)*amp*amp)
TAU_LIM = 250*GYR
pr(f"  J -> two neutrinos (Garcia-Cely & Heeck 2017): rate = m_J x sum(m_nu^2) / (16 pi f^2), sum(m_nu^2) = {SUM_M2:.3e} eV^2.")
pr(f"  Limit on dark matter decaying to invisible radiation: lifetime > about 250 Gyr = {TAU_LIM:.1e} s at 95%")
pr("  (Audren et al. 2014: 160 Gyr; Poulin, Serpico & Lesgourgues 2016: 170 Gyr; Nygaard, Tram & Hannestad 2021 and Simon et al. 2022: about 250 Gyr).")
pr(f"  J -> two photons (two loops; Garcia-Cely & Heeck 2017, Heeck & Patel 2019), with rule 5's links: tr(m_D m_D^+) = {ME**2+MMU**2+MTAU**2:.3f} GeV^2.")
pr("  X-ray line limits for 1 to 20 keV: lifetime to two photons above about 1e27 to 1e29 s (XMM-Newton, Chandra, NuSTAR).")
pr("   f (GeV)    angle   m_J (eV)     lifetime to neutrinos   over limit    lifetime to photons   over X-ray limit(1e28 s)")
for f, th in ((6.3e10, 1.0), (6.3e10, 1.81), (1.0e11, 1.0), (1.4e12, 1.0), (1.4e12, 0.1), (3e13, 1.0)):
    m = mJ_for_dm(f, th); t1, t2 = tau_nunu(m, f), tau_gamgam(m, f)
    pr(f"   {f:8.1e}   {th:4.2f}   {m:10.3e}   {t1:10.2e} s          {t1/TAU_LIM:8.1e}      {t2:10.2e} s       {t2/1e28:8.1e}")
pr("  Along the dark-matter line the neutrino lifetime grows as f^6: the shortest lifetime is at the smallest f.")
fl = 1e9
while tau_nunu(mJ_for_dm(fl, 2.5), fl) < TAU_LIM: fl *= 1.02
pr(f"  Lifetime limit alone (angle up to 2.5): f > {fl:.1e} GeV. The seesaw's own floor (link y_3 below sqrt(4 pi)) is higher, so the lifetime limit never bites.")
pr("  -> No X-ray line is expected: the photon lifetime is more than 13 powers of ten beyond today's reach. the earlier 'faint X-ray line' is withdrawn.")
pr("  -> The neutrino line sits at E = m_J/2 (keV or below), far under any neutrino detector's threshold (Borexino, KamLAND, Super-K reach MeV). Not testable.")
BISO, AS = 0.038, 2.1e-9
pr(f"  Isocurvature (Planck 2018: share below {BISO} at 95%): if the spot is settled before inflation, inflation's own jitter H_I/(2 pi) shakes the majoron wave.")
kiso = math.pi*math.sqrt(BISO*AS/(1-BISO))
for f, th in ((6.3e10, 1.0), (1.0e11, 1.0), (1.4e12, 0.1), (3e13, 1.0)):
    HI = kiso*f*th; r = (HI/(math.pi*MRED))**2*2/AS
    pr(f"   f = {f:.1e}, angle {th}: inflation's expansion rate H_I < {HI:.1e} GeV, i.e. primordial gravity waves r < {r:.0e} (LiteBIRD and CMB-S4 reach 1e-3)")
pr("  -> Prediction of the cold (pre-inflation) start: NO primordial gravity waves will be found. A detection of r above about 1e-3 would break it.")

# ================================================================= 4. seesaw, CP asymmetries, the non-thermal matter excess
head("4. The matter excess from a cold start: the hidden partners are handed over directly, then decay (Lazarides & Shafi 1991; Asaka et al. 1999)")
D = np.diag([ME, MMU, MTAU])
def pmns(s12, s13, s23, d, a21, a31):
    c12, c13, c23 = math.sqrt(1-s12), math.sqrt(1-s13), math.sqrt(1-s23); s12, s13, s23 = math.sqrt(s12), math.sqrt(s13), math.sqrt(s23)
    e = np.exp(1j*d)
    U = np.array([[c12*c13, s12*c13, s13/e], [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, s23*c13],
                  [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13]])
    return U @ np.diag([1, np.exp(1j*a21/2), np.exp(1j*a31/2)])
def takagi(A):
    U, S, Vh = np.linalg.svd(A); Z = U.conj().T @ Vh.T
    return S, U @ np.diag(np.sqrt(np.diag(Z)))
def herm(p):
    return np.array([[p[0], p[3]+1j*p[4], p[5]+1j*p[6]], [p[3]-1j*p[4], p[1], p[7]+1j*p[8]], [p[5]-1j*p[6], p[7]-1j*p[8], p[2]]])
def build(VL, d, a21, a31, s23, s12=0.307, s13=0.02195):          # same construction as matter_excess_rule5_exact_2026-11-14.py (rule 5 exact)
    U = pmns(s12, s13, s23, d, a21, a31); mnu = U @ np.diag(MNU) @ U.T
    MR = -D @ VL.conj() @ np.linalg.inv(mnu) @ VL.conj().T @ D
    S, W = takagi(MR); o = np.argsort(S); S = S[o]; W = W[:, o]
    return S, VL.conj().T @ D @ W.conj()/V174
def fvs(x): return math.sqrt(x)*(1/(1-x) + 1 - (1+x)*math.log((1+x)/x))
def cp_asym(S, h):
    """eps[i, alpha]: CP asymmetry of partner i into flavour alpha (Covi, Roulet & Vissani 1996; flavoured form Abada et al. 2006)."""
    hh = h.conj().T @ h; eps = np.zeros((3, 3))
    for i in range(3):
        for a in range(3):
            t = 0.0
            for j in range(3):
                if j == i: continue
                x = (S[j]/S[i])**2
                t += (np.conj(h[a, i])*h[a, j]*hh[i, j]).imag*fvs(x) + (np.conj(h[a, i])*h[a, j]*hh[j, i]).imag/(1-x)
            eps[i, a] = t/(8*math.pi*hh[i, i].real)
    return eps, hh
MSTAR = 1.07e-12                                                     # GeV, equilibrium neutrino mass
def nonthermal(S, h, TR_cap_frac=0.1):
    """For each source partner i: effective asymmetry after the lightest partner's later washout, and the handover number X = Br x share x T_R/m needed."""
    eps, hh = cp_asym(S, h); out = []
    K1a = np.abs(h[:, 0])**2*V174**2/S[0]/MSTAR                      # the lightest partner's washout per flavour (acts at T ~ M_1 if T_R > M_1)
    for i in (1, 2):
        e_raw = eps[i].sum(); e_eff = float(np.sum(eps[i]*np.exp(-3*math.pi/8*K1a)))
        X = YB/((28/79)*1.5*abs(e_eff)) if e_eff != 0 else float("inf")
        Xmax = TR_cap_frac*S[1]/(2*S[i])                              # T_R < M_2/10 (no heat-made middle partner), m_phi > 2 M_i
        out.append((i+1, e_raw, e_eff, X, Xmax))
    return out, K1a, eps
pr(f"  Needed: baryons per entropy {YB:.3e} (eta_B = {ETA_B:.2e} +- {ETA_B_ERR:.0e}). Handover formula: n_B/s = (28/79) x (3/2) x eps x X,")
pr("  where X = (branching to that partner) x (the source's share of the energy) x T_reheat / m_source.  X can never exceed 1/20 x M_2/M_i")
pr("  (the source must outweigh two partners, and the heat must stay 10x below the middle partner so nothing is made or erased by heat).")
I3 = np.eye(3, dtype=complex)
def tm1_delta(s23, s13=0.02195):
    c2 = 1-2*s23; cot = c2/(2*math.sqrt(s23*(1-s23)))
    cd = -cot*(1-5*s13)/(2*math.sqrt(2)*math.sqrt(s13)*math.sqrt(1-3*s13))
    return 360 - math.degrees(math.acos(cd))
rng = np.random.default_rng(82)
def scan(VL, n, tm1=False):
    res = []
    for _ in range(n):
        a, b = rng.uniform(0, 2*math.pi, 2); s23 = float(rng.uniform(0.435, 0.585))
        if tm1: d = math.radians(tm1_delta(s23)); s12 = 0.318
        else: d = float(rng.uniform(0, 2*math.pi)); s12 = 0.307
        S, h = build(VL, d, a, b, s23, s12=s12); o, K1a, eps = nonthermal(S, h)
        res.append((S, o, K1a, (d, a, b, s23)))
    return res
def report(name, res):
    for k, lab in ((0, "middle partner"), (1, "heaviest partner")):
        e = np.array([abs(r[1][k][2]) for r in res]); X = np.array([r[1][k][3] for r in res]); Xm = np.array([r[1][k][4] for r in res])
        ok = (X <= Xm).mean()*100; ib = int(np.argmin(X/Xm)); S = res[ib][0]
        pr(f"   {name:34s} source = {lab:16s}: |eps| median {np.median(e):.1e}, largest {e.max():.1e}; best X needed {X[ib]:.1e} against ceiling {Xm[ib]:.1e}"
           f" ({'REACHED' if X[ib] <= Xm[ib] else 'short by %.0fx' % (X[ib]/Xm[ib])}); settings that reach it: {ok:.0f}%")
    S = res[0][0]; pr(f"     partners (first setting): {S[0]:.2e}, {S[1]:.2e}, {S[2]:.2e} GeV; lightest partner's washout K(e, mu, tau) = {np.round(res[0][2], 2).tolist()}")
pr("\n  (a) Rule 5 read literally (same perspective), 400 random settings of the unmeasured phases, theta23 over its 3-sigma range:")
resA = scan(I3, 400); report("literal rule 5, phases free", resA)
resB = scan(I3, 400, tm1=True); report("literal rule 5, TM1 delta (~270)", resB)
pr("\n  (b) The six ways of pairing families (parameter-free), phases free:")
import itertools
names = ["e", "mu", "tau"]
for perm in itertools.permutations(range(3)):
    P = np.zeros((3, 3), dtype=complex)
    for i, j in enumerate(perm): P[i, j] = 1
    r = scan(P.T, 150); X = min(min(q[1][0][3]/q[1][0][4], q[1][1][3]/q[1][1][4]) for q in r)
    frac = np.mean([min(q[1][0][3]/q[1][0][4], q[1][1][3]/q[1][1][4]) <= 1 for q in r])*100
    pr(f"   nu_e~m_{names[perm[0]]}, nu_mu~m_{names[perm[1]]}, nu_tau~m_{names[perm[2]]}: partners {r[0][0][0]:.0e}, {r[0][0][1]:.0e}, {r[0][0][2]:.0e} GeV; best X/ceiling = {X:.2e}; reach it: {frac:.0f}%")
pr("\n  (c) Random perspectives (uniform over all rotations), random phases, 1500 tries: how common is a cold-start success?")
def haar(rng):
    Z = (rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))/math.sqrt(2); Q, R = np.linalg.qr(Z); return Q @ np.diag(np.diag(R)/np.abs(np.diag(R)))
okc = 0; okc_free = 0; best = None; Ms3 = []
for _ in range(1500):
    VL = haar(rng); d, a, b = rng.uniform(0, 2*math.pi, 3); s23 = float(rng.uniform(0.435, 0.585))
    S, h = build(VL, d, a, b, s23); o, K1a, eps = nonthermal(S, h); q = min(o[0][3]/o[0][4], o[1][3]/o[1][4]); Ms3.append(S[2])
    okc += q <= 1
pr(f"   reach the observed excess with an allowed handover: {okc/15:.1f}%  (thermal hot start, an earlier test: 0.0% reach it, 0.6% within 10x). Heaviest partner: median {np.median(Ms3):.1e} GeV.")
pr("\n  (d) the earlier small-tilt corner (the one thermal survivor), now with a cold start:")
pts = json.load(open(os.path.join(HERE, "matter_excess_best_points_2026-11-14.json")))
x = np.array(pts[1]["x"]); VLc = expm(1j*herm(x[:9])); Sc, hc = build(VLc, x[9], x[10], x[11], x[12]); oc, K1c, epsc = nonthermal(Sc, hc)
pr(f"   partners {Sc[0]:.2e}, {Sc[1]:.2e}, {Sc[2]:.2e} GeV; lightest partner's washout K(e, mu, tau) = {np.round(K1c, 2).tolist()}")
for (i, er, ee, X, Xm) in oc:
    pr(f"   source = partner {i}: eps = {er:+.2e} (after the lightest partner's washout {ee:+.2e}); X needed = {X:.2e}; ceiling {Xm:.2e}: {'REACHED with room %.0fx' % (Xm/X) if X <= Xm else 'short by %.1fx' % (X/Xm)}")
EPS_C = abs(oc[0][2]); X_C = oc[0][3]
from scipy.optimize import differential_evolution
def lit_obj(p):
    S, h = build(I3, p[0], p[1], p[2], p[3]); o, _, _ = nonthermal(S, h); return min(o[0][3]/o[0][4], o[1][3]/o[1][4])
de = differential_evolution(lit_obj, [(0, 2*math.pi)]*3 + [(0.435, 0.585)], seed=7, popsize=20, maxiter=80, tol=1e-8, polish=True)
S_o, h_o = build(I3, *de.x); o_o, K_o, _ = nonthermal(S_o, h_o)
pr(f"\n  Literal rule 5, directed search over the three phases and theta23 (differential evolution): best X needed / ceiling = {de.fun:.2f}"
   f" ({'REACHED' if de.fun <= 1 else 'still short by %.1fx: FAIL' % de.fun}) at delta, a21, a31 = {[round(math.degrees(v) % 360) for v in de.x[:3]]} deg, sin^2(th23) = {de.x[3]:.3f};"
   f" partners {S_o[0]:.1e}, {S_o[1]:.1e}, {S_o[2]:.1e} GeV")
LIT_BEST = de.fun
# literal best
ibA = int(np.argmin([min(r[1][0][3]/r[1][0][4], r[1][1][3]/r[1][1][4]) for r in resA])); bestA = resA[ibA]
kA = 0 if bestA[1][0][3]/bestA[1][0][4] < bestA[1][1][3]/bestA[1][1][4] else 1
EPS_L, X_L, XM_L, S_L = abs(bestA[1][kA][2]), bestA[1][kA][3], bestA[1][kA][4], bestA[0]
pr(f"\n  Literal rule 5, best setting: source partner {kA+2}, eps = {EPS_L:.2e}, X needed {X_L:.2e}, ceiling {XM_L:.2e}; phases delta, a21, a31 = "
   f"{[round(math.degrees(v) % 360) for v in bestA[3][:3]]} deg, sin^2(th23) = {bestA[3][3]:.3f}; partners {S_L[0]:.1e}, {S_L[1]:.1e}, {S_L[2]:.1e} GeV")
pr(f"  Davidson-Ibarra ceiling on eps for a partner of mass M: 3 M m_3/(16 pi v^2) = {3*MNU[2]/(16*math.pi*V174**2):.2e} x M/GeV.")
pr(f"  Share of partner decays that go to a majoron instead (N -> nu J): about v^2/(2 f^2) = {V174**2/(2*F0**2):.0e} at f = 6.3e10. Negligible.")

# ================================================================= 5. one event for both
head("5. ONE settling event for both numbers: three versions, each pushed until it works or fails")
pr("  Version A. Hot start + majoron wave.  Matter excess by heat, dark matter by the wave.")
pr(f"    f near the heaviest partner: hot majorons (section 2). FAIL.   f above about {F_HOT:.0e} GeV: no hot majorons, but the matter excess is then")
pr("    only 0.32 to 0.39 of observed (cold thermal start, since reading (i) has no right-side force to pre-fill the partners; earlier tests): FAIL by 2.6 to 3x.")
pr(f"    (eta_B is measured to {ETA_B_ERR/ETA_B*100:.1f}%; the rate equations are good to tens of per cent; a 2.6x gap is outside both.)")
pr("\n  Version B. The settling piece itself (rho) decays and makes BOTH: partners (-> matter excess) and majorons (-> dark matter).")
pr("    Rates: rho -> J J = m_rho^3/(32 pi f^2);  rho -> N_i N_i = M_i^2 m_rho (1 - 4M_i^2/m_rho^2)^(3/2)/(16 pi f^2).")
xs = np.linspace(0.01, 0.99, 981); R = xs/2*(1-xs)**1.5; RMAX = float(R.max()); XBEST = float(xs[np.argmax(R)])
pr(f"    The ratio (partners)/(majorons) is at most {RMAX:.4f} (at m_rho = {2/math.sqrt(XBEST):.2f} M_i): the settling piece makes at least {100/(1+RMAX):.0f}% majorons.")
pr(f"    Then the ratio of the two cosmic amounts is fixed, with the handover number cancelling out:")
kB = (OM_DM/OM_B)*MP*(28/79)
pr(f"      m_J = (Omega_DM/Omega_B) x m_proton x (28/79) x eps x (partners/majorons) = {kB:.3f} GeV x eps x ratio  <=  {kB*RMAX*1e3:.0f} eV x (eps/1e-6)")
def warm_floor(p_over_T, mwdm=5.3, g=106.75): return 7.2*(mwdm/3.0)**(4/3)*(p_over_T/3.15)*(106.75/g)**(1/3)      # keV (Kamada & Yanagi 2019 mapping)
pr(f"    Those majorons are born fast (momentum m_rho/2). Lyman-alpha (thermal relic > 5.3 keV) maps to m_J > {warm_floor(1.0):.1f} keV x (momentum/temperature at reheating).")
pr("     case                         eps       m_J (all DM)   momentum/T allowed   needs expansion between decay and reheating   wave must not overfill: angle <")
for lab, eps, M2, M3 in (("small-tilt corner, partner 2", EPS_C, Sc[1], Sc[2]), ("literal rule 5, best", EPS_L, S_L[kA+1], S_L[2]),
                         ("Davidson-Ibarra ceiling at 1.6e10", 3*MNU[2]/(16*math.pi*V174**2)*1.63e10, 1.63e10, 1.0e11)):
    mJ = kB*eps*RMAX*1e9; pT = mJ/1e3/warm_floor(1.0); Xn = YB/((28/79)*1.5*eps); BrN = RMAX/(1+RMAX)
    TR_over_m = Xn/(BrN*0.5)                                          # share r = 0.5 at most (the start field must dominate)
    stretch = (1/(2*TR_over_m))/pT                                    # a_R/a_dec needed
    fmin = math.sqrt(2)*M3
    thmax = (mJ_for_dm(fmin, 1.0)/mJ)**0.25
    pr(f"     {lab:28s} {eps:8.1e}   {mJ:9.1f} eV    {pT:10.2e}        {'IMPOSSIBLE (T_R/m_rho = %.1f > 1/20)' % TR_over_m if TR_over_m > 0.05 else 'stretch %.0e (%.0f e-folds)' % (stretch, math.log(stretch)):38s}     {('%.2f' % thmax) if thmax < 3.1 else 'any'} (at f = sqrt2 M_3 = {fmin:.1e})")
pr("    -> It can be made to work only with a tilted perspective (eps near 1e-7 or more), a long cold stretch (about 8 e-folds) between the settling piece's decay")
pr(f"       and reheating, and a settling piece holding about half of the energy. Then m_J is PREDICTED: {kB*EPS_C*RMAX*1e9:.0f} eV in the small-tilt corner, at most {kB*3*MNU[2]/(16*math.pi*V174**2)*1.63e10*RMAX*1e9:.0f} eV.")
pr("       Numbers used: eps (from the phases and perspective), m_rho/M_i, the share, the stretch, the angle: 5 for 2 data. LEAD at best, not a derivation.")
pr("       Without the long stretch the same majorons are hot: FAIL (this is the earlier failure, now with its exact reason).")
pr("\n  Version C. Cold start, spot settled before inflation (reading (i): 'the spot is exact always'). Partners handed over by the start; majoron wave = dark matter.")
pr("   f (GeV)    y_3 (M_3=1.0e11)   m_J (angle 1)   lifetime/limit   H_I below     T_reheat window (GeV)            hot-majoron share at the top of the window")
for f in (4.0e10, 6.3e10, 1.0e11, 2.0e11, 5e11, 1.4e12):
    m = mJ_for_dm(f); Tos = omega_mis(m, f)[1]; y3 = math.sqrt(2)*Sc[2]/f
    pr(f"   {f:8.1e}   {y3:6.2f}{' (too big)' if y3 > math.sqrt(4*math.pi) else '          '}   {m:10.3e} eV   {tau_nunu(m, f)/TAU_LIM:8.1e}      {kiso*f:8.1e}     {Tos:.1e} (wave starts) to {Sc[1]/10:.1e}   {Y_hot(f, Sc[1]/10, list(Sc))*m*1e-9/RHO_OVER_S_DM:.1e}")
pr(f"   Matter excess in every row: handover X = {X_C:.2e} (small-tilt corner), e.g. T_reheat/m_source = {X_C/0.1:.1e} if (branching x share) = 0.1. Ceiling 0.05.")
pr("   Every row passes every limit. The region is a 2-surface: one relation fixes m_J x angle^4 x f^4, one fixes X x eps. Nothing else is fixed.")

# ================================================================= 6. does TRPG fix any number?
head("6. Does TRPG fix any of these numbers by its own rules? With luck checks")
pr("  6a. f tied to the heaviest partner (link y_3 = 1 in M = y f, or sqrt2 in M = y f/sqrt2):")
for lab, M3 in (("an earlier test, generation by generation", 6.3e10), ("small-tilt corner", float(Sc[2])), ("literal rule 5 with the measured mixing", float(np.median([r[0][2] for r in resA])))):
    pr(f"    {lab:42s} M_3 = {M3:.2e}: f = M_3 gives m_J = {mJ_for_dm(M3)/1e3:9.3g} keV / angle^4;  smallest allowed f = {math.sqrt(2)*M3/math.sqrt(4*math.pi):.1e} gives {mJ_for_dm(math.sqrt(2)*M3/math.sqrt(4*math.pi))/1e3:9.3g} keV")
pr("    The three readings of 'the heaviest partner' differ by 16x in f, i.e. 65,000x in m_J. And nothing in the rules sets y_3 = 1. NOT FIXED.")
pr(f"    Robust statement that does follow: f cannot be below M_3/2.5, so m_J x angle^4 has a CEILING: {mJ_for_dm(math.sqrt(2)*6.3e10/math.sqrt(4*math.pi))/1e3:.0f} keV for the earlier M_3, "
   f"{mJ_for_dm(math.sqrt(2)*float(Sc[2])/math.sqrt(4*math.pi))/1e3:.0f} keV for the corner, {mJ_for_dm(math.sqrt(2)*1.44e12/math.sqrt(4*math.pi)):.1f} eV for the literal 1.4e12.")
pr("\n  6b. Gravity's leak redone with two corrections: the exact oscillator, and the fact that a power-n leak has n valleys,")
pr("      so the wave's swing is n times smaller (period 2 pi f/n). V = c (f/sqrt2)^n/M^(n-4) cos(n J/f), c = 1.")
def m_leak(n, f, M): return math.sqrt(n*n*(f/math.sqrt(2))**n/(M**(n-4)*f*f))*1e9     # eV
pr("      scale         n    m_J at f=6.3e10    wave's share of dark matter at f=6.3e10 (angle 1; 2.5)     f where it is ALL the dark matter (angle 1; 2.5)   m_J there      lifetime there")
sols = []
for Mn, M in (("Planck", MPL), ("reduced Planck", MRED)):
    for n in range(5, 13):
        m = m_leak(n, F0, M); sh1, sh2 = omega_mis(m, F0, 1.0, n)[0]/OM_DM, omega_mis(m, F0, 2.5, n)[0]/OM_DM
        fs = []
        for th in (1.0, 2.5):
            lo, hi = 6.0, 17.0
            for _ in range(60):
                mid = 0.5*(lo+hi)
                if omega_mis(m_leak(n, 10**mid, M), 10**mid, th, n)[0] < OM_DM: lo = mid
                else: hi = mid
            fs.append(10**hi)
        sols.append((Mn, n, fs))
        pr(f"      {Mn:14s} {n:2d}   {m:10.2e} eV      {sh1:9.2e}; {sh2:9.2e}                                {fs[0]:.2e}; {fs[1]:.2e}                         {m_leak(n, fs[0], M):.2e} eV    {tau_nunu(m_leak(n, fs[0], M), fs[0]):.1e} s{' (too short)' if tau_nunu(m_leak(n, fs[0], M), fs[0]) < TAU_LIM else ''}")
pr(f"      -> At f = 6.3e10 GeV the n = 8 leak (3.4 keV) gives a wave that is only {omega_mis(m_leak(8, F0, MPL), F0, 1.0, 8)[0]/OM_DM*100:.1f}% of the dark matter (angle 1): the earlier 'lands' does NOT survive.")
pr("      Each power picks its own f. Does any land on a partner scale the theory has? Luck check: draw f at random (log-uniform, 1e9 to 1e13 GeV);")
rngL = np.random.default_rng(5); hits = 0; NL = 4000
for _ in range(NL):
    lf = rngL.uniform(9, 13)
    hits += any(min(math.log10(fs[0]), math.log10(fs[1])) - 0.15 <= lf <= max(math.log10(fs[0]), math.log10(fs[1])) + 0.15 for (_, n, fs) in sols)
pr(f"      some power and some gravity scale reproduces all of the dark matter (angle 1 to 2.5, f within 1.4x) for {hits/NL*100:.0f}% of random f. A landing is not evidence.")
for lab, M3 in (("6.3e10", 6.3e10), ("1.0e11", 1.0e11), ("1.0e12", 1.0e12)):
    near = [(Mn, n) for (Mn, n, fs) in sols if min(fs)/1.4 <= M3 <= max(fs)*1.4 or min(fs)/1.4 <= math.sqrt(2)*M3 <= max(fs)*1.4]
    pr(f"      f = M_3 or sqrt2 M_3 with M_3 = {lab}: powers that land: {near if near else 'none'}")
pr("      (Planck n = 7: a 35 MeV majoron living 1e21 s; neutrino-line searches at Super-K and KamLAND reach lifetimes of that order there (Garcia-Cely & Heeck 2017): borderline.)")
pr("      And the lower powers must be switched off by hand, and post-inflation settling with n > 1 leaves domain walls (Sikivie 1982). NOT A RULE. FAIL as a derivation.")
pr("\n  6c. Mass-as-balance candidates for m_J (an earlier test style), with nothing fitted:")
cands = {"f^2/M_Planck": F0**2/MPL, "M_3 x m_nu3 / f (seesaw mirror)": 6.3e10*MNU[2]/F0, "sqrt(m_nu3 x m_e)": math.sqrt(MNU[2]*ME), "m_e^2/v": ME**2/V246,
         "v^2/f": V246**2/F0, "m_nu3 x f/M_Planck": MNU[2]*F0/MPL, "(v/f) x m_e": V246/F0*ME, "sqrt(m_nu3 x v)": math.sqrt(MNU[2]*V246)}
tgt = mJ_for_dm(F0, 1.0); lo_t, hi_t = mJ_for_dm(F0, 2.5), mJ_for_dm(F0, 0.5)
nl = 0
for k, v in cands.items():
    ok = lo_t <= v*1e9 <= hi_t; nl += ok
    pr(f"      {k:34s} = {v*1e9:10.3e} eV   {'inside' if ok else 'outside'} the window {lo_t:.1e} to {hi_t:.1e} eV (angle 0.5 to 2.5 at f = 6.3e10)")
wid = math.log10(hi_t/lo_t); span = math.log10(max(cands.values())/min(cands.values()))
pr(f"      {nl} of {len(cands)} land; the window is {wid:.1f} decades wide and the candidates span {span:.0f} decades, so about {len(cands)*wid/span:.1f} are expected by luck. NOT EVIDENCE.")
pr("\n  6d. The matter excess: does a rule fix eps x X?  eps follows from the phases and the perspective (not fixed: an earlier test). X is a start condition")
pr("      (how the universe's start handed over its energy). Scope: the theory does not have to solve for the bang. X is therefore outside it.")

# ================================================================= 7. neutrino-side consistency
head("7. Consistency of the neutrino numbers the picture leans on")
sm = MNU.sum()*1e9
pr(f"  Sum of masses {sm:.4f} eV. Smallest possible for normal ordering (NuFIT 6.0 splittings): {math.sqrt(7.49e-5)+math.sqrt(2.513e-3):.4f} eV.")
pr(f"  DESI DR1 + CMB (2024): below 0.072 eV (95%). DESI DR2 + CMB (2025): below 0.064 eV (95%). TRPG sits at {sm/0.064*100:.0f}% of the tighter limit: inside, close.")
pr("  (Those fits peak at or below zero, so ANY normal-ordering theory is in about 2 to 3 sigma tension with them; TRPG is as low as normal ordering allows. Inverted ordering, 0.10 eV, is out.)")
pr("  delta_CP from TM1 depends on theta23 (Albright & Rodejohann 2009): cos(delta) = -cot(2 theta23)(1 - 5 s13^2)/(2 sqrt2 s13 sqrt(1 - 3 s13^2))")
for s23 in (0.435, 0.470, 0.500, 0.561, 0.585): pr(f"    sin^2(theta23) = {s23:.3f}: delta = {tm1_delta(s23):.1f} deg")
pr("  NuFIT 6.0 (normal ordering): best fit 177 (+19, -20) deg without Super-K atmospheric data, 212 (+26, -41) deg with it; 3-sigma ranges 96 to 422 and 124 to 364 deg.")
pr(f"  TM1's band {tm1_delta(0.435):.0f} to {tm1_delta(0.585):.0f} deg is inside both 3-sigma ranges; about {(tm1_delta(0.470)-212)/26:.1f} error bars above the with-Super-K best fit,")
pr("  and more than 3 naive error bars above the without-Super-K best fit (the likelihood is far from Gaussian; the 3-sigma range still contains it). Allowed, not favoured.")
pr("  TM1 sin^2(theta12) = 0.318 against NuFIT 6.0 0.308 +- 0.012: +0.8 sigma. JUNO's first 59 days (2025) gave 0.3092 +- 0.0087: TM1 at +1.0 sigma.")
pr("  JUNO's full data set reaches about 0.5%: 0.318 +- 0.002 would confirm TM1, 0.308 +- 0.002 would break it at 5 sigma.")
def mee(d, a21, a31, s12=0.307, s13=0.02195):
    return abs((1-s12)*(1-s13)*MNU[0] + s12*(1-s13)*MNU[1]*np.exp(1j*a21) + s13*MNU[2]*np.exp(1j*(a31-2*d)))*1e12
g = np.linspace(0, 2*math.pi, 181); vals = np.array([[mee(0, a, b) for b in g] for a in g])
pr(f"  Neutrinoless double-beta decay: effective mass {vals.min():.1f} to {vals.max():.1f} meV (limit 28 to 122 meV, KamLAND-Zen 2024; nEXO and LEGEND-1000 reach about 5 to 20 meV).")
pr("  A signal above about 5 meV would break the picture; a null result at 10 meV leaves it standing. In a majoron picture a second mode exists (decay WITH a majoron),")
pr(f"  with coupling m_ee/f = {vals.max()*1e-12/F0:.0e}: 20 powers of ten under today's limit of about 1e-5. Not testable.")

# ================================================================= 8. count and verdict
head("8. Honest count and verdict")
pr("  Data: 2 numbers (eta_B, Omega_DM h^2).")
pr("  Free numbers in the working version (C): f, m_J, the starting angle (dark matter: 3 numbers, 1 relation);")
pr("     the handover X and the CP asymmetry eps (matter excess: eps needs the perspective and the 3 unmeasured phases; 1 relation). At least 5 for 2 data.")
pr("  So neither number is derived. Both are ALLOWED, in one consistent cold-start picture that passes every limit checked.")
pr("  What changed against earlier tests:  (1) the 1.4 keV is 8.7 keV with the exact oscillator (and it is a relation, not a mass);")
pr("     (2) gravity's leak at n = 8 does not give the dark matter at f = 6.3e10 once its 8 valleys are counted;")
pr(f"     (3) no X-ray line; (4) with a cold start {okc/15:.0f}% of random perspectives reach the matter excess (0% with a hot start); rule 5 read literally still fails, by {LIT_BEST:.1f}x;")
pr(f"     (5) the hot start is not excluded at f > {F_HOT:.0e} GeV, but then the matter excess is 2.6 to 3x short;")
pr(f"     (6) the settling piece decays at least 91% to majorons: as the single source it works only with a long cold stretch, and then m_J = {kB*EPS_C*RMAX*1e9:.0f} eV (corner).")
