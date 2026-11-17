"""Falsification pass, an earlier test: the matter excess (14 Nov 2026).
Attacks on the 'working corner' of an earlier test (rule 5 exact, perspective V_L free):
  A. the literal rule (no free perspective);  B. how rare the corner is;  C. neutrino data at the corner;
  D. start: partners present from the start (thermal) vs made only by their own links (cold start);
  E. the right side's own force carriers (W_R), which the standard code leaves out, added at rough strength over a wide band;
  F. which partner makes the excess;  G. what the corner demands of the partners' couplings.
Solver: ULYSSES 3.0.2 (Granelli, Moffat, Perez-Gonzalez, Schulz & Turner 2021/2023), density-matrix equations, with one added term (E).
Credit: Fukugita & Yanagida 1986; Sakharov 1967; Ma, Sarkar & Sarkar 1999 and Frere, Hambye & Vertongen 2009 (right-side carriers dilute leptogenesis);
Blanchet, Di Bari, Jones & Marzola 2013; NuFIT 6.0 (Esteban et al. 2024); DESI 2025; KamLAND-Zen 2024; Planck 2018."""
import math, json, os, sys, warnings, importlib.util
warnings.filterwarnings("ignore")
import numpy as np, ulysses
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("r5", os.path.join(HERE, "matter_excess_rule5_exact_2026-11-14.py")); r5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)
from ulysses.etab3DME import fast_RHS
from odeintw import odeintw
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
OBS = 6.1e-10
MI, A2R = 10**11.74, 1/51.5            # right side's scale and coupling there (falsification_pass, minimal running)
GR2 = (math.sqrt(2)*4*math.pi*A2R/(8*MI**2))**2      # (Fermi-type constant of the right-handed force)^2, GeV^-4

class WR(ulysses.EtaB_3DME):
    """Same equations, plus (i) an extra rate that pulls each heavy partner toward equilibrium through the right side's carriers,
    Gamma ~ c G_R^2 T^5 (T > M), c G_R^2 T^3 M^2 (T < M), capped above the carriers' own mass; and (ii) switches for each partner's source."""
    cmul = 0.0; src = (1.0, 1.0, 1.0)
    def RHS(self, y0, zzz, ETA, _C, K, _W):
        ETA = list(ETA)
        for i in range(3):
            for j in range(6): ETA[6*i+j] = ETA[6*i+j]*self.src[i]
        out = list(super().RHS(y0, zzz, ETA, _C, K, _W))
        if self.cmul > 0:
            Ms = [self.M1, self.M2, self.M3]; neq = [self._n1eq, self._n2eq, self._n3eq]; T = self.M1/zzz
            for i in range(3):
                M = float(np.real(Ms[i])); H = 1.66*math.sqrt(self.gstar)*T*T/self.MP
                if T > MI: G = GR2*MI**4*T
                elif T > M: G = GR2*T**5
                else: G = GR2*T**3*M*M
                out[i] = out[i] - self.cmul*G/(H*zzz)*(y0[i] - neq[i])
        return out

def eta(x, init=1, cmul=0.0, src=(1, 1, 1), zsteps=1500, model=None, full=False, Tmax=None):
    VL = expm(1j*r5.herm(np.array(x[:9]))); S, h = r5.build(VL, x[9], x[10], x[11], x[12])
    cls = WR if model is None else getattr(ulysses, model)
    L = cls(zmin=(S[0]/Tmax if Tmax else min(1e-3, S[0]/(30*S[2]))), zmax=1000, zsteps=zsteps, ordering=0, initial_abundance=init, which_param="manual")
    L.cmul = cmul; L.src = src
    pd = dict(m=math.log10(r5.MNU[0]*1e9), M1=math.log10(S[0]), M2=math.log10(S[1]), M3=math.log10(S[2]))
    for i in range(3):
        for j in range(3): pd[f"Y{i+1}{j+1}_mag"] = abs(h[i, j]); pd[f"Y{i+1}{j+1}_phs"] = float(np.angle(h[i, j]))
    try: e = float(L(pd))
    except Exception as ex: e = float("nan")
    return (e, S, h, L) if full else e

pts = json.load(open(os.path.join(HERE, "matter_excess_best_points_2026-11-14.json")))
names = ["corner 1 (largest found; perspective tilted 64 deg)", "corner 2 (small tilt, 25 deg)"]

head("A. The literal rule 5 (the neutrino's link is the electron family's, same perspective, nothing free but the unmeasured phases)")
rng = np.random.default_rng(34); vals = []
for _ in range(200):
    d, a, b = rng.uniform(0, 2*math.pi, 3); s23 = float(rng.uniform(0.435, 0.585)); vals.append(abs(r5.eta(np.eye(3, dtype=complex), d, a, b, s23)))
pr(f"  200 random settings of the unmeasured phases: largest matter excess {max(vals):.1e}; observed {OBS:.1e}: short by {OBS/max(vals):,.0f}x")
pr("  -> As a PREDICTION of the rules as written, the matter excess FAILS. Nothing below changes that; it only asks whether the excess is still ALLOWED.")

head("B. How rare is a perspective that works? (saved scan of an earlier test: 3,000 random perspectives, each with random phases)")
sc = json.load(open(os.path.join(HERE, "matter_excess_random_perspectives_2026-11-14.json")))
arr = np.abs(np.array(sc["eta"], dtype=float))
pr(f"  {len(arr)} random perspectives: median {np.median(arr):.1e}; largest {arr.max():.1e}; reach the observed value: {100*(arr >= OBS).mean():.2f}%; within 10x: {100*(arr >= OBS/10).mean():.2f}%")

head("C. The two working corners against neutrino data")
nufit = dict(d=(124, 364), s23=(0.435, 0.585))     # NuFIT 6.0, normal ordering, 3 sigma (with SK atmospheric)
for nm, p in zip(names, pts):
    x = p["x"]; e, S, h, L = eta(x, full=True); d = math.degrees(x[9]) % 360; s23 = x[12]
    U = r5.pmns(0.307, 0.02195, s23, x[9], x[10], x[11]); mee = abs((U[0]**2*r5.MNU).sum())*1e12
    dd = d if d >= 124 else d + 360
    pr(f"  {nm}: excess {e:.2e}; partners {S[0]:.2e}, {S[1]:.2e}, {S[2]:.2e} GeV; washout strengths K = {L.k1.real:.1f}, {L.k2.real:.1f}, {L.k3.real:.1f}")
    pr(f"     delta_CP = {d:.0f} deg (allowed 124-364 at 3 sigma: {'inside' if 124 <= dd <= 364 else 'OUTSIDE'}{', at the very edge' if min(abs(dd-124), abs(dd-364)) < 8 else ''}); "
       f"sin^2(theta23) = {s23:.3f} (allowed 0.435-0.585: {'inside' if 0.435 <= s23 <= 0.585 else 'OUTSIDE'}{', at the very edge' if min(abs(s23-0.435), abs(s23-0.585)) < 0.003 else ''})")
    pr(f"     double-beta mass m_ee = {mee:.1f} meV (limit 28-122 meV: fine)")
pr(f"  sum of neutrino masses = {r5.MNU.sum()*1e9:.4f} eV (cosmology: under about 0.064-0.072 eV: fine, but close)")
pr("  -> Corner 2 passes every neutrino test. Corner 1 sits on two 3-sigma edges at once: counted as failing. Only corner 2 is carried forward.")

head("D. Start: were the heavy partners already there?")
for nm, p in zip(names, pts):
    x = p["x"]; a = eta(x, init=1); b = eta(x, init=0); c = eta(x, init=1, model="EtaB_3DMEsct")
    pr(f"  {nm}: thermal start {a:.2e} ({a/OBS:.2f} of observed) | cold start {b:.2e} ({b/OBS:.2f}) | thermal, with Higgs/top scatterings (3DMEsct) {c:.2e} ({c/OBS:.2f})")
pr("  The thermal start is what the 4x4 picture gives, IF the universe was hot enough for the right side's carriers to make the partners.")

head("E. The right side's own carriers (not in the standard code). They keep heavy partners closer to equilibrium, which lowers the excess.")
x = pts[1]["x"]; e0, S, h, L = eta(x, full=True)
H = lambda T: 1.66*math.sqrt(106.75)*T*T/1.22e19
pr(f"  Rough rate over expansion rate at T = M:  middle partner ({S[1]:.1e} GeV): {GR2*S[1]**5/H(S[1]):.2f} x c;   heaviest ({S[2]:.1e} GeV): {GR2*S[2]**5/H(S[2]):.1f} x c")
pr("  c is an unknown number of order 1-10 (channel count); the carriers' mass is uncertain by about 3x, which moves the rate by about 80x. So a wide band is scanned.")
pr(f"  {'strength c':>12s} {'corner 2 excess':>16s} {'of observed':>12s}")
band = {}
for c in (0, 0.1, 1, 3, 10, 100):
    e = eta(x, cmul=c); band[c] = e; pr(f"  {c:12g} {e:16.2e} {e/OBS:12.2f}")

head("E2. How hot must the universe have been? Cold start (no partners at first), right-side carriers switched on, corner 2")
pr("  The partners are then made only by the right side's carriers and their own links. Monopoles (falsification_pass 1g) cap the temperature near 5e11 GeV.")
pr(f"  {'hottest temperature':>22s} {'c = 1':>12s} {'c = 10':>12s}   (fraction of the observed excess)")
for Tm in (3e12, 5e11, 2e11, 1e11, 5e10, 2e10):
    r = [eta(x, init=0, cmul=c, Tmax=Tm)/OBS for c in (1, 10)]
    pr(f"  {Tm:22.0e} {r[0]:12.2f} {r[1]:12.2f}")

head("F. Which partner makes the excess at corner 2?")
for lab, src in (("lightest only", (1, 0, 0)), ("middle only", (0, 1, 0)), ("heaviest only", (0, 0, 1))):
    e = eta(x, src=src); pr(f"  source from the {lab:14s}: {e:+.2e} ({e/OBS:+.2f} of observed)")

head("G. What corner 2 demands of the partners")
pr(f"  partner masses / right side's scale ({MI:.1e} GeV): {S[0]/MI:.1e}, {S[1]/MI:.2f}, {S[2]/MI:.2f}  (the lightest needs a coupling {S[2]/S[0]:.0e} times smaller than the heaviest)")
hh = np.abs(h)**2; share = hh[:, 0]/hh[:, 0].sum()
pr(f"  lightest partner's link shares (e, mu, tau): {share[0]:.1e}, {share[1]:.3f}, {share[2]:.3f}  (the corner works only because the electron share is almost zero)")
json.dump(dict(band={str(k): v for k, v in band.items()}), open(os.path.join(HERE, "falsification_matter_excess_band_2026-11-14.json"), "w"))
