"""Matter excess with rule 5 held EXACTLY (14 Nov 2026).
Rule 5: the neutrino's Dirac link has the same sizes as the electron family's (m_e, m_mu, m_tau).
Then the hidden partners' mass matrix follows from the seesaw with nothing to choose:
    M_R = - D V_L* m_nu^-1 V_L^dagger D ,   D = diag(m_e, m_mu, m_tau),   m_nu = U diag(m1, m2, m3) U^T
What is free: V_L, the neutrino link's 'perspective' relative to the electron family's (a rotation, 6 numbers), and the 3 unmeasured phases.
Solver: ULYSSES 3.0.2, EtaB_3DME (density-matrix equations, three partners, thermal start as the 4x4 model predicts; an earlier test).
Install: see ulysses_check_2026-11-14.py.
Usage: python3 matter_excess_rule5_exact_2026-11-14.py            (quick parts + summary of the saved scans)
       python3 matter_excess_rule5_exact_2026-11-14.py search 3.2 2  (rerun a directed search: max tilt in radians, seed; about 1 hour)
Credit: Fukugita & Yanagida 1986; Casas & Ibarra 2001; Akhmedov, Frigerio & Smirnov 2003; Di Bari & Riotto 2009; Blanchet, Di Bari, Jones & Marzola 2013;
Granelli, Moffat, Perez-Gonzalez, Schulz & Turner (ULYSSES) 2021/2023; Harrison, Perkins & Scott (tribimaximal); Cabibbo; Kobayashi & Maskawa;\nPati & Salam; NuFIT 6; Planck."""
import math, sys, json, os, itertools, warnings
warnings.filterwarnings("ignore")
import numpy as np, ulysses
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
D = np.diag([0.51099895e-3, 0.1056583755, 1.77693])      # GeV
MNU = np.array([0.00036, 0.0087, 0.0501])*1e-9             # GeV (record M2 / an earlier test)
V, OBS = 174.0, 6.1e-10
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
def build(VL, d, a21, a31, s23):
    U = pmns(0.307, 0.02195, s23, d, a21, a31); mnu = U @ np.diag(MNU) @ U.T
    MR = -D @ VL.conj() @ np.linalg.inv(mnu) @ VL.conj().T @ D
    S, W = takagi(MR); o = np.argsort(S); S = S[o]; W = W[:, o]
    return S, VL.conj().T @ D @ W.conj()/V
def eta(VL, d, a21, a31, s23, init=1, model="EtaB_3DME", zsteps=1500, full=False):
    S, h = build(VL, d, a21, a31, s23)
    L = getattr(ulysses, model)(zmin=min(1e-3, S[0]/(30*S[2])), zmax=1000, zsteps=zsteps, ordering=0, initial_abundance=init, which_param="manual")
    pd = dict(m=math.log10(MNU[0]*1e9), M1=math.log10(S[0]), M2=math.log10(S[1]), M3=math.log10(S[2]))
    for i in range(3):
        for j in range(3): pd[f"Y{i+1}{j+1}_mag"] = abs(h[i, j]); pd[f"Y{i+1}{j+1}_phs"] = float(np.angle(h[i, j]))
    try: e = float(L(pd))
    except Exception: e = 0.0
    if not np.isfinite(e): e = 0.0
    return (e, S, h, L) if full else e
def eta_x(x, **kw): return eta(expm(1j*herm(x[:9])), x[9], x[10], x[11], x[12], **kw)
def scan_phases(VL, n, rng):
    vals, best = [], (0, None, None, None)
    for _ in range(n):
        d, a, b = rng.uniform(0, 2*math.pi, 3); s23 = float(rng.uniform(0.435, 0.585))
        e, S, h, L = eta(VL, d, a, b, s23, full=True); vals.append(abs(e))
        if abs(e) >= abs(best[0]): best = (e, S, (L.k1.real, L.k2.real, L.k3.real), (d, a, b, s23))
    return np.array(vals), best
def search(scale, seed):
    from scipy.optimize import differential_evolution
    f = lambda p: -abs(eta_x(p))
    res = differential_evolution(f, [(-scale, scale)]*9 + [(0, 2*math.pi)]*3 + [(0.435, 0.585)], seed=seed, popsize=10, maxiter=45, polish=False, init="sobol",
                                 callback=lambda xk, convergence=None: print(f"  best so far {-f(xk):.3e}", flush=True))
    print("best", -res.fun, list(res.x))
def pr(*a): print(*a, flush=True)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "search": search(float(sys.argv[2]), int(sys.argv[3])); sys.exit()
    rng = np.random.default_rng(33); I3 = np.eye(3, dtype=complex)
    pr("="*96); pr("A. Self-check: the seesaw returns the light masses, and the Dirac links are exactly the lepton masses"); pr("="*96)
    VLr = expm(1j*herm(rng.normal(size=9))); S, h = build(VLr, 1.0, 2.0, 3.0, 0.47)
    mnu = -(V*h) @ np.diag(1/S) @ (V*h).T
    pr("  light masses back (eV):", np.round(np.sort(np.sqrt(np.linalg.eigvalsh(mnu @ mnu.conj().T)))*1e9, 6).tolist(), "| Dirac links (GeV):", np.round(np.sort(np.linalg.svd(V*h, compute_uv=False)), 6).tolist())

    pr("\n"+"="*96); pr("B. Parameter-free perspectives (only the 3 unmeasured phases and theta23 scanned). Observed: 6.1e-10"); pr("="*96)
    names = ["e", "mu", "tau"]
    for perm in itertools.permutations(range(3)):
        P = np.zeros((3, 3), dtype=complex)
        for i, j in enumerate(perm): P[i, j] = 1
        vals, best = scan_phases(P.T, 120, rng)
        lab = "links: " + ", ".join(f"nu_{names[i]}~m_{names[j]}" for i, j in enumerate(perm)) + ("   <- rule 5 read literally" if perm == (0, 1, 2) else ("   <- conjugate side (twist -2/9 swaps e and mu)" if perm == (1, 0, 2) else ""))
        pr(f"  {lab:78s} best {abs(best[0]):.1e}  typical {np.median(vals):.0e}  partners {', '.join('%.0e' % s for s in best[1])} GeV")
    w = np.exp(2j*math.pi/3); F = np.array([[1, 1, 1], [1, w, w*w], [1, w*w, w]])/math.sqrt(3)
    TBM = np.array([[math.sqrt(2/3), 1/math.sqrt(3), 0], [-1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)], [1/math.sqrt(6), -1/math.sqrt(3), 1/math.sqrt(2)]], dtype=complex)
    def ckm(s12=0.22501, s23=0.04183, s13=0.003732, d=1.147):      # PDG 2024
        c12, c23, c13 = (math.sqrt(1-x*x) for x in (s12, s23, s13)); e = np.exp(1j*d)
        return np.array([[c12*c13, s12*c13, s13/e], [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, s23*c13], [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13]])
    for name, VL in [("third-of-a-turn (Fourier) rotation", F), ("its inverse", F.conj().T), ("tribimaximal rotation", TBM), ("its inverse", TBM.T),
                     ("the quarks' own mixing (CKM), as the 4x4 would suggest", ckm()), ("its inverse", ckm().conj().T)]:
        vals, best = scan_phases(VL, 120, rng)
        pr(f"  {name:78s} best {abs(best[0]):.1e}  typical {np.median(vals):.0e}  partners {', '.join('%.0e' % s for s in best[1])} GeV")
    pr("  -> every parameter-free choice tried is 2,000x or more short.")

    pr("\n"+"="*96); pr("C. Random perspectives (uniform over all rotations), random phases"); pr("="*96)
    f = os.path.join(HERE, "matter_excess_random_perspectives_2026-11-14.json")
    if os.path.exists(f):
        e = np.abs(np.array(json.load(open(f))["eta"]))
        pr(f"  {len(e)} random perspectives: typical (median) {np.median(e):.1e}; top 10% above {np.percentile(e, 90):.1e}; top 1% above {np.percentile(e, 99):.1e}; best {e.max():.1e}")
        pr(f"  fraction reaching the observed 6.1e-10: {(e >= OBS).mean()*100:.2f}%   fraction within 10x of it: {(e >= OBS/10).mean()*100:.2f}%")
    else: pr("  (saved scan not found)")

    pr("\n"+"="*96); pr("D. Directed search (differential evolution, then local polish): the largest excess found with rule 5 exact"); pr("="*96)
    f = os.path.join(HERE, "matter_excess_best_points_2026-11-14.json")
    if os.path.exists(f):
        for pt in json.load(open(f)):
            x = np.array(pt["x"]); e, S, h, L = eta_x(x, full=True); VL = expm(1j*herm(x[:9])); hh = (h.conj().T @ h).real
            pr(f"\n  [{pt['label']}]")
            pr(f"   eta_B = {abs(e):.2e} ({abs(e)/OBS:.2f} x observed); with scatterings {abs(eta_x(x, model='EtaB_3DMEsct')):.2e}; zero start {abs(eta_x(x, init=0)):.2e}")
            pr(f"   partners: {S[0]:.2e}, {S[1]:.2e}, {S[2]:.2e} GeV (right-side scale 10^11.7 = 5.0e11); washout strengths K = {L.k1.real:.1f}, {L.k2.real:.1f}, {L.k3.real:.1f}")
            pr(f"   perspective |V_L| rows: {np.round(np.abs(VL), 2).tolist()}; largest tilt from 'literal' = {math.degrees(math.asin(min(1, np.abs(VL - np.diag(np.diag(VL))).max()))):.0f} deg")
            pr(f"   lightest partner's link shares (e, mu, tau): {np.round(np.abs(h[:, 0])**2/hh[0, 0], 3).tolist()}   phases delta, a21, a31 = {np.round(np.degrees(x[9:12]) % 360).tolist()} deg, sin^2(th23) = {x[12]:.3f}")
            pr(f"   numerical stability: 500 steps {abs(eta_x(x, zsteps=500)):.4e}; 6000 steps {abs(eta_x(x, zsteps=6000)):.4e}")
            r2 = np.random.default_rng(0)
            for pc in (0.01, 0.05):
                nv = np.array([abs(eta_x(x*(1 + pc*r2.normal(size=13)))) for _ in range(12)])
                pr(f"   nudging every setting by {int(pc*100)}%: stays between {nv.min():.1e} and {nv.max():.1e} (median {np.median(nv):.1e})")
    pr("\nReading: with rule 5 exact the observed excess CAN be reached, but only in a narrow corner (a 5% nudge loses most of it),")
    pr("and no parameter-free perspective works. It is allowed, not predicted. The missing piece is a rule for the perspective V_L.")
