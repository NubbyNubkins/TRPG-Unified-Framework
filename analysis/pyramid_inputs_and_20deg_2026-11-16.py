"""Is anything else possible for the pyramid near miss? Tilt inputs (global fit vs direct, inclusive vs exclusive) and the post-hoc 20-degree (pi/9) lean."""
import math, numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
np.set_printoptions(precision=4, suppress=True)
J1 = np.ones((3, 3)); mu = np.array([1.242, 627.8, 162500.0]); md = np.array([2.678, 53.36, 2873.0])
rel_u = np.array([0.034, 0.010, 0.006]); rel_d = np.array([0.018, 0.014, 0.005])     # PDG 2024 errors (+1% running)
Vm = np.array([0.2250, 0.04183, 0.003732, 0.00857]); Vs = np.array([0.00068, 0.00074, 0.000088, 0.00018])   # |Vus| |Vcb| |Vub| |Vtd|
F = np.array([[1, 1, 1], [-1, 1, 1], [0, -2, 1]], float)/np.array([math.sqrt(2), math.sqrt(6), math.sqrt(3)])  # hierarchical basis
SG = [(1, -1), (-1, 1)]; TOP = "top's energy"
base = lambda p: p[0]*(J1 + np.diag([-p[1], p[1], p[2]]))
def fit_base(m, signs, start=None):
    tgt = np.sort(np.array(m)*np.array([*signs, 1])); f = lambda p: np.sort(np.linalg.eigvalsh(base(p))) - tgt
    starts = [start] if start is not None else [[m[2]/3, d0, e0] for d0 in (0.0001, 0.001, 0.01, 0.05) for e0 in (-0.5, -0.2, 0.01, 0.2, 0.5)]
    for s in starts:
        p, info, ok, _ = fsolve(f, s, full_output=True)
        if ok == 1 and np.max(np.abs(f(p))) < 1e-7*m[2]: return p
    return None
def Hd(p, kap, u):
    m = base(p); return m @ (np.eye(3) - kap*np.outer(u, u.conj())) @ m
def obs(V):
    A = np.abs(V); t = np.array([A[0,1], A[1,2], A[0,2], A[2,0]])
    Jc = (V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1])).imag
    al = np.degrees(np.angle(-V[2,0]*np.conj(V[2,2])/(V[0,0]*np.conj(V[0,2]))))
    be = np.degrees(np.angle(-V[1,0]*np.conj(V[1,2])/(V[2,0]*np.conj(V[2,2]))))
    ga = np.degrees(np.angle(-V[0,0]*np.conj(V[0,2])/(V[1,0]*np.conj(V[1,2]))))
    return t, Jc, al, be, ga
def ckm(Uu, H): return Uu.conj().T @ np.linalg.eigh(H)[1]
def up_U(signs, m=mu):
    p = fit_base(m, signs); w, O = np.linalg.eigh(base(p)); return O[:, np.argsort(np.abs(w))], p
worst = lambda t, V=Vm: math.exp(max(abs(np.log(t/V))))

# one-loop running from m_t to the far-side scale (Machacek & Vaughn 1984; Arason and others 1992)
def rge(t, y):
    g1, g2, g3, yt = y; k = 1/(16*math.pi**2)
    return [k*41/10*g1**3, -k*19/6*g2**3, -k*7*g3**3, k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - 0.85*g1**2), k*yt**2]
y0 = [math.sqrt(4*math.pi/58.9), math.sqrt(4*math.pi/29.6), math.sqrt(4*math.pi*0.108), 162.5/174.1]
tfar = math.log(2.1e8/162.5); sol = solve_ivp(lambda t, y: rge(t, y[:4]), [0, tfar], y0 + [0.0], rtol=1e-9, atol=1e-12)
I_far = sol.y[4, -1]; ef = math.exp(1.5*I_far)
mu_far = mu*np.array([1, 1, ef]); md_far = md*np.array([1, 1, 1/ef]); Vm_far = Vm*np.array([1, 1/ef, 1/ef, 1/ef]); Vs_far = Vs*np.array([1, 1/ef, 1/ef, 1/ef])

shape = {}
for name, uu, dd in (("top's energy", mu, md), ("far-side 2.1e8 GeV", mu_far, md_far)):
    for s in SG:
        pu, pd = fit_base(uu, s), fit_base(dd, s); shape[(name, s)] = (pu, pd)
ups = {s: up_U(s) for s in SG}
def mk_u(basis, ph, a, b):
    s = np.array([math.cos(a)*math.cos(b)*np.exp(1j*ph[0]), math.cos(a)*math.sin(b)*np.exp(1j*ph[1]), math.sin(a)*np.exp(1j*ph[2])])
    return F @ s if basis == "hier" else s
def fit_nb(basis, ph, su, n, seed, kapfix=None, mdn=md, Uu=None, Vt=Vm, Vsg=Vs, allmin=False, dirmap=None):
    """3 free numbers (kappa, 2 angles) + the down base (3, pinned by the 3 down masses). Staged least squares from n random starts
    for each of the two down-base sign branches. Returns (chi2 on the four tilts, x, (kappa, u, H_d))."""
    r_ = np.random.default_rng(seed); Uu = ups[su][0] if Uu is None else Uu; best = None; mins = []
    def parts(x):
        kap = kapfix if kapfix is not None else 1/(1 + math.exp(-x[3])); a_, b_ = (x[4], x[5]) if dirmap is None else dirmap(x[4]); u = mk_u(basis, ph, a_, b_); return kap, u, Hd(x[:3], kap, u)
    def res(x, wt):
        kap, u, H = parts(x); w = np.sqrt(np.abs(np.linalg.eigvalsh(H))); o = obs(ckm(Uu, H))
        return np.concatenate([wt*(np.log(w) - np.log(mdn)), (o[0] - Vt)/Vsg])
    for sd in SG:
        p0 = fit_base(mdn, sd)
        for trial in range(n):
            x0 = [*p0*(1 + 0.3*r_.normal(size=3)), r_.uniform(2, 8), r_.uniform(0.9, 1.5)*r_.choice([-1, 1]), r_.uniform(0, 2*math.pi)]
            try: r = least_squares(res, x0, args=(30.0,), max_nfev=300); r = least_squares(res, r.x, args=(1e3,), max_nfev=300)
            except Exception: continue
            c = float(np.sum(r.fun**2)); mins.append(round(c, 1))
            if best is None or c < best[0]: best = (c, r.x, parts(r.x))
    return (best, sorted(set(mins))) if allmin else best
QT = (math.pi/2, 0.0, 0.0)
def run(step, Vt, Vsg, n=12, seed=900):
    s=np.array([step/3,1/3,1]); s/=np.linalg.norm(s); a=math.asin(s[2]); b=math.atan2(s[1],s[0]); best=None
    for su in SG:
        r=fit_nb("hier",QT,su,n,seed,Vt=Vt,Vsg=Vsg,dirmap=lambda x:(a,b))
        if best is None or r[0]<best[0]: best=(r[0],su,r)
    c,su,r=best; t,J,al,be,ga=obs(ckm(ups[su][0],r[2][2])); return c,t,abs(be),abs(ga),abs(J)
P=1/(2*math.sqrt(2))
sets=[("global-fit tilt values used so far (CKMfitter-style; these already lean on beta)",Vm,Vs),
 ("direct measurements only (PDG 2024: Vus 0.2243(8), Vcb 41.1(12)e-3, Vub 3.82(20)e-3, Vtd 8.6(2)e-3)",np.array([0.2243,0.0411,0.00382,0.0086]),np.array([0.0008,0.0012,0.00020,0.0002])),
 ("exclusive B decays (Vcb 39.8(6)e-3, Vub 3.70(16)e-3)",np.array([0.2243,0.0398,0.00370,0.0086]),np.array([0.0008,0.0006,0.00016,0.0002])),
 ("inclusive B decays (Vcb 42.2(5)e-3, Vub 4.13(26)e-3)",np.array([0.2243,0.0422,0.00413,0.0086]),np.array([0.0008,0.0005,0.00026,0.0002]))]
print("Pyramid (light step 1/(2 sqrt2)): predicted beta against measured 22.6 +- 0.4 deg, for each choice of tilt inputs")
for lab,Vt,Vsg in sets:
    c,t,be,ga,J=run(P,Vt,Vsg); print(f"  {lab}\n     tilts chi2 {c:5.2f} | beta {be:.2f} ({(be-22.6)/0.42:+.1f} sigma) | gamma {ga:.1f} ({(ga-65.9)/3.4:+.1f} sigma) | J {J:.2e}")
print("Post hoc (3-6-9 angle, found after seeing the band): face lean 20 deg = pi/9, step tan(20) =",round(math.tan(math.pi/9),4))
for lab,Vt,Vsg in sets[:2]:
    c,t,be,ga,J=run(math.tan(math.pi/9),Vt,Vsg); print(f"  {lab[:40]:40s}: tilts chi2 {c:5.2f} | beta {be:.2f} ({(be-22.6)/0.42:+.1f} sigma) | gamma {ga:.1f} | J {J:.2e}")
