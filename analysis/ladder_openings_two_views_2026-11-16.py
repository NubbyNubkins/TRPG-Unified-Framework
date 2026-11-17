"""Are the two quark ladder openings the one quark pair seen from two perspectives?"""
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
def Q(m): m=np.abs(m); return m.sum()/np.sqrt(m).sum()**2
def opening(m): e=math.sqrt(2*(3*Q(m)-1)); return math.degrees(math.atan(e/math.sqrt(2)))
pr("openings (top footing): up %.3f  down %.3f  lepton %.3f"%(opening(mu),opening(md),opening(np.array([0.51099895,105.6583755,1776.93]))))
Qu=Q(mu)
# (1) pure change of facing: rotate the down matrix by a random unitary -> masses unchanged
rng=np.random.default_rng(1); A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)); U,_=np.linalg.qr(A)
H=np.diag(md**2); H2=U@H@U.conj().T; pr("rotation check: masses before",md," after",np.sqrt(np.linalg.eigvalsh(H2)))
# (2) quark pair seen two ways, press on the down side only: unpressed down base forced to the up ladder's opening (and optionally its full shape)
def fitc(mode, su, n=40, seed=3):
    r_=np.random.default_rng(seed); Uu=ups[su][0]; best=None
    def parts(x):
        kap=1/(1+math.exp(-x[3])); u=mk_u("hier",QT,x[4],x[5]); return kap,u,Hd(x[:3],kap,u)
    def res(x,wt):
        kap,u,H=parts(x); w=np.sqrt(np.abs(np.linalg.eigvalsh(H))); o=obs(ckm(Uu,H)); eb=np.sort(np.abs(np.linalg.eigvalsh(base(x[:3]))))
        r=[wt*(np.log(w)-np.log(md)), (o[0]-Vm)/Vs]
        if mode in ("opening","shape"): r.append([wt*(Q(eb)-Qu)/0.002])
        if mode=="shape": r.append(wt*(np.log(eb[:2]/eb[2])-np.log(mu[:2]/mu[2]))/0.03)
        return np.concatenate(r)
    for sd in SG:
        p0=fit_base(md,sd)
        if p0 is None: continue
        for t in range(n):
            x0=[*p0*(1+0.5*r_.normal(size=3)), r_.uniform(2,8), r_.uniform(0.9,1.5)*r_.choice([-1,1]), r_.uniform(0,2*math.pi)]
            try: r=least_squares(res,x0,args=(30.0,),max_nfev=400); r=least_squares(res,r.x,args=(1e3,),max_nfev=400)
            except Exception: continue
            kap,u,H=parts(r.x); w=np.sqrt(np.abs(np.linalg.eigvalsh(H))); o=obs(ckm(Uu,H)); eb=np.sort(np.abs(np.linalg.eigvalsh(base(r.x[:3]))))
            chi_t=float(np.sum(((o[0]-Vm)/Vs)**2)); chi_m=float(np.sum(((np.log(w)-np.log(md))/rel_d)**2))
            dev_open=opening(eb)-opening(mu)
            c=chi_t+chi_m+(0 if mode=="free" else (dev_open/0.018)**2)+(0 if mode!="shape" else float(np.sum(((np.log(eb[:2]/eb[2])-np.log(mu[:2]/mu[2]))/rel_u[:2])**2)))
            if best is None or c<best[0]: best=(c,chi_t,chi_m,dev_open,1-kap,o,eb)
    return best
for mode in ("free","opening","shape"):
    for su in SG:
        b=fitc(mode,su)
        if b is None: pr(mode,su,"no fit"); continue
        c,ct,cm,dv,omk,o,eb=b
        pr(f"{mode:8s} up-signs {su}: total chi2 {c:9.2f} (tilts {ct:8.2f}, masses {cm:8.2f}, opening dev {dv:+.3f} deg), 1-kappa {omk:.2e}, tilts {o[0]}, unpressed down base {eb}")
def Q(m): m=np.abs(m); return m.sum(-1)/np.sqrt(m).sum(-1)**2
def th(q): return np.degrees(np.arctan(np.sqrt(2*(3*q-1))/math.sqrt(2)))
rng=np.random.default_rng(0); N=200000
for lab,uu,dd in (("top footing",mu,md),("Z footing (top 168.3 GeV)",np.array([1.23,620.,168300.]),np.array([2.67,53.2,2839.]))):
    U=uu*(1+rel_u*rng.normal(size=(N,3))); D=dd*(1+rel_d*rng.normal(size=(N,3)))
    qu,qd=Q(U),Q(D); tu,td=th(qu),th(qd)
    tests={"R1 Q_down x Q_up = 2/3 (two views multiply to the lepton view)":(qu*qd,2/3),
           "R2 tan(down) x tan(up) = sqrt2 (two views multiply to the cube edge)":(np.tan(np.radians(tu))*np.tan(np.radians(td)),math.sqrt(2)),
           "R3 mean opening = midpoint of 45 and 54.74 (49.87)":((tu+td)/2,(45+math.degrees(math.atan(math.sqrt(2))))/2)}
    print(lab, f"up {tu.mean():.3f}+-{tu.std():.3f}  down {td.mean():.3f}+-{td.std():.3f}")
    for k,(x,t) in tests.items(): print(f"   {k}: {x.mean():.5f} +- {x.std():.5f} vs {t:.5f}  -> {(x.mean()-t)/x.std():+.1f} sigma")
