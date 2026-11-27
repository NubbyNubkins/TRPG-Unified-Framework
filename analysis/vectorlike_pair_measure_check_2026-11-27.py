"""Lemma A check: on a 4D lattice, the left-handed overlap projector in the
conjugate representation R* is the image of the one in R under one fixed
antiunitary map, v -> A conj(v) with A = gamma5 C^{-1}.  Then a left-handed pair
R + R* (TRPG's far-side 5 + 5bar) has a canonical, phase-free, gauge-invariant
measure, and its return phases cancel exactly: T[R*] = conj(T[R]).
Control: in 2D no such map exists (there L(q)+L(-q) is anomalous), so the
cancellation must FAIL there.  Overlap set-up as in Luscher (1999) and
Craig (2026): H = g5 (D_W - 1), P_minus = (1 + sign H)/2.
"""
import numpy as np, itertools
rng = np.random.default_rng(20261127)
s1 = np.array([[0,1],[1,0]],complex); s2 = np.array([[0,-1j],[1j,0]]); s3 = np.diag([1,-1]).astype(complex)
I2 = np.eye(2); Z2 = np.zeros((2,2))

def gammas(d):
    if d == 2:
        return [s1, s2], s3
    g = [np.block([[Z2, -1j*s], [1j*s, Z2]]) for s in (s1, s2, s3)]
    g.append(np.block([[Z2, I2], [I2, Z2]]))
    return g, np.diag([1,1,-1,-1]).astype(complex)

def su_n(n, eps):
    X = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)); X = (X + X.conj().T)/2
    X -= np.trace(X)/n*np.eye(n)
    w, v = np.linalg.eigh(X); return v @ np.diag(np.exp(1j*eps*w)) @ v.conj().T, X

def expm_h(X, t):
    w, v = np.linalg.eigh(X); return v @ np.diag(np.exp(1j*t*w)) @ v.conj().T

def wilson(U, L, d, gam, nc):
    sites = list(itertools.product(range(L), repeat=d)); idx = {s:i for i,s in enumerate(sites)}
    ns = len(gam[0]); N = len(sites)*ns*nc
    D = np.zeros((N,N), complex)
    blk = lambda i: slice(i*ns*nc, (i+1)*ns*nc)
    for s in sites:
        i = idx[s]; D[blk(i), blk(i)] += d*np.eye(ns*nc)
        for mu in range(d):
            sp = list(s); sp[mu] = (sp[mu]+1) % L; j = idx[tuple(sp)]
            u = U[mu][i]
            D[blk(i), blk(j)] += -0.5*np.kron(np.eye(ns) - gam[mu], u)
            D[blk(j), blk(i)] += -0.5*np.kron(np.eye(ns) + gam[mu], u.conj().T)
    return D

def pminus(U, L, d, nc):
    gam, g5 = gammas(d)
    D = wilson(U, L, d, gam, nc)
    G5 = np.kron(np.eye(L**d), np.kron(g5, np.eye(nc)))
    H = G5 @ (D - np.eye(len(D))); H = (H + H.conj().T)/2
    w, v = np.linalg.eigh(H)
    return (np.eye(len(D)) + v @ np.diag(np.sign(w)) @ v.conj().T)/2, np.min(np.abs(w))

def return_phase(Ufun, L, d, nc, steps):
    P0, gap = pminus(Ufun(0.0), L, d, nc)
    w, v = np.linalg.eigh(P0); V0 = v[:, w > 0.5]; V = V0.copy(); mingap = gap
    for k in range(1, steps+1):
        P, gap = pminus(Ufun(k/steps), L, d, nc); mingap = min(mingap, gap)
        W = P @ V                                    # parallel transport: polar part only,
        u, _, vh = np.linalg.svd(W, full_matrices=False)   # so no spurious phase or sign enters
        V = u @ vh
    T = np.linalg.det(V0.conj().T @ V)
    return T/abs(T), mingap

def loop(L, d, nc, eps0, eps1):
    nsite = L**d
    base = [[su_n(nc, eps0)[0] for _ in range(nsite)] for _ in range(d)]
    X = [[su_n(nc, 1)[1] for _ in range(nsite)] for _ in range(d)]
    Y = [[su_n(nc, 1)[1] for _ in range(nsite)] for _ in range(d)]
    def U(t, conj=False):
        c, s = np.cos(2*np.pi*t), np.sin(2*np.pi*t)
        out = [[expm_h(c*X[m][i] + s*Y[m][i], eps1) @ base[m][i] for i in range(nsite)] for m in range(d)]
        return [[u.conj() for u in row] for row in out] if conj else out
    return U

out = []; p = out.append
# --- exact identity in 4D -------------------------------------------------
gam, g5 = gammas(4)
Cc = gam[1] @ gam[3]                                   # C with C g C^-1 = -g^T
ok = all(np.allclose(Cc @ g @ np.linalg.inv(Cc), -g.T) for g in gam)
A = g5 @ np.linalg.inv(Cc)
p(f"4D: C gamma C^-1 = -gamma^T for all mu: {ok}")
p(f"4D: A conj(gamma_mu) A^-1 = gamma_mu: {all(np.allclose(A @ g.conj() @ np.linalg.inv(A), g) for g in gam)};"
  f"  A conj(g5) A^-1 = +g5 (keeps handedness): {np.allclose(A @ g5.conj() @ np.linalg.inv(A), g5)}")
L, d, nc = 3, 4, 3
U = loop(L, d, nc, 0.25, 0.15)
P, gap = pminus(U(0.3), L, d, nc); Pc, _ = pminus(U(0.3, conj=True), L, d, nc)
AA = np.kron(np.eye(L**d), np.kron(A, np.eye(nc)))
p(f"4D SU(3), 3^4 lattice: || P-[U*] - A conj(P-[U]) A^-1 || = {np.linalg.norm(Pc - AA @ P.conj() @ np.linalg.inv(AA)):.2e}"
  f"  (|| P-[U] || = {np.linalg.norm(P):.1f}; H gap {gap:.3f})")
p(f"   rank P-[U] = {round(np.trace(P).real)}, rank P-[U*] = {round(np.trace(Pc).real)}")
# --- return phases along a closed loop -----------------------------------
for steps in (24,):
    T, g1 = return_phase(lambda t: U(t), L, d, nc, steps)
    Tc, g2 = return_phase(lambda t: U(t, True), L, d, nc, steps)
    p(f"4D loop ({steps} steps, min H gap {min(g1,g2):.3f}): T[R] = {np.angle(T):+.6f} rad, "
      f"T[R*] = {np.angle(Tc):+.6f} rad, T[R]T[R*] phase = {np.angle(T*Tc):+.2e}")
# --- 2D control: must fail ------------------------------------------------
L2 = 6
U2 = loop(L2, 2, 3, 0.3, 0.6)
T2, g = return_phase(lambda t: U2(t), L2, 2, 3, 40)
T2c, g = return_phase(lambda t: U2(t, True), L2, 2, 3, 40)
p(f"2D control, SU(3) on 6^2: T[R] = {np.angle(T2):+.6f}, T[R*] = {np.angle(T2c):+.6f}, "
  f"product phase = {np.angle(T2*T2c):+.4f}  (nonzero expected: no handedness-keeping map in 2D)")
print("\n".join(out))
