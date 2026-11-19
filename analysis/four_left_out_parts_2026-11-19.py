"""Relationships, stage, quantum rules, starting conditions (19 Nov 2026). Exploration only."""
import math, numpy as np
pi=math.pi; rng=np.random.default_rng(7)
print("="*72); print("1. RELATIONSHIPS: can the link strengths be reduced?")
print("  Link strengths in the Standard Model: 6 quark + 3 lepton weights, 4 quark mixings = 13; neutrinos add 3 weights + 3-6 mixings.")
# Gatto-Sartori-Tonin: ratio md/ms runs identically, so scale-consistent
md_ms=1/19.9
print(f"  Cabibbo from down/strange: sqrt(md/ms) = {math.sqrt(md_ms):.4f} vs |Vus| = 0.2243 +/- 0.0008 -> HOLDS (scale-independent; TRPG already uses it)")
print(f"  |Vcb| from up/charm, mixed scales (mu@2GeV, mc@mc): sqrt(2.16/1273) = {math.sqrt(2.16/1273):.4f} vs 0.0411 +/- 0.0012 (looks perfect)")
mc2=1.273*(0.301/0.391)**(12/25)
print(f"     same scale (both at 2 GeV, mc(2GeV) = {mc2:.3f} GeV): {math.sqrt(2.16e-3/mc2):.4f} -> {(math.sqrt(2.16e-3/mc2)-0.0411)/0.0012:+.1f} sigma. FAILS when done consistently.")
print(f"  |Vub| from up/top (both at MZ: 1.27 MeV, 169 GeV): sqrt = {math.sqrt(1.27e-3/169):.5f} vs 0.00382 +/- 0.00020 -> {(math.sqrt(1.27e-3/169)-0.00382)/0.0002:+.0f} sigma. FAILS.")
print(f"  Fritzsch pattern |Vcb| = |sqrt(ms/mb) - sqrt(mc/mt)| = {abs(math.sqrt(93.5/4183)-math.sqrt(1273/172570)):.4f} vs 0.0411 -> FAILS (known, ruled out ~2000s)")
print("  Lepton weights: Koide/45 deg rule fixes tau from e, mu (TRPG): 3 -> 2. Neutrinos: TM1 + 3-2sqrt2 (TRPG leads).")
print("  Net: 13 quark/lepton strengths -> about 10 honest (tau, down/Cabibbo, one neutrino) ; the rest stay measured.")
print("  Charges as relationships: each charge = which links touch that house (Noether 1918). 0 extra numbers. HOLDS.")

print("="*72); print("2. THE STAGE: can space and time come out of relationships?")
# 2x2 Hermitian matrix <-> spacetime point; det = interval; SL(2,C) = Lorentz
s=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.array([[1,0],[0,-1]])]
x=rng.normal(size=4); X=sum(x[i]*s[i] for i in range(4))
A=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); A/=np.sqrt(np.linalg.det(A))
Y=A@X@A.conj().T
print(f"  A 2x2 matrix of complex relationships: its determinant = t^2-x^2-y^2-z^2 = {np.linalg.det(X).real:.6f} (direct: {x[0]**2-x[1]**2-x[2]**2-x[3]**2:.6f})")
print(f"  after an arbitrary 'change of perspective' A: determinant {np.linalg.det(Y).real:.6f} -> unchanged = special relativity, 3 space + 1 time")
print("  -> 3+1 dimensions and light-speed geometry ARE the geometry of a 2x2 complex relationship matrix (Penrose spinors). HOLDS.")
print("     (octonion 2x2 gives 9+1, the superstring stage; quaternion 2x2 gives 5+1) -> complex numbers pick 3+1.")
# Ehrenfest: stable orbits only for 3 space dims
for n in (2,3,4,5):
    p=n-1  # force ~ 1/r^p; circular orbit stable iff p<3
    print(f"  {n} space dims: pull ~ 1/r^{p}: orbits {'stable' if p<3 else 'UNSTABLE (fall in or fly off)'}{'; but no inverse-square, atoms differ' if n==2 else ''}")
print("  -> only 3 space dimensions give both stable planets and stable atoms (Ehrenfest 1917). HOLDS as a selection rule.")
print("  Gravity from relationships: Einstein's equation follows from heat flow + 'entanglement entropy ~ area' (Jacobson 1995). HOLDS (derivation).")
print("  Space from entanglement: cutting entanglement pulls space apart (Van Raamsdonk 2010; Ryu-Takayanagi). HOLDS in model universes; not yet shown for ours.")
print("  Time from relationships: a still universe still shows change to any part entangled with a clock (Page-Wootters 1983; Moreva et al. 2014 lab). HOLDS.")
print("  Grain of space: Fermi GRB 090510 photons arrive together -> no light-speed change below ~7x Planck energy (Fermi-LAT 2009). Any grain must hide from light.")
print("  Units: c, hbar, G are conversion factors (SI 2019 fixes c, h exactly). Only pure ratios are real numbers: m_p/M_Planck = 7.7e-20.")

print("="*72); print("3. QUANTUM RULES: which follow from which?")
x=np.linspace(-40,40,2**14); dx=x[1]-x[0]
for name,f in [("bell curve",np.exp(-x**2/2)),("square pulse",(abs(x)<1.0)*1.0),("two-peak wave",np.exp(-(x-3)**2)+np.exp(-(x+3)**2))]:
    f=f/np.sqrt((abs(f)**2).sum()*dx); px=(abs(f)**2)*dx
    sx=np.sqrt((px*x**2).sum()-(px*x).sum()**2)
    F=np.fft.fftshift(np.fft.fft(f)); k=np.fft.fftshift(np.fft.fftfreq(len(x),dx))*2*pi; pk=abs(F)**2; pk/=pk.sum()
    sk=np.sqrt((pk*k**2).sum()-(pk*k).sum()**2)
    print(f"  {name:13s}: spread in place x spread in wave-number = {sx*sk:.3f}  (never below 0.5)")
print("  -> Uncertainty is not a separate rule; it is what any wave does (Fourier). Folds into 'everything is a wave'. HOLDS.")
print(f"  Bell/CHSH: no-dice local limit 2, quantum 2*sqrt2 = {2*math.sqrt(2):.3f}; loophole-free 2.42 +/- 0.20 (Hensen 2015), many later at >10 sigma.")
print("     -> local hidden 'no dice' FAILS; non-local 'touching through geometry' (Bohm 1952 type) ALLOWED, as TRPG says.")
print("  Born rule (chance = amplitude^2): forced by the state space itself in 3+ dimensions (Gleason 1957). Three-path test of extra rules: kappa ~ 0 to 1e-2 (Sinha 2010). HOLDS.")
print("  Complex numbers needed: real-number quantum theory ruled out experimentally (Renou et al. 2021; Chen et al. 2022). Matches TRPG rule 4 (one i = one facing). HOLDS.")
print("  Spin-statistics (Pauli 1940) and CPT (Luders-Pauli) are theorems of stage + quantum + locality. Fold away. HOLDS.")
print("  Entanglement-first derivation: quantum theory derived from information rules incl. 'every mixed state is part of an entangled whole'")
print("     (Chiribella, D'Ariano, Perinotti 2011; Hardy 2001). HOLDS as a published derivation.")
print("  Gravity-caused collapse (Diosi-Penrose, simplest form): ruled out by Gran Sasso radiation search (Donadi et al. 2021). FAILS.")

print("="*72); print("4. STARTING CONDITIONS: how few?")
print("  List: matter excess 6.1e-10; ripple size A_s 2.1e-9; ripple tilt n_s 0.9649+/-0.0042; flatness |Omega_k|<0.002;")
print("        same temperature everywhere (1 in 1e5); low-entropy start; dark-matter amount; dark-energy amount; light species 3.")
for N in (55,57,60):
    ns=1-2/N; r=12/N**2
    print(f"  One field with a flat top (Starobinsky 1980 / Higgs as the field, Bezrukov-Shaposhnikov 2008), N={N}: n_s = {ns:.4f}, r = {r:.4f}")
print("     measured n_s 0.9649+/-0.0042, r < 0.036 -> HOLDS. Pending: LiteBIRD/CMB-S4 will reach r ~ 0.001-0.003 (not 100% until then).")
print("     Using the HIGGS (the resistance) as the starter folds 3 start facts (flat, same temperature, ripples) into 1 existing field + 1 number.")
print("     Catch: with Standard-Model values the Higgs self-pull turns negative near 1e10-1e11 GeV, which blocks Higgs-start; TRPG's far-side")
print("     states change that running (not checked here). Starobinsky (space's own stiffness, R^2) has no such catch.")
S_now=3e104; S_max=1e123
print(f"  Low-entropy start: entropy now ~{S_now:.0e} (mostly black holes) vs ~{S_max:.0e} possible -> chance ~ 10^-(10^123) (Penrose). One condition: 'smooth start' (Weyl curvature = 0).")
print("     Arrow of time = this one condition. Folds 'time direction' into the start. HOLDS as standard reasoning.")
print("  Matter excess from hidden-partner decays (leptogenesis): TRPG lead (~2e-10 vs 6.1e-10 in solid runs). Would turn a start number into link strengths.")
print("  Dark-matter amount: TRPG majoron start angle (measured). Dark-energy amount: the tension (one number).")
print("  Net: about 9 start facts -> 1 condition (smooth start) + 2 ripple numbers (or 1 if the Higgs/stiffness starter holds) + 2 amounts.")
