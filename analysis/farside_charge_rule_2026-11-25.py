"""Independent recheck of the Z6 charge claims (anomalies, operators, generation count),
the domain-wall temperatures, and the gravity-broken majoron, written from scratch.

Every fermion is written as a left-handed Weyl field. A right-handed field psi_R of charge q
enters as its conjugate psi^c with charge -q. Instanton zero-mode counting gives, for a Z_N
symmetry, the mixed anomaly with SU(M): A = sum_i q_i * 2T(r_i) (2T = 1 for a fundamental),
and for gravity A_grav = sum_i q_i * dim(r_i). Two conventions are reported:
  strict  : A_SU3 = A_SU2 = 0 mod N, A_grav = 0 mod N/2   (Ibanez-Ross with the usual relaxation for gravity)
  paper   : 2*A_SU3 = 2*A_SU2 = 0 mod N, A_grav = 0 mod N  (the half-strength form sometimes used)
"""
import math

N = 6
qQ, qu, qd, qH = 0, 0, 0, 0
ql, qe, qN = 3, 3, 3
qPhi = 1
dD, dL = -1, -2          # q(D_R) - q(D_L), q(L_R) - q(L_L)

def weyl_list(ngen, far=True, dL_=dL, dD_=dD):
    """(name, charge, colour_dim, su2_dim) for left-handed Weyl fields; colour and SU(2) multiplicity explicit."""
    F = []
    for g in range(ngen):
        F += [("Q", qQ, 3, 2), ("u^c", -qu, 3, 1), ("d^c", -qd, 3, 1),
              ("l", ql, 1, 2), ("e^c", -qe, 1, 1), ("N^c", -qN, 1, 1)]
        if far:
            qDL = 0; qDR = qDL + dD_
            qLL = 0; qLR = qLL + dL_
            F += [("D_L", qDL, 3, 1), ("D^c", -qDR, 3, 1), ("L_L", qLL, 1, 2), ("L^c", -qLR, 1, 2)]
    return F

def anomalies(F):
    a3 = sum(q * c2 for (_, q, c, s) in F for c2 in [1] if c == 3) * 1  # each colour triplet: 2T = 1, counted once per SU(2) component
    a3 = sum(q * s for (_, q, c, s) in F if c == 3)       # SU(3)^2: triplets, times SU(2) multiplicity
    a2 = sum(q * c for (_, q, c, s) in F if s == 2)       # SU(2)^2: doublets, times colour multiplicity
    ag = sum(q * c * s for (_, q, c, s) in F)             # grav^2: every component
    return a3, a2, ag

def verdict(a3, a2, ag):
    strict = (a3 % N == 0, a2 % N == 0, ag % (N // 2) == 0)
    paper = ((2 * a3) % N == 0, (2 * a2) % N == 0, ag % N == 0)
    return strict, paper

print("Z6 ANOMALY RECHECK (three generations)")
for far in (False, True):
    a = anomalies(weyl_list(3, far))
    s, p = verdict(*a)
    print(f"  far-side states {'included' if far else 'absent  '}: A_SU3={a[0]:3d}  A_SU2={a[1]:3d}  A_grav={a[2]:3d}"
          f" | strict pass={s} | paper's convention pass={p}")
print("  A count that treats the lepton doublet and e^c as cancelling gets 9 and 30. Field by field, the lepton doublet,")
print("  e^c and N^c give 2*3 - 3 - 3 = 0 per generation, so the SM part is 0, not 9, and the far part is 21.")

print("\nGENERATION SCAN (does N_g = 3 come out uniquely?)")
for ng in range(1, 9):
    a = anomalies(weyl_list(ng, True)); s, p = verdict(*a)
    b3 = 11 - (4 / 3) * ng - (2 / 3) * ng   # SM quarks plus one vector-like D per generation (Dirac: 2/3 each)
    print(f"  N_g={ng}: strict all-pass={all(s)}  paper-convention all-pass={all(p)}  b3={b3:+.2f} (asymptotically free: {b3 > 0})")

print("\nCAN ANY far-side charge split give both a dimension-5 L_F mass and an SU(2)-safe count (strict)?")
for dl in range(6):
    dim4 = any((k * qPhi + dl) % N == 0 for k in (1, -1)) or dl % N == 0
    dim5 = any((k * qPhi + dl) % N == 0 for k in (2, -2))
    a = anomalies(weyl_list(3, True, dL_=dl - 6 if dl else 0))
    print(f"  dq_L={dl}: lowest mass operator dim-4 or bare={dim4}, dim-5={dim5}, A_SU2 mod 6 = {a[1] % 6}")

print("\nOPERATOR CHARGES")
print(f"  QQQl: {(3*qQ+ql)%N}  u^c u^c d^c e^c: {(-3*qu - qe)%N}  N_R N_R: {(2*qN)%N}  l H N: {(ql+qH-qN)%N}"
      "  (charge 3 = lepton parity; any Delta B = 1, Delta L = 1 operator is odd. The 81 flavour copies are one condition.)")

# ---------------------------------------------------------------- domain walls
MPl, MPlr, E, gs = 1.2209e19, 2.435e18, 8.63e18, 106.75
v = 4.6e13
sig = v ** 3
T_of_H = lambda H: math.sqrt(H * MPl / (1.66 * math.sqrt(gs)))
Tdom = T_of_H(sig / MPlr ** 2)
print("\nZ6 DOMAIN WALLS")
print(f"  T_dom = {Tdom:.3e} GeV")
for n in (5, 6, 7):
    breaks = (n * qPhi) % N != 0
    Tann = T_of_H((v ** n / E ** (n - 4)) / sig)
    print(f"  Phi^{n}/E^{n-4}: Z6 charge {n % 6} -> {'breaks Z6' if breaks else 'Z6-INVARIANT, cannot bias the vacua'};"
          f" T_ann = {Tann:.3e} GeV")
print("  Note: an anomaly-free (gauged) discrete symmetry cannot be broken by gravity (Krauss and Wilczek 1989);")
print("  the walls argument needs the opposite. The two claims cannot both hold.")

# ---------------------------------------------------------------- majoron
print("\nMAJORON")
f = 2.1e16
H0, T0, h = 1.437642e-42, 2.348654e-13, 0.674
rho_c = 3 * MPlr ** 2 * H0 ** 2
for m in (3e-19, 1e-18, 1e-17, 1e-16):
    mG = m * 1e-9
    # (a) oscillation starts at 3H = m with g* = 3.36, entropy dilution
    g = 3.36; gs_ = 3.91
    Tosc = math.sqrt(mG * MPlr / (3 * math.sqrt(math.pi ** 2 * g / 90)))
    rho = 0.5 * mG ** 2 * f ** 2 * (gs_ * T0 ** 3) / (gs_ * Tosc ** 3)
    th_a = math.sqrt(0.12 / (rho / rho_c * h ** 2))
    # (b) Hui, Ostriker, Tremaine and Witten (2017), eq. 7: Omega ~ 0.1 (f/1e17)^2 (m/1e-22)^(1/2) theta^2
    th_b = math.sqrt(0.12 / h ** 2 / (0.1 * (f / 1e17) ** 2 * (m / 1e-22) ** 0.5))
    print(f"  m={m:.0e} eV: T_osc={Tosc*1e6:6.1f} keV  theta_i (3H=m) = {th_a:.3f}   theta_i (HOTW 2017 fit) = {th_b:.3f}")
# mass from gravity breaking of U(1)_{B-L} at dimension n: m^2 ~ f^(n-2)/E^(n-4)
for m in (3e-19, 1e-16):
    target = m * 1e-9 / f
    n = 4 + 2 * math.log(target) / math.log(f / E)
    print(f"  to get m_J = {m:.0e} eV by gravity, breaking must start near operator dimension {n:.0f} (dimension 5 gives"
          f" m_J ~ {math.sqrt(f**3/E):.1e} GeV)")
print("  Without inflation the angle theta_i differs from one causal patch to the next; the patch-averaged")
print(f"  <theta^2> = pi^2/3 = {math.pi**2/3:.2f}, so a single tunable theta_i = 0.1-0.6 has no justification.")
for m in (3e-19, 1e-18, 1e-17, 1e-16):
    print(f"  superradiance alpha=0.05..0.45 -> M_BH = {0.05*1.3357e-10/m:.2e} .. {0.45*1.3357e-10/m:.2e} Msun for m = {m:.0e} eV")
