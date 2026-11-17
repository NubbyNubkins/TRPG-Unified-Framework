"""What stops a relationship going all the way to nothing? Test whether the
electron-gravity balance point sqrt(m_e x M_Planck) sets the far-side mass that closes the strength gap, now with two-loop running.
Main line: one strength at the meeting (strong = weak = hypercharge, sin^2 = 3/8 there); three spin-1/2 far-side weak pairs below the meeting,
coloured far-side parts at the meeting. Credit: Machacek & Vaughn 1983-85; Jones 1982; Arvanitaki, Davis, Graham & Wacker 2004 and
Giudice & Romanino 2004 (two-loop terms for spin-1/2 weak pairs, as in split supersymmetry); PDG 2024."""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def pr(*a): print(*a, flush=True)
mZ, ainv, a3, s2m = 91.1876, 127.951, 0.1180, 0.23122
b1 = np.array([41/10, -19/6, -7.0]); B2 = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12], [11/10, 9/2, -26]]); ct = np.array([17/10, 3/2, 2])
db1 = np.array([2/5, 2/3, 0]); dB2 = np.array([[9/50, 9/10, 0], [3/10, 49/6, 0], [0, 0, 0]])     # one spin-1/2 weak pair (both charges)
def run(s2, MF, npair=3, two=True):
    a = 1/ainv; g1 = math.sqrt(4*math.pi*a/(1 - s2)*5/3); g2 = math.sqrt(4*math.pi*a/s2); g3 = math.sqrt(4*math.pi*a3)
    yt = 0.9369   # MS-bar top Yukawa at the Z (approx.)
    def rhs(t, y):
        g = y[:3]; Y = y[3]; L = t > math.log(MF/mZ) if MF else False
        bb = b1 + (npair*db1 if L else 0); BB = B2 + (npair*dB2 if L else 0)
        g2v = g**2
        dg = g**3*bb/(16*math.pi**2)
        if two: dg += g**3/(16*math.pi**2)**2*(BB @ g2v - ct*Y**2)
        dY = Y/(16*math.pi**2)*(4.5*Y**2 - 8*g2v[2] - 2.25*g2v[1] - 17/20*g2v[0])
        return [*dg, dY]
    ev = lambda t, y: y[1] - y[2]; ev.terminal = True
    sol = solve_ivp(rhs, [0, 60], [g1, g2, g3, yt], events=ev, rtol=1e-9, atol=1e-12, max_step=0.5)
    if not sol.t_events[0].size: return None
    tU = sol.t_events[0][0]; gU = sol.y_events[0][0]
    return gU[0] - gU[1], mZ*math.exp(tU), 4*math.pi/gU[1]**2
def predict(MF, npair=3, two=True):
    f = lambda s: run(s, MF, npair, two)[0]
    return brentq(f, 0.19, 0.26, xtol=1e-7)
pr("Check: Standard Model alone, one loop and two loops")
for two in (False, True):
    s = predict(None, 0, two); pr(f"  {'two' if two else 'one'} loop(s): sin^2 = {s:.4f} ({100*(s/s2m - 1):+.1f}%)")
pr("\nThree spin-1/2 far-side weak pairs at mass M (coloured far-side parts at the meeting):")
for M in (1e3, 1e5, 1e6, 1e7, 3e7, 1e8, 1e9):
    s = predict(M); r = run(s, M); pr(f"  M = {M:8.1e} GeV: sin^2 = {s:.5f} ({100*(s/s2m - 1):+.2f}%), meeting {r[1]:.2e} GeV, 1/strength there {r[2]:.1f}")
Mneed1 = 10**brentq(lambda x: predict(10**x, 3, False) - s2m, 4, 11)
Mneed2 = 10**brentq(lambda x: predict(10**x, 3, True) - s2m, 4, 11)
pr(f"\nMass that closes the gap exactly: one loop {Mneed1:.3e} GeV; two loops {Mneed2:.3e} GeV")
me = 0.51099895e-3
for name, M in (("sqrt(electron x Planck mass)", math.sqrt(me*1.220890e19)), ("sqrt(electron x reduced Planck mass)", math.sqrt(me*2.435e18)),
                ("top^2/electron", 172.57**2/me), ("Higgs mass^2/electron", 125.25**2/me), ("Z^2/electron", mZ**2/me)):
    s = predict(M); pr(f"  {name:38s} = {M:.3e} GeV -> sin^2 = {s:.5f} ({100*(s/s2m - 1):+.2f}%, {(s - s2m)/0.00004:+.0f} sigma of the measurement)")
pr("  Remaining unknowns that move these numbers: threshold details at the meeting (coloured far-side parts, carriers), roughly +-0.3% in sin^2.")
