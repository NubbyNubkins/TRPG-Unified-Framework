"""Borrowed test 2: loop quantum gravity's smallest area as TRPG's grain of space.
Area spectrum: A = 8 pi gamma l_P^2 sqrt(j(j+1)); smallest step (j = 1/2): A_min = 4 sqrt(3) pi gamma l_P^2
(Rovelli & Smolin 1995; Ashtekar & Lewandowski 1997). gamma from black-hole entropy: 0.2375 (Domagala & Lewandowski 2004; Meissner 2004).
TEST FIXED BEFORE RUNNING:
  grain length = side of a square of area A_min, so grain energy = M_Planck / c with c = sqrt(4 sqrt3 pi gamma), gamma = 0.2375.
  Put that c into TRPG's gravity-leak strength calculation (three-loop, leak_establish_three_loop_2026-11-16.py, unchanged)
  and find the strong strength alpha_s(MZ) for which all three forces close exactly.
  PASS: within 2 sigma of PDG 0.1180 +- 0.0009.  FAIL: beyond 3 sigma.
  Other conventions (gamma = 0.274 from the SU(2) count, Ghosh & Mitra 2005, Engle, Noui & Perez 2010; disc radius instead of square side)
  are shown for completeness only. Picking whichever fits after the fact is not evidence.
"""
import math, io, contextlib
src=open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'leak_establish_three_loop_2026-11-16.py')).read()
p1=src.split('print("A. grain')[0]
p2=src[src.index('nG=3; nt=1'):src.index('# sanity')]
ns={}; exec(p1,ns); exec(p2,ns)
exec('''
def solve(c=1.0, MP=MPL, n=3, **kw):
    def meet(ML):
        s_=predict_s2(thr=thrL(ML,n),**kw); return mismatch(s2=s_,thr=thrL(ML,n),**kw)[1]
    f=lambda x: math.log10(c*meet(10**x)**2/MP)-x
    M=10**brentq(f,3,13.5,xtol=1e-6); return M,meet(M),predict_s2(thr=thrL(M,n),**kw)
def a3_close(c): return brentq(lambda a3: solve(c=c,a3=a3)[2]-S2,0.110,0.126,xtol=1e-6)
''',ns)
A=lambda g: 4*math.sqrt(3)*math.pi*g
cases=[("MAIN: square side, gamma 0.2375",math.sqrt(A(0.2375))),
       ("check: c = 1 (grain = Planck)",1.0),("check: c = sqrt2 (TRPG main grain)",math.sqrt(2)),
       ("other: square side, gamma 0.274",math.sqrt(A(0.274))),
       ("other: disc radius, gamma 0.2375",math.sqrt(A(0.2375)/math.pi)),
       ("other: disc radius, gamma 0.274",math.sqrt(A(0.274)/math.pi))]
for lab,c in cases:
    a=ns['a3_close'](c); M,F,s=ns['solve'](c=c)
    z=(a-0.1180)/0.0009; v="PASS" if abs(z)<2 else ("FAIL" if abs(z)>3 else "PARTIAL")
    print(f"{lab:38s} c = {c:.3f}  grain = {ns['MPL']/c:.3e} GeV  MMG {M:.2e}  closes at alpha_s = {a:.4f}  ({z:+.1f} sigma vs PDG)  -> {v}")
