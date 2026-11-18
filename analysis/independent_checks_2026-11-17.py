import numpy as np
from scipy.optimize import brentq
# 1. Koide tau from e, mu (PDG 2024)
me,mmu,mtau,dtau=0.51099895,105.6583755,1776.86,0.12
def Q(a,b,c,s=(1,1,1)): return (a+b+c)/(s[0]*np.sqrt(a)+s[1]*np.sqrt(b)+s[2]*np.sqrt(c))**2
t=brentq(lambda x:Q(me,mmu,x)-2/3,1000,3000); print(f"tau from Koide: {t:.3f} MeV  ({(t-mtau)/dtau:+.2f} sigma)")
# 2. neutrinos: NuFIT 6.0 NO splittings
d21,d31=7.49e-5,2.513e-3
def nus(m1,io=False):
    if not io: return m1,np.sqrt(m1**2+d21),np.sqrt(m1**2+d31)
    m3=m1; m1_=np.sqrt(m3**2+d31-d21); m2=np.sqrt(m3**2+d31); return m3,m1_,m2
for label,s in [("all +",(1,1,1)),("lightest -",(-1,1,1))]:
  for io in (False,True):
    f=lambda m: Q(*nus(m,io),s=s)-2/3
    ms=np.logspace(-7,0,4000); v=[f(m) for m in ms]; roots=[]
    for i in range(len(ms)-1):
        if np.sign(v[i])!=np.sign(v[i+1]): roots.append(brentq(f,ms[i],ms[i+1]))
    for r in roots:
        m=nus(r,io); print(f"  {label:10s} {'IO' if io else 'NO'}: masses {[f'{x*1e3:.2f}' for x in m]} meV sum {sum(m):.4f} eV")
    if not roots: print(f"  {label:10s} {'IO' if io else 'NO'}: no solution")
# 3. 0vbb range for NO with lightest- solution
f=lambda m: Q(*nus(m),s=(-1,1,1))-2/3
m1=brentq(f,1e-6,2e-3); m=nus(m1)
s12,s13=0.307,0.02195
c12,c13=1-s12,1-s13
U2=[c12*c13,s12*c13,s13]
vals=[abs(U2[0]*m[0]+U2[1]*m[1]*np.exp(1j*a)+U2[2]*m[2]*np.exp(1j*b)) for a in np.linspace(0,2*np.pi,181) for b in np.linspace(0,2*np.pi,181)]
print(f"0vbb m_bb range: {min(vals)*1e3:.2f} to {max(vals)*1e3:.2f} meV (m1={m1*1e3:.3f} meV)")
# 4. TM1: sin2 th12 = 1 - 2/(3(1-s13)) ... TM1: |U_e1|^2=2/3 -> c12^2 c13^2 = 2/3
s12tm1=1-2/(3*(1-s13)); print(f"TM1 sin2th12 = {s12tm1:.4f} (NuFIT 0.307 +0.012-0.011 -> {(s12tm1-0.307)/0.012:+.1f} sigma)")
# TM1 delta from th23: cos d = (1-5 s13^2)(... ) use standard formula (Albright-Rodejohann)
for s23 in (0.470,0.561,0.50):
    t12=np.sqrt(s12tm1/(1-s12tm1)); th13=np.arcsin(np.sqrt(s13)); th23=np.arcsin(np.sqrt(s23))
    # cos delta = - (cot 2th23 (1 - 5 s13^2)) / (2 sqrt2 s13 sqrt(1-3 s13^2))
    cd=-(1/np.tan(2*th23))*(1-5*s13)/(2*np.sqrt(2)*np.sqrt(s13)*np.sqrt(1-3*s13))
    d=np.degrees(np.arccos(np.clip(cd,-1,1))); print(f"  s23={s23}: delta = {d:.1f} or {360-d:.1f} deg")
