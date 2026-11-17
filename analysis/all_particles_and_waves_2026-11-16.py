"""All particles laid out; the lepton V (charged vs neutrino partner) against the quark V; each particle's own wave."""
import math, numpy as np
from scipy.optimize import brentq
hbarc=197.3269804e-15  # MeV*m
me,mmu,mtau=0.51099895,105.6583755,1776.93
D21,D31=7.42e-5,2.51e-3; dD21,dD31=0.21e-5,0.027e-3   # eV^2 (NuFIT, normal order)
def nus(m1): return np.array([m1,math.sqrt(m1**2+D21),math.sqrt(m1**2+D31)])
m1T=0.00036
def rate(lo,hi):
    w=np.log(np.array(hi)/np.array(lo)); return (w[2]-w[1])/(w[1]-w[0]), w
ml=np.array([me,mmu,mtau])*1e6   # eV
r,w=rate(nus(m1T),ml); print("lepton V (charged / neutrino partner) with TRPG m1=0.00036 eV: ratios",np.exp(w),"rate",round(r,4))
# quark V at top footing and Z
mu=np.array([1.242,627.8,162500.]); md=np.array([2.678,53.36,2873.])
print("quark V rate (top footing)",round(rate(md,mu)[0],4))
# m1 that gives lepton V rate = 1/2, and = quark V rate at Z 0.504
for t in (0.5,0.4856,0.5035,0.5294):
    f=lambda x: rate(nus(10**x),ml)[0]-t
    xs=np.linspace(-7,-1.5,400); v=[f(x) for x in xs]; roots=[brentq(f,xs[i],xs[i+1]) for i in range(len(xs)-1) if v[i]*v[i+1]<0]
    print(f"  lepton V rate = {t}: m1 = {[f'{10**x:.2e}' for x in roots]} eV, sum = {[round(nus(10**x).sum(),4) for x in roots]}")
# error on lepton V rate at m1T from splittings
rng=np.random.default_rng(0); R=[]
for _ in range(20000):
    a=D21+dD21*rng.normal(); b=D31+dD31*rng.normal(); n=np.array([m1T,math.sqrt(m1T**2+a),math.sqrt(m1T**2+b)]); R.append(rate(n,ml)[0])
print("lepton V rate at m1=0.00036: %.4f +- %.4f (splitting errors only)"%(np.mean(R),np.std(R)))
# how rate depends on m1
for m1 in (1e-5,1e-4,2e-4,3.6e-4,5e-4,1e-3,3e-3):
    print(f"   m1={m1:.1e}: rate {rate(nus(m1),ml)[0]:.4f}")
# waves: reduced Compton wavelength hbar/(mc)
P={"neutrino 1 (TRPG)":m1T*1e-6,"neutrino 2":nus(m1T)[1]*1e-6,"neutrino 3":nus(m1T)[2]*1e-6,"electron":me,"up":2.16,"down":4.70,"strange":93.5,"muon":mmu,
   "charm":1273,"tau":mtau,"bottom":4183,"W":80369,"Z":91188,"Higgs":125200,"top":172570,"proton":938.272,"pion":139.57}
print("\nreduced Compton wavelength (the particle's own wave length), metres:")
for k,m in sorted(P.items(),key=lambda kv:kv[1]): print(f"   {k:18s} {m:12.5g} MeV  ->  {hbarc/m:.3e} m")
lP=1.616255e-35; print("grain (Planck length) %.3e m; grain 2 (M_Pl/sqrt2) wave %.3e m"%(lP, lP*math.sqrt(2)))
print("proton bubble radius 0.841e-15 m = %.2f x its own wave"%(0.841e-15/(hbarc/938.272)))
print("proton radius / own wave: %.4f +- %.4f"%(0.84075e-15/(hbarc/938.272088), 0.00064e-15/(hbarc/938.272088)))
