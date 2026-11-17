# How the weights of the particles enter the Hubble value; speeds of things at one temperature
# Same simple three-level atom as hubble_late_resistance (differences from its own baseline only). Sky angle of the glow held fixed.
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
c=2.99792458e10; kB=8.617333e-5; hbar=6.582120e-16; mec2=510998.95; sigT0=6.6524587e-25; mH=1.6735e-24; Mpc=3.0857e24
T0=2.7255; Yp=0.245; wg=2.4729e-5; Tnu=1.9454*kB   # eV
P=dict(wb=0.02237,wc=0.1200,mnu=[0.36e-3,8.62e-3,50.1e-3],d=0.0)
def wnu(z,p):  # three neutrinos, each: massless share x sqrt(1+(m/(3.15 T))^2)
    return sum(0.2271*1.0132*wg*(1+z)**4*np.sqrt(1+(m/(3.15*Tnu*(1+z)))**2) for m in p["mnu"])
def Hz(z,h,p):
    f=lambda zz: wg*(1+zz)**4+wnu(zz,p)+(p["wb"]+p["wc"])*(1+zz)**3
    return 100*1e5/Mpc*np.sqrt(f(z)+(h*h-f(0.0)))
def zstar(h,p):
    d=p["d"]; s=1+d; E1,E2=13.6057*s,3.4014*s; E21=E1-E2; lam=1.215670e-5/s; Lam=8.2246*s; sig=sigT0/s**2; nH0=(1-Yp)*p["wb"]*1.8788e-29/mH/p.get("mp",1.0)
    alpha=lambda T: 1.14*4.309e-13*(T/s/1e4)**-0.6166/(1+0.6703*(T/s/1e4)**0.53)/s**2
    beta=lambda T: alpha(T)*(mec2*s*kB*T/(2*np.pi))**1.5/(hbar*c)**3*np.exp(-E2/(kB*T))
    def saha(z):
        T=T0*(1+z); n=nH0*(1+z)**3; A=(mec2*s*kB*T/(2*np.pi))**1.5/(hbar*c)**3*np.exp(-E1/(kB*T))/n; return (-A+np.sqrt(A*A+4*A))/2
    def rhs(z,x):
        x=min(max(x[0],1e-10),1.0); T=T0*(1+z); n=nH0*(1+z)**3; H=Hz(z,h,p); K=lam**3/(8*np.pi*H); n1=(1-x)*n; b=beta(T); C=(1+K*Lam*n1)/(1+K*(Lam+b)*n1)
        return [C*(alpha(T)*x*x*n-b*(1-x)*np.exp(-E21/(kB*T)))/(H*(1+z))]
    zs=np.linspace(2200,200,4001); z0=1650*s; sol=solve_ivp(rhs,[z0,zs[-1]],[saha(z0)],t_eval=zs[zs<=z0],method="Radau",rtol=1e-6,atol=1e-10)
    x=np.ones_like(zs); x[zs<=z0]=sol.y[0]; x[zs>z0]=[min(1,saha(z)) for z in zs[zs>z0]]
    dt=x*nH0*(1+zs)**2*sig*c/np.array([Hz(z,h,p) for z in zs]); tau=np.array([np.trapezoid(dt[i:][::-1],zs[i:][::-1]) for i in range(0,len(zs))]); return zs[np.argmax(dt*np.exp(-tau))]
def rs(z,h,p): R=lambda zz: 3*p["wb"]/(4*wg)/(1+zz); return quad(lambda zz: c/np.sqrt(3*(1+R(zz)))/Hz(zz,h,p),z,1e7,limit=400)[0]/Mpc
def DA(z,h,p): return quad(lambda zz: c/Hz(zz,h,p),0,z,limit=400)[0]/Mpc
h0=0.674; z0=zstar(h0,P); th0=rs(z0,h0,P)/DA(z0,h0,P)
def H0for(p):
    zp=zstar(h0,p); return 100*brentq(lambda h: rs(zp,h,p)/DA(zp,h,p)-th0,0.5,0.95)
print("(1) Change one weight by 1% and keep the sky the same: how far does the Hubble value move? (baseline 67.4)")
tests=[("electron heavier by 1% (atoms form earlier)",dict(P,d=0.01)),
       ("proton and neutron heavier by 1% (ordinary matter weighs 1% more, same count)",dict(P,wb=P["wb"]*1.01,mp=1.01)),
       ("dark matter 1% more (weight x count; the theory cannot yet say which)",dict(P,wc=P["wc"]*1.01)),
       ("all three neutrinos heavier by 1%",dict(P,mnu=[m*1.01 for m in P["mnu"]]))]
for n,p in tests: print(f"   {n:78s} {H0for(p)-67.4:+.2f}")
print("\n(2) The neutrinos' weights (the only weights still unmeasured)")
for lab,m in (("weightless",[0,0,0]),("the theory: 0.36, 8.6, 50.1 meV (sum 0.059 eV)",P["mnu"]),("sum 0.12 eV (survey upper limit)",[0.032,0.033,0.059]),("sum 0.30 eV",[0.095,0.096,0.109])):
    print(f"   {lab:52s} H0 = {H0for(dict(P,mnu=m)):.1f}")
print("\n(3) Speeds at one temperature: v = share of light speed, for movement energy 3/2 kT (slow things) or momentum 3.15 kT (neutrinos)")
def v_th(m,kT): return np.sqrt(3*kT/m)*299792.458
for lab,kT in (("when atoms formed (0.26 eV)",0.26),("today's glow (0.000235 eV)",2.349e-4)):
    print(f"   {lab}: electron {v_th(0.511e6,kT):.3g} km/s; proton {v_th(938.27e6,kT):.3g} km/s; light 299,792 km/s")
pn=3.15*1.676e-4
for i,m in enumerate(P["mnu"]): print(f"   relic neutrino {i+1} today: {299792.458*pn/np.sqrt(pn*pn+m*m):,.0f} km/s ({100*pn/np.sqrt(pn*pn+m*m):.1f}% of light speed)")
print("\n(4) When each kind stopped behaving like light (weight = 3 x temperature) and began to slow the stretch as held matter")
for lab,m in (("proton",938.27e6),("electron",0.511e6),("neutrino 3",50.1e-3),("neutrino 2",8.62e-3),("neutrino 1",0.36e-3)):
    print(f"   {lab:11s} light at {m/3:.3g} eV")
