import math, numpy as np
def cone(m):
    r=np.sqrt(np.array(m,float)); s=r.mean(); c=(r/s-1)        # c_k = eps*cos(delta + 2 pi k/3), k=1,2,3 (e, mu, tau)
    eps=math.sqrt(2*np.mean(c**2)); best=None
    for d in np.linspace(-math.pi,math.pi,200001):
        v=eps*np.cos(d+2*math.pi*np.arange(1,4)/3); e=np.sum((v-c)**2)
        if best is None or e<best[0]: best=(e,d)
    return eps, best[1], math.degrees(math.atan(eps/math.sqrt(2)))
sets={"each weighed at its own settling (pole masses)":[0.51099895,105.6583755,1776.93],
      "all weighed at one temperature: Z (91 GeV), running":[0.486570161,102.7181359,1746.24],
      "all at the top's energy (173 GeV), running":[0.485289396,102.4467419,1741.63],
      "all at 1 TeV, running (approx.)":[0.48,101.4,1725.0]}
for k,m in sets.items():
    eps,d,op=cone(m); print(f"{k:55s} opening {op:.4f} deg, twist {d:+.6f} rad ({math.degrees(d):+.4f} deg); 2/9 = 0.222222")
