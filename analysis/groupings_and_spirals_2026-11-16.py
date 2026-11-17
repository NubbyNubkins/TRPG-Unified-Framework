# Do values land in groups more than chance? any log-periodic (spiral/fractal) spacing?
import numpy as np
rng=np.random.default_rng(1)
H={ # log10(GeV), independent measured scales (no derived duplicates)
 "nu1":np.log10(3.6e-13),"nu2":np.log10(8.6e-12),"nu3":np.log10(5.0e-11),"dark energy":np.log10(2.24e-12),"CMB today":np.log10(6.3e-13),
 "e":np.log10(5.11e-4),"u":np.log10(2.16e-3),"d":np.log10(4.70e-3),"s":np.log10(0.0935),"mu":np.log10(0.10566),
 "p":np.log10(0.93827),"c":np.log10(1.27),"tau":np.log10(1.77686),"b":np.log10(4.18),
 "W":np.log10(80.37),"Z":np.log10(91.19),"H":np.log10(125.2),"t":np.log10(172.6),"v":np.log10(246.2),
 "Planck":np.log10(1.221e19)}
names=list(H); x=np.array([H[k] for k in names])
def close_pairs(x,d): 
    x=np.sort(x); return sum(1 for i in range(len(x)) for j in range(i+1,len(x)) if x[j]-x[i]<d)
lo,hi=x.min(),x.max()
for d in (0.05,0.1,0.2):
    obs=close_pairs(x,d); sims=np.array([close_pairs(rng.uniform(lo,hi,len(x)),d) for _ in range(5000)])
    print(f"pairs closer than {d} decades: observed {obs}, random mean {sims.mean():.1f}, chance of >= observed {np.mean(sims>=obs):.4f}")
# fairer: merge groups whose closeness is already explained (Higgs family: W,Z,H,t,v set by v; light-neutrino/CMB/DE)
merged=dict(H); 
for k in ("W","Z","H","t"): merged.pop(k)
x2=np.array(list(merged.values()))
# compare against a fairer null: random points but same 'clumpy' density (bootstrap of gaps)
gaps=np.diff(np.sort(x2))
for d in (0.05,0.1,0.2):
    obs=close_pairs(x2,d); sims=[]
    for _ in range(5000):
        g=rng.permutation(gaps); sims.append(close_pairs(np.concatenate([[0],np.cumsum(g)]),d))
    sims=np.array(sims); u=np.array([close_pairs(rng.uniform(x2.min(),x2.max(),len(x2)),d) for _ in range(3000)])
    print(f"[Higgs group counted once] < {d}: observed {obs}; uniform-random chance {np.mean(u>=obs):.3f}")
sx=np.sort(x2); print("closest independent pairs:",sorted([(round(sx[i+1]-sx[i],3),[k for k,v in merged.items() if abs(v-sx[i])<1e-9][0],[k for k,v in merged.items() if abs(v-sx[i+1])<1e-9][0]) for i in range(len(sx)-1)])[:6])
# outliers: biggest gaps
print("biggest empty stretches (decades):",sorted([(round(sx[i+1]-sx[i],2),[k for k,v in merged.items() if abs(v-sx[i])<1e-9][0],[k for k,v in merged.items() if abs(v-sx[i+1])<1e-9][0]) for i in range(len(sx)-1)])[-3:])
# log-periodicity (spiral that climbs by a fixed factor each turn): Rayleigh power over periods 0.3-6 decades
xs=x[(x>-4)&(x<3)]   # charged fermions + bosons (the dense region)
P=np.linspace(0.3,6,3000)
def maxpow(v): return max(abs(np.exp(2j*np.pi*v/p).sum())**2/len(v) for p in P)
obs=maxpow(xs); best=P[np.argmax([abs(np.exp(2j*np.pi*xs/p).sum())**2/len(xs) for p in P])]
sims=np.array([maxpow(rng.uniform(xs.min(),xs.max(),len(xs))) for _ in range(1000)])
print(f"\nspiral test (charged fermions + W,Z,H,v; {len(xs)} values): best repeat every {best:.3f} decades (x{10**best:.2f}), power {obs:.2f}; chance with random values {np.mean(sims>=obs):.3f}")
