"""Are the up and down openings set TOGETHER? (an earlier test question: the two ladders share one state.)
Opening from Koide-type ratio Q = sum m / (sum sqrt m)^2 ; tan(opening) = sqrt2 * link/self.
Error bars by Monte Carlo on PDG 2024 masses; second footing = running masses at M_Z (Xing, Zhang & Zhou 2020 style values).
Every simple joint relation is tried against a list of simple targets, and a luck check counts how often random
angle pairs give a hit as good."""
import numpy as np, itertools, math
rng=np.random.default_rng(11)
def ang(m):
    m=np.asarray(m,float); Q=m.sum(-1)/np.sqrt(m).sum(-1)**2; e=np.sqrt(np.clip(2*(3*Q-1),0,None)); return np.degrees(np.arctan(e/np.sqrt(2)))
def asym(c,lo,hi,n): # asymmetric gaussian
    z=rng.standard_normal(n); return np.where(z<0,c+z*lo,c+z*hi)
N=200000
pdg_d=np.stack([asym(4.70,0.17,0.48,N),asym(93.5,3.4,8.6,N),asym(4183,7,7,N)],-1)
pdg_u=np.stack([asym(2.16,0.26,0.49,N),asym(1273,4.6,4.6,N),asym(172570,290,290,N)],-1)
mz_d=np.stack([asym(2.67,0.19,0.19,N),asym(53.2,4.6,4.6,N),asym(2839,26,26,N)],-1)
mz_u=np.stack([asym(1.23,0.21,0.21,N),asym(620,17,17,N),asym(168260,750,750,N)],-1)
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
foot={}
for nm,(D,U) in {'PDG own-scale masses':(pdg_d,pdg_u),'all at M_Z':(mz_d,mz_u)}.items():
    td,tu=ang(D),ang(U); foot[nm]=(td,tu)
    P(f"{nm}: down {np.median(td):.2f} (+{np.percentile(td,84)-np.median(td):.2f}/-{np.median(td)-np.percentile(td,16):.2f})  up {np.median(tu):.2f} (+{np.percentile(tu,84)-np.median(tu):.2f}/-{np.median(tu)-np.percentile(tu,16):.2f})")
# joint quantities
def joint(td,tu):
    a,b=np.radians(td),np.radians(tu)
    return {'tan d * tan u':np.tan(a)*np.tan(b),'tan u / tan d':np.tan(b)/np.tan(a),'d + u (deg)':td+tu,'u - d (deg)':tu-td,
            'cos d * cos u':np.cos(a)*np.cos(b),'sin d * sin u':np.sin(a)*np.sin(b),'tan^2 d + tan^2 u':np.tan(a)**2+np.tan(b)**2,
            '(u-45)/(d-45)':(tu-45)/(td-45),'r_d^2 + r_u^2':(np.tan(a)**2+np.tan(b)**2)/2}
targets={'1':1,'sqrt2':math.sqrt(2),'sqrt3':math.sqrt(3),'2':2,'3':3,'1/2':0.5,'1/3':1/3,'2/3':2/3,'3/2':1.5,'4/3':4/3,'5/4':1.25,'9/4':2.25,
         'pi/2':math.pi/2,'e':math.e,'phi':(1+5**0.5)/2,'90':90,'100':100,'7':7,'5':5,'4':4,'8/3':8/3,'5/3':5/3,'3/4':0.75}
P("\nJoint relations within 1 error bar of a simple number (both footings must agree to count):")
hits=[]
for q in joint(np.array([48.]),np.array([52.])).keys():
    ok=[]
    for nm,(td,tu) in foot.items():
        v=joint(td,tu)[q]; med=np.median(v); sd=(np.percentile(v,84)-np.percentile(v,16))/2
        best=min(targets.items(),key=lambda kv: abs(kv[1]-med)); ok.append((nm,med,sd,best,abs(best[1]-med)/sd))
    line=" | ".join(f"{nm}: {med:.4f}±{sd:.4f} (nearest {b[0]}, {pull:.1f}σ)" for nm,med,sd,b,pull in ok)
    P(f"  {q:>18}: {line}")
    if all(o[4]<1 for o in ok) and ok[0][3][0]==ok[1][3][0]: hits.append((q,ok[0][3][0]))
P(f"\nrelations holding at BOTH footings within 1σ, same target: {hits if hits else 'none'}")
# luck: random angle pairs in same range, same widths -> chance of >=len(hits) hits
sdd=np.std(foot['PDG own-scale masses'][0]); sdu=np.std(foot['PDG own-scale masses'][1])
cnt=0; M=4000
for _ in range(M):
    d0=rng.uniform(45,55); u0=rng.uniform(45,60); dd=d0+rng.standard_normal(4000)*sdd; uu=u0+rng.standard_normal(4000)*sdu
    J=joint(dd,uu); k=0
    for q,v in J.items():
        med=np.median(v); sd=(np.percentile(v,84)-np.percentile(v,16))/2
        if min(abs(t-med) for t in targets.values())<sd: k+=1
    if k>=max(1,len(hits)): cnt+=1
P(f"luck check: random angle pairs (down 45-55, up 45-60, same error bars) give at least {max(1,len(hits))} such hit(s) {100*cnt/M:.0f}% of the time (one footing only, so this is generous to the real case)")
open('joint_quark_openings_2026-11-16.txt','w').write('\n'.join(out)+'\n')
