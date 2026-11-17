# Max and min: how close to zero the lightest grip can get. For each family, the grips (sqrt-weights) lie on a cone
# round the cube's body diagonal: s_i = A(1 + e cos(d + 120 i)), e = sqrt2 tan(opening). The lightest grip is exactly ZERO
# at a wall twist d_wall = acos(-1/e) - 120 deg. d = 0 means the two light ones are equal. Where does each family sit?
import numpy as np
def ana(m):
    s=np.sqrt(np.sort(np.array(m,float))); n=np.ones(3)/np.sqrt(3)
    par=s@n; perp=s-par*n; op=np.degrees(np.arctan(np.linalg.norm(perp)/par)); e=np.sqrt(2)*np.tan(np.radians(op))
    # twist: angle of perp from the direction of the heaviest
    e3=np.array([0,0,1.0]); b=e3-(e3@n)*n; b/=np.linalg.norm(b)
    d=np.degrees(np.arccos(np.clip(perp@b/np.linalg.norm(perp),-1,1)))
    wall=np.degrees(np.arccos(-1/e))-120 if e>1 else None
    return op,e,d,wall,s[0]/s[2]
F={"up (M_Z)":[1.27e-3,0.619,171.7],"up (own scale)":[2.16e-3,1.273,162.5],"down (M_Z)":[2.67e-3,0.0550,2.86],"down (own scale)":[4.70e-3,0.0935,4.183],
   "leptons":[0.51099895,105.6583755,1776.93],"neutrinos (m1 0.36 meV)":[0.36,(0.36**2+74.2)**0.5,(0.36**2+2510)**0.5],
   "neutrinos (m1 = 0)":[0,74.2**0.5,2510**0.5]}
print(f"{'family':26s} opening  e     twist   wall    twist/wall  lightest grip / heaviest grip")
for k,m in F.items():
    op,e,d,w,r=ana(m); print(f"{k:26s} {op:6.2f} {e:6.3f} {d:6.2f} {('%6.2f'%w) if w else '  none'}  {('%.3f'%(d/w)) if w else '  -  '}      {r:.2e}")
print("\nNo wall exists below opening 35.26 deg (face diagonal): there the lightest can never reach zero.")
print("Maximum opening 54.74 deg (edge): wall at twist 0, both light ones exactly zero.")
# lepton: what twist would put electron at zero, and how far is 2/9
print("lepton wall = 15 deg exactly = pi/12 =",np.degrees(np.pi/12),"; 2/9 rad =",np.degrees(2/9),"; ratio",(2/9)/(np.pi/12))
# random check: fraction of random families with twist/wall >= 0.84 given a wall
rng=np.random.default_rng(1); c=n_=0
for _ in range(200000):
    m=10**rng.uniform(-3,3,3); op,e,d,w,r=ana(m)
    if w and w>0: n_+=1; c+= (d/w>=0.84)
print("random log-uniform weights over 6 decades: share with twist/wall >= 0.84:",round(c/n_,3))
