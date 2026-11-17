"""Reading tested: 45 deg is the equilibrium (physical perspective, leptons);
measure every family by its signed distance from it; top half one type, bottom half another; the extremes are the sources.
Opening = angle between (sqrt m) and the cube's body diagonal. Extremes: edge 54.7356, face diagonal 35.2644, body diagonal 0."""
import numpy as np, math
from scipy.optimize import brentq
rng=np.random.default_rng(2)
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
def op(m):
    r=np.sqrt(np.asarray(m,float)); return math.degrees(math.acos(min(1,r.sum()/math.sqrt(3)/np.linalg.norm(r))))
EDGE,FACE=math.degrees(math.acos(1/math.sqrt(3))),math.degrees(math.acos(math.sqrt(2/3)))
def nu(m1): return [m1,math.sqrt(m1**2+7.42e-5),math.sqrt(m1**2+2.51e-3)]  # eV
fam=[('up quarks (at M_Z)',[1.23,620,168260],'coloured, upper of weak pair'),('up quarks (own scale)',[2.16,1273,172570],''),
     ('down quarks (at M_Z)',[2.67,53.2,2839],'coloured, lower of weak pair'),('down quarks (own scale)',[4.70,93.5,4183],''),
     ('charged leptons',[0.51099895,105.658,1776.93],'colourless, lower of weak pair'),
     ('neutrinos (TRPG m1 = 0.36 meV)',nu(3.6e-4),'colourless, upper of weak pair'),('neutrinos (m1 = 0)',nu(0),''),
     ('hidden partners N1-N3',[7.3e5,1.3e9,6.3e10],'far side'),('MMG x3 (one shared weight)',[1,1,1],'far side'),('D partners x3 (at the meeting)',[1,1,1],'far side')]
P("signed distance from the 45 deg equilibrium (+ = toward the edge 'one heavy' at +9.74; - = toward the face diagonal 'two heavy' at -9.74; body diagonal 'all equal' at -45)")
for nm,m,tag in fam:
    t=op(m); P(f"  {nm:<32} opening {t:6.2f}  distance {t-45:+6.2f} deg  ({(t-45)/9.7356*100:+5.0f}% of the way to the extreme)  {tag}")
# mirror test: up vs neutrinos
P("\nmirror test across 45: is the upper partner of each weak pair equally far out on opposite sides?")
for lab,u in (('at M_Z',52.21),('own scale',51.20)):
    tgt=90-u
    try: m1=brentq(lambda x: op(nu(x))-tgt,0,5e-3); s=f"neutrino mirror needs m1 = {m1*1e3:.2f} meV"
    except ValueError: s="no m1 gives the mirror"
    P(f"  up {lab}: +{u-45:.2f}; mirror opening {tgt:.2f}; {s}  (TRPG's own m1 = 0.36 meV gives {op(nu(3.6e-4)):.2f})")
# error on neutrino opening from splittings (NuFIT-like 1 sigma: 7.42+-0.21e-5, 2.51+-0.027e-3)
ops=[op([3.6e-4,math.sqrt(3.6e-4**2+rng.normal(7.42e-5,0.21e-5)),math.sqrt(3.6e-4**2+rng.normal(2.51e-3,0.027e-3))]) for _ in range(20000)]
P(f"  neutrino opening error from the measured splittings alone: +-{np.std(ops):.2f} deg; from m1 (0 to 1 meV): {op(nu(0)):.2f} to {op(nu(1e-3)):.2f}")
# classification chance
P("\nclassification: both coloured families above 45, the colourless ones at or below.")
P("  chance that the two quark families are the top two of four by luck: 1 in 6 (17%). Suggestive, not evidence.")
P("  within each weak pair the UPPER partner sits farther from equilibrium: up +7.2 vs down +3.1; neutrino -7.9 vs charged lepton 0. Chance of that in both pairs: 1 in 4.")
# 'energy' reading: distance vs the heaviest member's Higgs link
P("\n'distance = energy away from equilibrium' against the heaviest member's Higgs link (weight / 174 GeV):")
for nm,m3,d in (('up',168.26,7.21),('down',2.839,3.12),('charged lepton',1.777,0.0),('neutrino',5.0e-11,-7.86)):
    P(f"  {nm:<15} link {m3/174.1:.2e}  distance {d:+.2f}")
P("  order is kept (stronger link = higher), but down and the tau have almost the same link (0.016 vs 0.010) and sit 3.1 deg apart: the link's size alone does not set the distance.")
open('distance_from_equilibrium_2026-11-16.txt','w').write('\n'.join(out)+'\n')
