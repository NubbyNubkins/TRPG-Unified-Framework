"""Idea tested: the pulls are opposing but not opposite, which is why a twist appears.
Geometry: opening = angle of (sqrt m) from the cube's body diagonal; TWIST = the turn of (sqrt m) AROUND the body diagonal,
measured from the nearest edge direction (0 deg = pointing at an edge: two lighter weights equal; 60 deg = pointing at a
face diagonal: two heavier weights equal). Leptons: 2/9 rad (Koide/Brannen form sqrt m_i ~ 1 + sqrt2 cos(twist + 2 pi i/3))."""
import numpy as np, math
out=[]
def P(*x):
    s=' '.join(map(str,x)); print(s); out.append(s)
def geo(m):
    r=np.sqrt(np.sort(np.asarray(m,float))); n=np.ones(3)/math.sqrt(3)
    par=r@n; perp=r-par*n; op=math.degrees(math.atan2(np.linalg.norm(perp),par))
    e=np.array([0,0,1.0]); ep=e-(e@n)*n; ep/=np.linalg.norm(ep)
    f=np.cross(n,ep); tw=math.degrees(math.atan2(perp@f,perp@ep))
    return op,abs(tw)
def nu(m1): return [m1,math.sqrt(m1**2+7.42e-5),math.sqrt(m1**2+2.51e-3)]
fam=[('up quarks, M_Z',[1.23,620,168260]),('up quarks, own scale',[2.16,1273,172570]),('down quarks, M_Z',[2.67,53.2,2839]),('down quarks, own scale',[4.70,93.5,4183]),
     ('charged leptons',[0.51099895,105.658,1776.93]),('neutrinos, m1 0.36 meV',nu(3.6e-4)),('neutrinos, m1 0',nu(0)),('hidden partners N',[7.3e5,1.3e9,6.3e10])]
P("family: distance from 45 (opening - 45), twist around the diagonal (deg; 0 = at an edge line, 60 = at a face-diagonal line)")
rows={}
for nm,m in fam:
    o,t=geo(m); rows[nm]=(o-45,t); P(f"  {nm:<26} distance {o-45:+6.2f}   twist {t:6.2f} deg = {math.radians(t):.4f} rad")
P(f"\n  lepton twist check: {math.radians(rows['charged leptons'][1]):.5f} rad vs 2/9 = {2/9:.5f}")
# sideways displacement: twist x sin(opening) = how far off the straight edge-to-diagonal line the family sits
P("\nsideways offset (how far each family sits off the straight line from the diagonal to the edge), deg of arc = sin(opening) x twist:")
for nm,m in fam:
    o,t=geo(m); P(f"  {nm:<26} {math.sin(math.radians(o))*t:6.2f}")
# two pulls: one along the line (toward edge / toward face diagonal), one across (the twist). ratio across/along
P("\nratio twist / |distance| (0 when the two pulls are exactly opposite, i.e. on one line):")
for nm in rows:
    d,t=rows[nm]
    P(f"  {nm:<26} "+("on the equilibrium (distance 0): pure twist" if abs(d)<0.01 else f"{t/abs(d):.2f}"))
open('twist_as_azimuth_2026-11-16.txt','w').write('\n'.join(out)+'\n')
