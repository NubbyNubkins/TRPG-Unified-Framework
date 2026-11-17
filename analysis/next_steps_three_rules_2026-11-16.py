# Next steps on the three far-side rules
import numpy as np, math, itertools
# --- octonions (Cayley-Dickson) ---
def cd(a,b):
    n=len(a)
    if n==1: return a*b
    h=n//2; a1,a2,b1,b2=a[:h],a[h:],b[:h],b[h:]
    conj=lambda x: np.concatenate([x[:1],-x[1:]])
    return np.concatenate([cd(a1,b1)-cd(conj(b2),a2), cd(b2,a1)+cd(a2,conj(b1))])
E=np.eye(8); u=E[1]
print("=== Rule 1: the facing u acts on the far side by commutator [u,x] ===")
for k in range(8):
    c=cd(u,E[k])-cd(E[k],u); print(f"  e{k}: |[u,e{k}]| = {np.linalg.norm(c):.0f}")
x=np.random.default_rng(3).normal(size=8)
xc=x.copy(); xc[:2]=0
print(f"  random x: |[u,x]|^2 = {np.linalg.norm(cd(u,x)-cd(x,u))**2:.4f} = 4 x |coloured part|^2 = {4*np.linalg.norm(xc)**2:.4f}")
print("  -> zero on exactly the colourless parts {1,u} (the MMG parts, plus the real diagonal b,c), equal non-zero on all 6 coloured parts (D).")
print("  Exact (octonion algebra; colour = stabiliser of u, Gunaydin-Gursey 1973). Same pattern as matter number (D +-2/3, MMG 0).")
print("  The weak force forbids any term weighing only some MMG parts, so the only far-side weights are |D|^2 and |MMG|^2.")
print("  The spot is a STATE (the hidden partner, an earlier test), not a setting: the determinant term a(bc-|x|^2) gives no weight -> |MMG|^2 coefficient 0.")
print("  Result: D weighs F (facing term), MMG only via the grain: F^2/grain = %.3e GeV." % (4.56e13**2/8.63e18))
print("\n=== Rule 2: does the same setting link D to the down quarks? ===")
print("  The facing commutator maps the 10 into itself; D-to-down-single (10 to 16) needs a spinor-direction value (E6 cubic has 10.16.16).")
F=4.56e13; G=8.63e18; N3=6.3e10
for lam in (1.0, math.sqrt(4*math.pi/38.6)):
    v16=34.5*F/lam
    print(f"  coupling {lam:.2f}: spinor value needed {v16:.2e} GeV; it would also weigh a hidden partner ~ v16^2/grain = {v16**2/G:.2e} GeV (N3 = {N3:.1e}, ratio {v16**2/G/N3:.1f})")
print("  -> the link must come from the hidden-partner (matter-number-breaking) direction, not the facing: same SETTING fails; 34.5 still a ratio of two values. HINT: v16^2/grain lands within a factor 2-5 of N3.")
print("\n=== Rule 3: angle vs partner type ===")
def op(m):
    m=np.array(m,float); s=np.sqrt(m); return math.degrees(math.acos(s.sum()/math.sqrt(3*m.sum())))
print("  mass-term types (Georgi-Glashow): up = pair.pair (10.10); down and charged leptons = pair.single (10.5bar); neutrinos = single.singlet + seesaw")
dn_base=[15.2,111.3,8347.8]; dn=[2.678,53.366,2872.693]
print(f"  prediction if 'down = lepton term + D press': unpressed down opening = 45.000")
print(f"  the unpressed down base {dn_base} MeV -> opening {op(dn_base):.2f} deg (pressed {op(dn):.2f}); leptons 45.00, up {op([1.242,627.8,162500]):.2f}")
