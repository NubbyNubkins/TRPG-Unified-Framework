# Can a stable push/pull balance (energy valley) fix D's reach y = link/own ~ 34.5?
# Press: 1-kappa = 1/(1+y^2). Profile chi2(1-kappa) from reverse_engineer_far_side_2026-11-16.txt (2 dof).
import numpy as np, itertools
prof=np.array([[2.00e-4,1103.21],[2.51e-4,1090.62],[3.16e-4,905.24],[3.98e-4,309.79],[5.01e-4,27.64],[6.31e-4,46.01],
 [7.94e-4,3.41],[1.00e-3,9.51],[1.26e-3,21.00],[1.58e-3,34.25],[2.00e-3,86.64],[2.51e-3,122.92],[3.16e-3,128.12],[3.98e-3,111.36]])
def chi2(omk):
    if omk<prof[0,0] or omk>prof[-1,0]: return float('inf')
    return float(np.interp(np.log(omk),np.log(prof[:,0]),prof[:,1]))
def yband():
    lo,hi=7.82e-4,9.22e-4; return np.sqrt(1/hi-1),np.sqrt(1/lo-1)
print("band for y (chi2 within 4 of best):",[round(v,2) for v in yband()])
def show(name,y):
    omk=1/(1+y*y); c=chi2(omk); print(f"  {name:55s} y={y:8.3f} 1-k={omk:.3e} chi2={c:9.2f} {'IN BAND' if c<=3.41+4 else 'FAIL'}")
print("\n1. Even balances (no scales)")
show("link = own (even 1:1)",1.0); show("rule-7 balance link/own = 1/sqrt2",1/np.sqrt(2)); show("rule-7 inverted sqrt2",np.sqrt(2))
# 2. power-law valley E = A r^a - B r^-b between two far-side scales: r* ~ ((b/a) R)^(1/(a+b))
F=4.56e13; MMG=2.2e8; G=8.63e18; Rs={"F/MMG":F/MMG,"grain/F":G/F,"grain/MMG":G/MMG}
print("\n2. Power-law valley (push r^a vs pull r^-b, a,b=1..4) with each pair of far-side scales")
hits=0; tot=0
for (rn,R),a,b in itertools.product(Rs.items(),range(1,5),range(1,5)):
    y=((b/a)*R)**(1/(a+b)); tot+=1; c=chi2(1/(1+y*y))
    if c<=7.41: hits+=1; print(f"  in band: {rn}, a={a}, b={b}: y={y:.2f}, chi2={c:.2f}")
print(f"  {hits} of {tot} simple valleys land in band")
for rn,R in Rs.items():
    print(f"  exponent needed for {rn}: {np.log(34.5)/np.log(R):.4f}")
# 3. Coleman-Weinberg valley in the mixing: heavy mass^2 = e*mu^2 -> 1+y^2 = e (mu/own)^2
print("\n3. One-loop (Coleman-Weinberg) valley: heavy D mass = sqrt(e)*mu")
for mn,mu in {"mu = F (own=F)":F,"mu = grain":G,"mu = MMG":MMG}.items():
    y2=np.e*(mu/F)**2-1; print(f"  {mn:20s}", "no valley (y^2<0)" if y2<=0 else ""); 
    if y2>0: show(mn,np.sqrt(y2))
print("  the valley exists for any y; mu needed for y=34.5:",f"{34.5/np.sqrt(np.e)*F:.2e} GeV (no scale of the shape sits there)")
# 4. Luck: fraction of log-uniform y in [1,1000] landing in band
lo,hi=yband(); print(f"\n4. Band width: a random log-uniform y in 1..1000 lands in band {np.log(hi/lo)/np.log(1000)*100:.1f}% of the time")
print("   simple constants inside band (NOT evidence):", {k:round(v,2) for k,v in {"4*pi*e":4*np.pi*np.e,"11*pi":11*np.pi,"1/alpha_meeting^(..)":39.0}.items() if lo<=v<=hi})
