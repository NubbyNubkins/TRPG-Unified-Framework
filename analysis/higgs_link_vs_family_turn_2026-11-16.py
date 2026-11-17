# Does the Higgs link (top's own pull) or the family turn set the places along the side?
# One-loop running (Machacek & Vaughn 1984), third-family links only.
import numpy as np
MZ, F = 91.19, 4.56e13
g1,g2,g3 = 0.4615,0.6519,1.2177          # GUT-normalised g1 at M_Z
yt,yb = 0.968, 0.0164
rc_t, rs_b = 0.619/171.7, 0.0550/2.86    # link2/link3 at M_Z, up ladder and down ladder
L = np.log(F/MZ); n=20000; dt=L/n; k=1/(16*np.pi**2)
lu=np.log(rc_t); ld=np.log(rs_b)
for _ in range(n):
    bt = yt*k*(4.5*yt**2+1.5*yb**2-8*g3**2-2.25*g2**2-17/20*g1**2)
    bb = yb*k*(4.5*yb**2+1.5*yt**2-8*g3**2-2.25*g2**2-1/4*g1**2)
    lu += -1.5*k*(yt**2-yb**2)*dt      # ln(y_c/y_t)
    ld += -1.5*k*(yb**2-yt**2)*dt      # ln(y_s/y_b)  : opposite sign
    g1 += k*41/10*g1**3*dt; g2 += -k*19/6*g2**3*dt; g3 += -k*7*g3**3*dt
    yt += bt*dt; yb += bb*dt
tu0,td0 = np.sqrt(rc_t),np.sqrt(rs_b); tu1,td1 = np.exp(lu/2),np.exp(ld/2)
print(f"top link at the meeting: {yt:.3f}")
print(f"up place   t: {tu0:.4f} at M_Z -> {tu1:.4f} at the meeting ({(tu1/tu0-1)*100:+.1f}%)")
print(f"down place t: {td0:.4f} at M_Z -> {td1:.4f} at the meeting ({(td1/td0-1)*100:+.1f}%)")
print(f"down/up step: {td0/tu0:.3f} at M_Z -> {td1/tu1:.3f} at the meeting; needed to close the gap: 1.000")
print("Family turn 2-3 (|Vcb| measured 0.0411) vs difference of places:")
print(f"  M_Z: {td0-tu0:.4f} to {td0+tu0:.4f};  meeting: {td1-tu1:.4f} to {td1+tu1:.4f}  (Vcb at meeting ~ {0.0411*np.exp(-(lu-np.log(rc_t))/1)**-0.5:.4f})")
print("Photons per unit of matter today: %.2e (eta = 6.1e-10)" % (1/6.1e-10))
