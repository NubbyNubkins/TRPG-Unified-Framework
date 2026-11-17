# Can the one-handed leftover from the hidden partners' fall stir up
# the magnetic field the lumps need? Order-of-magnitude checks.
import math
MPl=1.22e19; g=106.75; v=174.0; aY=0.0102; eta_obs=6.1e-10
print("1. Published falsifier: does a one-handed 0.05 nG field make too much matter at the 160 GeV step?")
B_need=5e-11
for B_ok in (1e-17,1e-16):
    print(f"   Kamada & Long 2016: right matter excess at B0={B_ok:.0e} G; excess grows as B0^3 on the endpoint line;"
          f" 0.05 nG overshoots by {(B_need/B_ok)**3:.1e}")
print("   Boyer & Neronov 2025 (redo, slower field decay): strong fields give excess = 0.6/(grain x T), independent of strength;")
lam=0.6/eta_obs
H_radius=MPl/(1.66*math.sqrt(g)*130.0)  # in units of 1/T at T=130 GeV
print(f"   observed excess needs grain = {lam:.1e}/T = {lam/H_radius:.1e} of the horizon at 130 GeV; relic today 0.2 to 100 pG (needed: 50 pG)")
print("   So the two published calculations disagree; the idea is not dead in general.\n")
print("2. Can TRPG's own one-handed leftover twist up the field? (Joyce & Shaposhnikov 1997 mechanism)")
T=8e4  # below this the electron's handedness is no longer kept
HoverT=1.66*math.sqrt(g)*T/MPl
for x in (1e-8,1e-5,1e-3,1e-2):
    rate=(aY*x/(2*math.pi))**2/100.0
    print(f"   imbalance {x:.0e} of the light: twisting rate / growth of space at 80 TeV = {rate/HoverT:.1e}")
x_need=math.sqrt(HoverT*100)/(aY/(2*math.pi))
print(f"   imbalance needed for the twisting to keep up: {x_need:.1e}\n")
print("3. Largest imbalance the hidden partners can make (Davidson & Ibarra 2002 cap, m3 = 0.050 eV):")
for M in (1e6,2.41e8):
    eps=3*M*0.050e-9/(16*math.pi*v*v)
    print(f"   partner weight {M:.2e} GeV: at most {eps:.1e}; short of {x_need:.0e} by {x_need/eps:.1e}")
M_need=x_need*16*math.pi*v*v/(3*0.050e-9)
print(f"   partner weight that could do it: {M_need:.1e} GeV (start band tops out at 2.41e8)")
