# Can the theory's own start, or a later step, leave the weak magnetic field the lumps need (about 0.05 billionths of a gauss today)?
# Order-of-magnitude scaling only. A tangled field made at temperature T has a grain no larger than the horizon then, and it decays as it untangles:
#   untwisted field: B^4 L^5 stays fixed (Hosking & Schekochihin 2023);  fully one-handed (twisted) field: B^2 L stays fixed (helicity; Frisch et al. 1975)
#   it ends, when atoms form, on the line B = 1e-8 gauss x (L / Mpc) (Banerjee & Jedamzik 2004)
import numpy as np
T0=2.349e-4; c=2.998e10; Mpc=3.0857e24
def Li(T,g): t=2.42/np.sqrt(g)/(T/1e6)**2; return 2*c*t*(T/T0)*(g/3.91)**(1/3)/Mpc   # horizon then, stretched to today, Mpc
Brad=3.24e-6    # gauss: field holding ALL of today's light energy
need=5e-11
print("needed today: 0.05 billionths of a gauss (Jedamzik & Pogosian 2020)\n")
for lab,T,g in (("the light begins, bottom of the band (1e6 GeV)",1e15,106.75),("the light begins, top of the band (2.4e8 GeV)",2.4e17,106.75),("the surface's resistance switches on (160 GeV)",160e9,106.75),("quarks lock into protons (155 MeV)",155e6,17.25)):
    print(lab)
    for share in (1.0,0.01):
        Bi=Brad*np.sqrt(share)*(3.91/g)**(1/6); L0=Li(T,g)
        Ln=(Bi*L0**1.25/1e-8)**(4/9); Bn=1e-8*Ln          # untwisted
        Lh=(Bi*L0**0.5/1e-8)**(2/3);  Bh=1e-8*Lh          # fully twisted
        print(f"   share of the energy in the field {100*share:3.0f}%: untwisted ends at {Bn/1e-9:.1e} billionths (short by x{need/Bn:,.0f}); fully twisted ends at {Bh/1e-9:.1e} billionths ({'short by x%.0f'%(need/Bh) if Bh<need else 'ENOUGH, x%.1f over'%(Bh/need)}), grain {Lh*3.262e6:,.0f} light-years")
