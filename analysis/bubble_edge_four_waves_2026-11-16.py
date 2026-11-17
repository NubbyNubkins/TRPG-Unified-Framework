# Is "edge at 4 own waves" a rule for every bubble (4 = 3D + time), or the proton only?
hc=197.3269804
B=[("proton (charge radius, muonic H)",938.272,0.84075,0.00064),("neutron (magnetic radius)",939.565,0.864,0.009),
   ("charged pion",139.570,0.659,0.004),("charged kaon",493.677,0.560,0.031),("deuteron",1875.613,2.12799,0.00074),("helium-4 nucleus",3727.379,1.6755,0.0028)]
for n,m,r,dr in B:
    w=hc/m; x=r/w; print(f"{n:34s} radius {r:.4f} fm, own wave {w:.4f} fm, ratio {x:.3f} +- {dr/w:.3f}  ({(x-4)/(dr/w):+.1f} sigma from 4)")
