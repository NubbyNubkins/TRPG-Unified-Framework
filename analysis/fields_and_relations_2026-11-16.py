# (A) other fields of physics from the same few numbers; (B) pairs/triples/quads search with the new far side
import numpy as np, math, itertools
me=0.51099895e-3; mp=0.93827; a=1/137.036; MPl=1.22089e19   # GeV
GeV_to_eV=1e9; Msun_GeV=1.989e30/1.7827e-27
print("=== A. Other fields from electron, proton, electric strength, grain ===")
tests=[("Rydberg (atom binding)", a*a*me/2*1e9, 13.606, "eV", "Bohr 1913"),
 ("Molecular vibration step (H2)", a*a*me/2*1e9*math.sqrt(me/mp)*2, 0.516, "eV", "Born-Oppenheimer 1927 (x2 for H2 reduced-mass/bond)"),
 ("Molecular rotation step (H2)", a*a*me/2*1e9*(me/mp)*2, 0.0151, "eV", "Born-Oppenheimer 1927"),
 ("Chandrasekhar limit (white dwarf, mu_e=2)", 3.098/4*(MPl**3/mp**2)/Msun_GeV, 1.44, "Suns", "Chandrasekhar 1931"),
 ("Largest cold planet ~ alpha^1.5 (M_Pl/m_p)^3 m_p", a**1.5*(MPl/mp)**3*mp/Msun_GeV/9.55e-4, 1.0, "Jupiters (order)", "Press & Lightman 1983"),
 ("Quark-soup melting ~ pion weight", 0.1396*1e3, 156.5, "MeV", "lattice QCD (HotQCD 2019)"),
 ("Hydrogen ionisation temperature ~ Rydberg/ (k ln(...)) ~ alpha^2 me/40", a*a*me/2*1e9/40*11604, 3740, "K (recombination ~3000-4000)", "Saha"),
]
for n,v,m,u,src in tests: print(f"  {n:52s} {v:12.4g} vs {m:<8g} {u:28s} ratio {v/m:6.3f}  [{src}]")
print("=== B. Relations among pairs, triples, quads (log10 GeV), with look-elsewhere control ===")
M={"nu1":3.6e-13,"nu2":8.6e-12,"nu3":5.0e-11,"e":me,"mu":0.10566,"tau":1.77686,"u":2.16e-3,"d":4.70e-3,"s":0.0935,"c":1.27,"b":4.18,"t":172.6,
   "W":80.37,"Z":91.19,"H":125.2,"v":246.2,"p":mp,"majoron":8.8e-6,"N1":7.3e5,"N2":1.3e9,"N3":6.3e10,"MMG":2.2e8,"D":4.56e13,"grain":8.63e18,"DE":2.24e-12}
K=list(M); L=np.array([math.log10(M[k]) for k in K])
built={frozenset(x) for x in [("nu1","e","N1"),("nu2","mu","N2"),("nu3","tau","N3"),("MMG","D","grain"),("W","Z","v"),("H","v","t")]}
def triples(L,tol):
    out=[]
    for i,j,k in itertools.permutations(range(len(L)),3):
        if i<j and abs(L[i]+L[j]-2*L[k])<tol: out.append((i,j,k))   # k = geometric mean of i,j (a mirror/plank)
    return out
def quads(L,tol):
    out=[]
    for i,j,k,l in itertools.combinations(range(len(L)),4):
        for (p,q,r,s) in ((i,j,k,l),(i,k,j,l),(i,l,j,k)):
            if abs(L[p]+L[q]-L[r]-L[s])<tol: out.append((p,q,r,s))   # a*b = c*d (balanced pairs)
    return out
rng=np.random.default_rng(7)
for tol in (0.02,):
    T=triples(L,tol); Q=quads(L,tol)
    simT=[len(triples(rng.uniform(L.min(),L.max(),len(L)),tol)) for _ in range(300)]
    simQ=[len(quads(rng.uniform(L.min(),L.max(),len(L)),tol)) for _ in range(150)]
    print(f" tolerance {tol} decades (~{(10**tol-1)*100:.0f}%):")
    print(f"  mirrors (k = geometric mean of i,j): {len(T)} found vs {np.mean(simT):.1f} random (chance >= {np.mean(np.array(simT)>=len(T)):.2f})")
    for i,j,k in T:
        tag=" [built in]" if frozenset((K[i],K[j],K[k])) in built else ""
        print(f"    {K[k]} = sqrt({K[i]} x {K[j]})  off {10**(L[i]+L[j]-2*L[k])**0.5 if False else (10**((L[i]+L[j])/2-L[k])-1)*100:+.1f}%{tag}")
    print(f"  balanced pairs a x b = c x d: {len(Q)} found vs {np.mean(simQ):.1f} random (chance >= {np.mean(np.array(simQ)>=len(Q)):.2f})")
