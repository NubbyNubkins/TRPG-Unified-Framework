"""Home-condition stability (14 Nov 2026): lifetime measured on the particle's own clock, and where it is stable.
Q = m c^2 / Gamma = (lifetime) x (m c^2 / hbar): how many times the particle's own wave 'rings' before it fades.
Data PDG 2024. Adds sheet Home_Stability to TRPG_Quark_Gluon_Matrices.xlsx."""
import math, pandas as pd
F = "TRPG_Quark_Gluon_Matrices.xlsx"
hbar = 6.582119569e-22  # MeV s
rows = [  # name, mass MeV, lifetime s (or None), width MeV (or None), home where it is stable / in balance, what opens the door here
 ("proton", 938.272, None, 0, "everyday matter", "nothing here; but inside a neutron star p + e -> n, so there it turns into neutrons"),
 ("electron", 0.511, None, 0, "everywhere", "nothing (lightest charged particle)"),
 ("neutron", 939.565, 878.4, None, "inside nuclei and neutron stars: STABLE there", "free in space: n -> p e nu is an open door (1.29 MeV to spare)"),
 ("muon", 105.658, 2.1969811e-6, None, "neutron-star cores (electron energy > 105.7 MeV blocks the decay); early universe at ~100 MeV", "empty space: mu -> e nu nu is open"),
 ("charged pion", 139.570, 2.6033e-8, None, "early universe ~150 MeV (made as fast as it decays); proposed pion condensate in dense matter", "empty space: pi -> mu nu"),
 ("neutral pion", 134.977, 8.43e-17, None, "early universe ~150 MeV", "pi0 -> 2 photons"),
 ("kaon (charged)", 493.677, 1.2380e-8, None, "early universe; kaon condensate proposed in neutron stars", "K -> mu nu"),
 ("Lambda (strange baryon)", 1115.683, 2.617e-10, None, "neutron-star cores above ~2-3x nuclear density (hyperons expected stable there)", "Lambda -> p pi"),
 ("Sigma0", 1192.642, 7.4e-20, None, "early universe", "Sigma0 -> Lambda photon"),
 ("Omega-", 1672.45, 8.21e-11, None, "quark-gluon plasma / early universe", "Omega -> Lambda K"),
 ("D0 (charm meson)", 1864.84, 4.103e-13, None, "early universe ~1 GeV", "c -> s W"),
 ("B0 (bottom meson)", 5279.66, 1.519e-12, None, "early universe ~5 GeV", "b -> c W"),
 ("J/psi", 3096.9, None, 0.0926, "quark-gluon plasma melts it; early universe ~3 GeV", "c c-bar -> 3 gluons"),
 ("rho", 775.26, None, 149.1, "inside hadron gas at ~150 MeV", "rho -> pi pi"),
 ("Delta", 1232.0, None, 117.0, "hot nuclear matter (heavy-ion collisions); early universe", "Delta -> N pi"),
 ("tau", 1776.93, 2.903e-13, None, "early universe ~2 GeV", "tau -> nu + W*"),
 ("W", 80369.2, None, 2085.0, "early universe above ~100 GeV (before the Higgs switched on, massless and stable as a carrier)", "W -> f f'"),
 ("Z", 91188.0, None, 2495.5, "same as W", "Z -> f f-bar"),
 ("Higgs", 125200.0, None, 4.1, "the Higgs FIELD is always 'on' everywhere; the particle is a ripple of it", "H -> b b-bar"),
 ("top quark", 172570.0, None, 1.42, "early universe above ~170 GeV", "t -> W b"),
]
out = []
for n, m, tau, G, home, door in rows:
    if G == 0: Q = float("inf"); tau_s = "stable"
    else:
        if G is None: G = hbar / tau
        else: tau = hbar / G
        Q = m / G; tau_s = tau
    out.append({"particle": n, "mass (MeV)": m, "lifetime on OUR clock (s)": tau_s,
                "own clock: rings before fading Q = m/Gamma": Q if Q != float("inf") else "forever",
                "log10 Q": round(math.log10(Q), 2) if Q != float("inf") else "inf",
                "home temperature ~ m/20 (MeV): hotter than this, made as fast as it decays": round(m / 20, 3),
                "home where it is stable or in balance": home, "door that is open in our everyday conditions": door})
df = pd.DataFrame(out)
with pd.ExcelWriter(F, mode="a", engine="openpyxl", if_sheet_exists="replace") as w:
    df.to_excel(w, sheet_name="Home_Stability", index=False)
print(df[["particle", "log10 Q", "home temperature ~ m/20 (MeV): hotter than this, made as fast as it decays"]].to_string(index=False))
