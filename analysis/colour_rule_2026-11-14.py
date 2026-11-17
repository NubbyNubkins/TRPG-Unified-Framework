"""Colour rule test (14 Nov 2026): can the bottom-mass input be removed, and do quarks share the 2/9 twist?
Idea tested: a coloured family carries the lepton family's relationship matrix once per colour.
Energy and phase both ADD, so scale -> 3 x scale and twist -> 3 x twist (3 = number of colours).
The numerical observation M_q = 3 M_l, delta_q = 3 delta_l for (s, c, b) is Rivero's (arXiv:1111.7232, 2011);
he gave no reason for the 3. Reading it as the colour count, and the checks below, are what is added here.
Credit: Koide 1981; Foot 1994; Brannen 2006; Rivero & Gsponer 2005; Rivero 2011, 2026 (arXiv:2606.10060);
Gatto, Sartori & Tonin 1968; Leutwyler 1996; Colangelo, Lanz, Leutwyler & Passemar 2018; Shulga 2026 (arXiv:2605.10245);
Pendleton & Ross 1981; Hill 1981; PDG 2024."""
import math, numpy as np
rng = np.random.default_rng(33)
me, mmu, mtau, dtau = 0.51099895, 105.6583755, 1776.93, 0.09
PDG = dict(s=(93.5, 0.8), c=(1273.0, 4.6), b=(4183.0, 7.0), t=(172570.0, 290.0), d=(4.70, 0.07), u=(2.16, 0.07))
R2 = math.sqrt(2)

def params(roots):
    """signed sqrt-masses -> (scale M, link strength r, twist delta). r = sqrt2 means 45 degrees (Q = 2/3)."""
    s = np.array(roots, float); a = s.mean(); x = s/a - 1
    r = math.sqrt((x**2).sum()*2/3); return a*a, r, math.acos(max(x)/r)
def family(M, delta, r=R2):
    return [math.sqrt(M)*(1 + r*math.cos(delta + 2*math.pi*k/3)) for k in range(3)]
def koide_third(m1, m2, heavy=True):
    a, b = math.sqrt(m1), math.sqrt(m2); B = -4*(a+b); C = 3*(m1+m2) - 2*(a+b)**2
    return ((-B + (1 if heavy else -1)*math.sqrt(B*B - 4*C))/2)**2
def pr(*a): print(*a)

pr("="*86); pr("1. The lepton family (inputs: electron and muon; everything else below follows)"); pr("="*86)
Ml, rl, dl = params([math.sqrt(x) for x in (me, mmu, mtau)])
dl_hi = params([math.sqrt(x) for x in (me, mmu, mtau+dtau)])[2]; ddl = abs(dl_hi - dl)
pr(f"  scale M_l = {Ml:.4f} MeV; link strength r = {rl:.6f} (sqrt2 = {R2:.6f}); twist = {dl:.7f} +- {ddl:.7f}  (2/9 = {2/9:.7f})")
pr(f"  3 x M_l = {3*Ml:.3f} MeV  (the record's M2 step already noted this is 0.35% from the proton, 938.272)")

pr("\n"+"="*86); pr("2. The colour rule: scale x3, twist x3  ->  (s, c, b), then top, down, up.  NO quark mass is put in."); pr("="*86)
roots = family(3*Ml, 3*(2/9)); s_p, c_p, b_p = sorted(x*x for x in roots)
t_p = koide_third(c_p, b_p); Vus = 0.22431; d_p = s_p*Vus**2
rsd = 1/Vus**2; u_p = {Q: d_p*math.sqrt(1 - (rsd/Q)**2) for Q in (21.4, 22.1, 22.8)}
pred = dict(t=t_p, b=b_p, c=c_p, s=s_p, d=d_p, u=u_p[22.1])
how = dict(t="45-degree balance of (c, b, t)", b="colour rule", c="colour rule", s="colour rule (its root is negative)",
           d="Cabibbo cell: m_d = m_s |V_us|^2", u="eta -> 3 pi cell (Q = 22.1 +- 0.7)")
pr(f"  signed roots of the tripled family: {[round(x, 3) for x in roots]}")
pr(f"  {'quark':7s}{'how':38s}{'predicted':>12s}{'measured':>20s}{'off':>8s}{'sigma':>7s}")
for q in "tbcsdu":
    m, e = PDG[q]; p = pred[q]
    pr(f"  {q:7s}{how[q]:38s}{p:12.2f}{m:12.2f} +- {e:<5.2f}{100*(p-m)/m:+7.2f}%{(p-m)/e:+7.1f}")
pr(f"  up-quark range from the eta data alone: {u_p[21.4]:.2f} to {u_p[22.8]:.2f} MeV")
pr(f"  cross-check on the top from a different rule (full Higgs link, y_t = 1): {246219.65/R2/1000:.2f} GeV vs {t_p/1000:.2f} GeV here: {100*(246219.65/R2-t_p)/t_p:+.2f}% apart")

pr("\n"+"="*86); pr("3. What the measured (s, c, b) say on their own, with errors"); pr("="*86)
vals = []
for _ in range(20000):
    s, c, b = (rng.normal(*PDG[q]) for q in "scb"); vals.append(params([-math.sqrt(s), math.sqrt(c), math.sqrt(b)]))
vals = np.array(vals); Mq, rq, dq = vals.mean(0); eM, er, ed = vals.std(0)
pr(f"  scale  M_q = {Mq:7.2f} +- {eM:.2f} MeV   -> M_q / M_l = {Mq/Ml:.3f}   (rule: 3)")
pr(f"  twist      = {dq:7.4f} +- {ed:.4f}       -> twist / (2/9) = {dq/(2/9):.3f}  (rule: 3)")
pr(f"  strength r = {rq:7.4f} +- {er:.4f}       (45-degree value: {R2:.4f})")
pr("  All three are 1-3% from the rule, and the whole offset is one number, the charm mass:")
c_from_sb = koide_third(PDG['s'][0], PDG['b'][0], heavy=False)   # needs the signed version; do it directly
a = -math.sqrt(93.5) + math.sqrt(4183.0); disc = (4*a)**2 - 4*(3*(93.5+4183.0) - 2*a*a); c_sb = ((4*a - math.sqrt(disc))/2)**2
c_bt = koide_third(4183.0, 172570.0, heavy=False)
pr(f"  charm wanted by three separate relations: from (b, t) {c_bt:.0f} MeV; from (s, b) {c_sb:.0f} MeV; from leptons x3 {c_p:.0f} MeV")
pr(f"  They agree with each other within 1.5%, and all sit 6-8% above the measured m_c(m_c) = 1273.0 +- 4.6 MeV.")
pr("  Charm's mass depends strongly on its definition: 0.99 GeV (at 3 GeV), 1.27 (own scale), 1.55 (half the J/psi), 1.67 (pole).")

pr("\n"+"="*86); pr("4. Independent support that the bottom is tied to the lighter quarks (Rivero 2026)"); pr("="*86)
b_r = (3*math.sqrt(93.5) + math.sqrt(1273.0))**2; eb = 2*math.sqrt(b_r)*math.hypot(3*0.8/(2*math.sqrt(93.5)), 4.6/(2*math.sqrt(1273.0)))
pr(f"  sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)  ->  m_b = {b_r:.1f} +- {eb:.0f} MeV vs measured 4183 +- 7   ({100*(b_r-4183)/4183:+.2f}%)")
pr("  This uses the measured strange and charm, so it does not remove an input, but it shows the bottom is not a free number.")

pr("\n"+"="*86); pr("5. Odds that the colour rule's match is luck (look-elsewhere included)"); pr("="*86)
rules = [(k, j) for k in range(1, 10) for j in range(1, 10)]
preds = np.array([sorted(x*x for x in family(k*Ml, j*dl)) for k, j in rules])
def n_match(tr, tol_tight=0.004, tol_loose=0.015):
    rel = np.abs(preds/np.array(tr) - 1)             # (rules, 3)
    ok_loose = rel < tol_loose; ok_tight = rel < tol_tight
    return int(((ok_loose.sum(1) >= 2) & (ok_tight.sum(1) >= 1)).sum())
real = sorted(PDG[q][0] for q in "scb")
pr(f"  rules tried: scale x k, twist x j, for k, j = 1..9 (81 rules). Real (s, c, b) matched by {n_match(real)} of them (the k = j = 3 one).")
import itertools
quarks = dict(u=2.16, d=4.70, s=93.5, c=1273.0, b=4183.0, t=172570.0)
hits = [(tr, n_match(sorted(quarks[q] for q in tr))) for tr in itertools.combinations("udscbt", 3)]
pr(f"  over all 20 quark triples, triples matched by any rule: {[''.join(t) for t, n in hits if n]}")
N = 200000; lo, hi = math.log(1.0), math.log(2e5); cnt = 0
for _ in range(N):
    tr = np.sort(np.exp(rng.uniform(lo, hi, 3)))
    if n_match(tr): cnt += 1
p1 = cnt/N
pr(f"  random triples (log-uniform, 1 MeV to 200 GeV) matched this well by ANY of the 81 rules: {cnt} in {N} = {p1:.1e}")
pr(f"  times 20 triples to choose from: about {20*p1:.1e}, i.e. roughly 1 chance in {1/(20*p1):.0f}.")
same = [i for i, (k, j) in enumerate(rules) if k == j]; only3 = [i for i, (k, j) in enumerate(rules) if k == j == 3]
N2 = 2000000; trs = np.sort(np.exp(rng.uniform(lo, hi, (N2, 3))), axis=1)
def frac(idx):
    hit = np.zeros(N2, bool)
    for i in idx:
        rel = np.abs(preds[i][None, :]/trs - 1); hit |= ((rel < 0.015).sum(1) >= 2) & ((rel < 0.004).sum(1) >= 1)
    return hit.mean()
f_same, f_3 = frac(same), frac(only3)
pr(f"  stricter counting: same multiplier for scale and twist (9 rules): 1 chance in {1/(20*f_same):.0f}; only 'x3 = colours' counted: 1 chance in {1/(20*f_3):.0f}.")
pr("  So between about 1 in 25 and 1 in 2,000 depending on how many rules one says would have counted. Suggestive, not proof.")

pr("\n"+"="*86); pr("6. The twist 2/9 as a ring number (three placements a third of a turn apart)"); pr("="*86)
n = np.arange(1, 2000001); gsum = ((1 - np.cos(2*math.pi*n/3))/n**2).sum()/math.pi**2
pr(f"  response of a ring between two seats a third of a turn apart: (1/pi^2) sum (1 - cos(2 pi n/3))/n^2 = {gsum:.7f}  (2/9 = {2/9:.7f})")
pr("  same thing in plain terms: two seats split the ring into arcs of 1/3 and 2/3; the two paths in parallel give")
pr(f"  (1/3 x 2/3)/(1/3 + 2/3) = {1/3*2/3:.7f} = 2/9.  For N seats: (1/N)(1 - 1/N).")
pr(f"  measured lepton twist / (2/9) = {dl/(2/9):.7f} +- {ddl/(2/9):.7f}; so the unit ('one full ring = one radian') holds to about 1 part in 10^5,")
pr("  but that unit is NOT derived (Shulga 2026 calls the same normalization an undetermined matching coefficient).")
pr(f"  quark triple: three colour links in a row add: 3 x 2/9 = {3*2/9:.4f}; measured {dq:.4f} +- {ed:.4f} (the charm offset again).")

pr("\n"+"="*86); pr("7. Why 3 generations makes the two lepton rules agree"); pr("="*86)
for Ngen in (2, 3, 4, 5):
    pr(f"  N = {Ngen}: 'self = cross' gives Q = 2/N = {2/Ngen:.4f};  ring fraction left over 1 - 1/N = {1-1/Ngen:.4f}  {'<- equal only here' if Ngen == 3 else ''}")
pr("  With 3 generations the balance number (2/3) and the ring's remaining fraction (2/3) coincide, so twist = Q/3 = 2/9 is one statement, not two.")
