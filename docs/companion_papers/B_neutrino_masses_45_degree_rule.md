---
title: "Normal-ordered neutrino masses from the Koide relation with one negative root, and a reason for Q = 2/3"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "25 November 2026"
---

## Abstract

The charged-lepton masses satisfy Koide's relation Q = (Σm)/(Σ√m)² = 2/3 to about 10⁻⁵. Brannen proposed applying the same relation to neutrinos with one negative square root. Fixed by the two measured mass splittings, it then has a normal-ordered solution:

- with JUNO's 2026 splittings: m₁ = 0.359, m₂ = 8.60, m₃ = 50.09 meV, summing to 59.05 meV;
- with the NuFIT 6.0 global fit (Esteban et al. 2024): m₂ = 8.66 meV and a sum of 59.16–59.36 meV.

The inverted-ordering solution sums to about 0.103 eV, which today's cosmological bound disfavours. With TM1 mixing and free Majorana phases, the effective neutrinoless double-beta mass is 1.34–4.02 meV.

We give a reason for Q = 2/3 itself. If the family amplitude matrix is invariant under a 120° turn, it is circulant: a·1 + b·P + b*·P². Maximum uncertainty between the self term and the relationship terms then fixes the self share of its weight at ½, which is exactly Q = 2/3. This principle was adopted knowing that Koide's relation holds, so we state it as a choice, not a derivation. Its value is that one statement covers the charged-lepton and neutrino cases together. The predictions can fail through inverted ordering, a sum clearly different from 0.059 eV, or a lightest mass far from 0.36 meV.

## 1. The relation and its geometry

Write the square roots of the three masses as the eigenvalues of a Hermitian circulant matrix, √m~k~ = a + 2|b| cos(δ + 2πk/3). The angle between the root vector and the democratic direction (1,1,1) satisfies tan θ = √2 |b|/a. Koide's Q = 2/3 is θ = 45°, that is |b| = a/√2 (Koide 1983; Foot 1994; Brannen 2006).

For the charged leptons, the masses of e and μ predict m~τ~ = 1776.969 MeV, against 1776.93 ± 0.09 MeV (+0.4σ). The twist δ comes out as 0.222222, which matches 2/9 to about one part in a million but is not derived here. Exactly 2/9 together with exactly 45° would put m~e~ 9.8 ppm low.

## 2. Why the self share is one half

The self share of the matrix weight is p = 3a²/(3a² + 6|b|²). Every symmetric measure of uncertainty between the two parts, self and relationship, is largest at p = ½ (Shannon 1948). And p = ½ is exactly |b|² = a²/2, the 45° condition. The measured charged leptons give p = 0.50000.

For quarks, p = 0.456 (down-type) and 0.393 (up-type), and no choice of signs reaches ½. The principle is therefore stated for colour-neutral families only. Running from a high scale should spoil the lepton relation by about 0.35% (Sumino 2009). It holds for pole masses, and we do not explain why.

## 3. Neutrinos

With one negative root, the conditions Q = 2/3, m₂² − m₁² = Δm²₂₁ and m₃² − m₁² = Δm²₃₁ have a unique normal-ordered solution:

| Splittings used | m₁ (meV) | m₂ (meV) | m₃ (meV) | Σ (meV) |
|---|---|---|---|---|
| JUNO 2026 (Δm²₂₁ = 7.388 × 10⁻⁵, Δm²₃₁ = 2.509 × 10⁻³ eV²) | 0.359 | 8.60 | 50.09 | 59.05 |
| NuFIT 6.0, without SK atmospheric data | 0.364 | 8.66 | 50.13 | 59.16 |
| NuFIT 6.0, with SK atmospheric data | 0.362 | 8.66 | 50.34 | 59.36 |

The sum sits just above the normal-ordering minimum (58.8 meV), so the most distinctive number is the lightest mass. The inverted solution (Σ ≈ 0.103 eV) is disfavoured by the cosmological bound, not by the relation.

With TM1 mixing (Xing and Zhou 2007; Albright and Rodejohann 2009) and free Majorana phases, the effective double-beta mass is 1.34–4.02 meV, or 1.82 meV with zero phases and δ = 270°. This is below the KamLAND-Zen limit of 28–122 meV (KamLAND-Zen Collaboration 2024) and below planned reach. The masses themselves do not depend on TM1.

## 4. Tests

- **Ordering:** after its first measurement of the solar parameters (JUNO Collaboration 2025), JUNO favours normal ordering at 2.1σ (JUNO Collaboration 2026) and expects 3σ after about six years.
- **Sum:** DESI, Euclid and next-generation CMB experiments approach the needed sensitivity.
- **Majorana nature:** proof that neutrinos are Dirac particles would remove the double-beta prediction.

Brannen's proposal and the present statement give the same numbers. This test can falsify both but cannot separate them. What this note adds is one principle behind both the charged-lepton and neutrino cases, and an up-to-date evaluation.

## Data and code

`analysis/recalibrate_to_2026_values_2026-11-18.py` and `analysis/derive_45_and_sqrt2_2026-11-18.py` in the TRPG code archive, doi:10.5281/zenodo.24329239. This result is one component of the framework in doi:10.5281/zenodo.24330763.

## Acknowledgements

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.

## References

- Albright, C. H. and Rodejohann, W. (2009). Comparing trimaximal mixing and its variants with deviations from tri-bimaximal mixing. *Eur. Phys. J. C* 62, 599.
- Brannen, C. A. (2006). The lepton masses. Online note, brannenworks.com.
- Esteban, I., Gonzalez-Garcia, M. C., Maltoni, M., Martinez-Soler, I., Pinheiro, J. P. and Schwetz, T. (2024). NuFit-6.0: updated global analysis of three-flavor neutrino oscillations. *JHEP* 12, 216.
- Foot, R. (1994). A note on Koide's lepton mass relation. arXiv:hep-ph/9402242.
- JUNO Collaboration (2025). First measurement of reactor neutrino oscillations at JUNO. arXiv:2511.14593.
- JUNO Collaboration (2026). Measurement of neutrino oscillation parameters at JUNO and indication of the neutrino mass ordering. arXiv:2609.37895.
- KamLAND-Zen Collaboration (2024). Search for Majorana neutrinos with the complete KamLAND-Zen dataset. arXiv:2406.11438.
- Koide, Y. (1983). A fermion–boson composite model of quarks and leptons. *Phys. Lett. B* 120, 161.
- Shannon, C. E. (1948). A mathematical theory of communication. *Bell Syst. Tech. J.* 27, 379.
- Sumino, Y. (2009). Family gauge symmetry as an origin of Koide's mass formula and charged lepton spectrum. *Phys. Lett. B* 671, 477.
- Xing, Z.-z. and Zhou, S. (2007). A partial μ–τ symmetry and its prediction for leptonic CP violation. *Phys. Lett. B* 653, 278.
