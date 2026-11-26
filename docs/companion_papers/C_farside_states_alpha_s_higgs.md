---
title: "Vector-like states of the 27 at two thresholds: a dated prediction α~s~(M~Z~) = 0.1178 and the Higgs mass from a vanishing quartic at M~Pl~/√2"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "25 November 2026"
---

## Abstract

Each 27 of E₆ contains, besides one Standard Model generation and a neutral singlet, a vector-like colour triplet D(3,1,−1/3) and a vector-like weak doublet L(1,2,−1/2) (Gürsey, Ramond and Sikivie 1976). We place three generations of each at two thresholds. The D states sit at M~D~ = 4.6 × 10¹³ GeV, and the doublets at M~L~ = M~D~²/E~g~ = 2.4 × 10⁸ GeV, where E~g~ = M~Pl~/√2. With the measured electroweak inputs, the three gauge couplings then cross at M~D~ to 0.02 in 1/α at two loops. Closing the crossing at three loops requires α~s~(M~Z~) = 0.1178. That is 0.2σ from the world average (0.1180 ± 0.0009) and 1.7σ below the most precise lattice determination (0.11873 ± 0.00056; Dalla Brida et al. 2025), so lattice results at about 0.3% will decide it.

The same states change the running of the Higgs quartic. With the boundary condition λ(E~g~) = 0 of Shaposhnikov and Wetterich (2010), the Higgs mass is 127.2 ± 1.5 GeV, against 125.13 ± 0.11 GeV (+1.4σ). The Standard Model alone gives 128.8 ± 1.5 GeV (+2.6σ). Read the other way, the rule needs a top pole mass of 171.6 GeV, 3.7σ below today's 172.60 ± 0.27 GeV.

Three choices were made with the data in view: E~g~ (the best of five equally natural readings), the M~L~ rule, and equality of the couplings at M~D~. We state each one, together with what would falsify the construction.

## 1. Field content and thresholds

Each generation adds one Dirac colour triplet D with hypercharge −1/3 and one Dirac weak doublet L with hypercharge −1/2. No coloured scalars are added, so baryon number is conserved at the renormalisable level, consistent with the proton-lifetime limit above 2.4 × 10³⁴ years (Takenaka et al. 2020). The D states mix only with d~R~. Their mixing with left-handed quarks is suppressed by (v/M~D~)² ≈ 10⁻²³, so first-row CKM unitarity holds to that level.

In GUT normalisation, each Dirac D contributes Δb = (4/15, 0, 2/3) at one loop, and each Dirac L contributes Δb = (2/5, 2/3, 0). The full two-loop matrices are used (Machacek and Vaughn 1984; Luo, Wang and Xiao 2003).

## 2. The crossing and α~s~

With sin²θ~W~ = 0.23122, 1/α(M~Z~) = 127.951 and α~s~(M~Z~) = 0.1178, the inverse couplings (GUT-normalised) run as follows:

| μ (GeV) | 1/α₁ | 1/α₂ | 1/α₃ |
|---|---|---|---|
| 2.4 × 10⁸ | 49.31 | 36.92 | 25.33 |
| 4.6 × 10¹³ | 38.97 | 38.99 | 38.99 |
| 8.63 × 10¹⁸ | 27.10 | 41.08 | 48.59 |

The couplings cross at M~D~ and separate again above it. No placement of the two thresholds makes them equal at E~g~: the smallest spread found is 16.8 in 1/α. The construction therefore takes equality at M~D~, where each 27 becomes complete, as a boundary condition, because there is no unified group to enforce it.

The prediction moves with the choice of E~g~. With the unreduced Planck energy it is α~s~ = 0.1196, which is not the value claimed here. Because E~g~ was chosen as the best of five readings, the prediction is tested only by measurements made after it was fixed (17 November 2026).

## 3. The Higgs mass

The two-loop running starts from the NNLO matching of Buttazzo et al. (2013), and the Standard Model's three-loop shift is added as an estimate. The uncertainty includes the top pole mass (±0.27 GeV plus 0.5 GeV for the pole/Monte-Carlo ambiguity), α~s~ and missing higher orders. The condition λ(E~g~) = 0 then gives:

- Standard Model only: m~H~ = 128.8 ± 1.5 GeV.
- With the D and L states: m~H~ = 127.2 ± 1.5 GeV. Without the 0.5 GeV allowance the tension is 2.1σ.

The pass line, fixed before running, was 2σ. The result passes only with the allowance.

## 4. What the construction does not do

The rule M~L~ = M~D~²/E~g~ was found after the needed value (2.2 × 10⁸ GeV) was known. A companion note shows that no flavour-blind discrete charge, of any order, can enforce it without a gauge anomaly. The rule is therefore an assumption: the doublets couple to the field that gives D its mass only through gravitational-strength terms.

If the universe reheated above M~L~, the doublets must decay through a small mixing with the lepton doublets, of order 10⁻¹⁶ or larger depending on the reheating temperature. That adds one coupling.

## 5. Falsifiers

The construction fails if any of the following is observed:

- α~s~(M~Z~) more than 3σ from 0.1178 at about 0.3% precision;
- a top mass that puts the Higgs condition more than 3σ off;
- any proton decay (today's limit: Takenaka et al. 2020);
- a confirmed shortfall in first-row CKM unitarity (today 0.9984 ± 0.0007, a 2.3σ shortfall);
- coloured scalars or extra vector-like matter below M~L~.

## Data and code

`analysis/leak_establish_three_loop_2026-11-16.py`, `analysis/crossing_vs_grain_2026-11-18.py`, `analysis/higgs_at_the_grain_2026-11-17.py` and `analysis/higgs_three_loop_check_2026-11-18.py` in the TRPG code archive, doi:10.5281/zenodo.24329239. This result is one component of the framework in doi:10.5281/zenodo.24330763.

## Acknowledgements

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.

## References

- Buttazzo, D., Degrassi, G., Giardino, P. P., Giudice, G. F., Sala, F., Salvio, A. and Strumia, A. (2013). Investigating the near-criticality of the Higgs boson. *JHEP* 12, 089.
- Dalla Brida, M., Höllwieser, R., Knechtli, F., Korzec, T., Ramos, A., Sint, S. and Sommer, R. (ALPHA) (2025). The strength of the interaction between quarks and gluons. arXiv:2501.06633.
- Gürsey, F., Ramond, P. and Sikivie, P. (1976). A universal gauge theory model based on E6. *Phys. Lett. B* 60, 177.
- Luo, M., Wang, H. and Xiao, Y. (2003). Two-loop renormalization group equations in general gauge field theories. *Phys. Rev. D* 67, 065019.
- Machacek, M. E. and Vaughn, M. T. (1984). Two-loop renormalization group equations in a general quantum field theory. II. Yukawa couplings. *Nucl. Phys. B* 236, 221.
- Shaposhnikov, M. and Wetterich, C. (2010). Asymptotic safety of gravity and the Higgs boson mass. *Phys. Lett. B* 683, 196.
- Takenaka, A. et al. (Super-Kamiokande Collaboration) (2020). Search for proton decay via p → e⁺π⁰ and p → μ⁺π⁰ with an enlarged fiducial volume in Super-Kamiokande I–IV. *Phys. Rev. D* 102, 112011.
