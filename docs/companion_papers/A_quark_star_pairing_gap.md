---
title: "The colour-superconducting gap required in self-bound quark stars when the bag constant is fixed by the nucleon and Δ(1232)"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "25 November 2026"
---

## Abstract

In most studies of strange quark stars, the bag constant B is adjusted to the observed stars. Here B is instead fixed beforehand by fitting the nucleon and Δ(1232) masses in a cloudy-bag model. That model has a centre-of-mass correction and no free zero-point term, and the fit gives B = 60.1 MeV fm⁻³ (57.1–63.7 across the measured pion–nucleon coupling). With this B and a strong coupling run from α~s~(M~Z~), unpaired strange-quark matter gives a heaviest star of only 1.54–1.80 M☉, and strange matter is not the ground state. Pairing in the colour–flavour-locked (CFL) phase repairs both. For m~s~ = 100 MeV, the stars meet the observed heaviest masses (2.0–2.15 M☉) and the GW170817 tidal limit only if the gap is Δ ≈ 82–97 MeV. Equivalently, the stars fix 4Δ² − m~s~² ≈ 1.7–2.8 × 10⁴ MeV². An independent solver reproduces the band at the central B: the heaviest star is 2.00–2.13 M☉ and Λ₁.₄ = 459–575. Within the band, strange matter is self-bound (E/A = 832–864 MeV) and paired two-flavour matter stays above 930 MeV, so ordinary nuclei are stable. Larger gaps (90–125 MeV at fixed a₄) push Λ₁.₄ to 652–1124, above the GW170817 limit. The observed stars set the band, so they are not evidence for it. The prediction is the gap itself, to be tested by an independent QCD determination.

## 1. Bag constant from hadron spectroscopy

The MIT bag energy (Chodos et al. 1974; DeGrand et al. 1975) is E = (4π/3)BR³ + (3x₀ − Z₀)/R plus a colour-magnetic term, with x₀ = 2.0428 for massless quarks. We keep three established ingredients and set Z₀ = 0:

- the centre-of-mass correction M = √(E² − ⟨P²⟩) (Donoghue and Johnson 1980);
- the pion cloud of the cloudy bag model (Théberge, Thomas and Miller 1980), with the measured pion–nucleon coupling f²/4π = 0.0778;
- the colour-magnetic interaction.

The two unknowns, B and α~c~, are fitted to M~N~ = 938.9 MeV and M~Δ~ = 1232.0 MeV, and each radius is set by its own pressure balance. As a check, the 1975 setup (Z₀ = 1.84, no correction, no cloud) gives B^1/4^ = 146 MeV, close to the original 145 MeV.

The full model gives B = 60.1 MeV fm⁻³ (B^1/4^ = 147 MeV), α~c~ = 1.55 and a proton bag radius of 0.995 fm. Moving the pion–nucleon coupling by ±5% gives B = 57.1–63.7 MeV fm⁻³. No compact-star input enters.

## 2. Quark matter and stellar structure

The pressure of three-flavour quark matter is (Alford, Braby, Paris and Reddy 2005)

P(μ) = (3/4π²)[a₄ μ⁴ − (m~s~² − 4Δ²) μ²] − B,

where a₄ = 1 − 2α~s~/π. α~s~ is run at two loops from α~s~(M~Z~) = 0.1180 to the scale 3μ. The quark density is n = dP/dμ and the energy density is ε = μn − P, so the equation of state is thermodynamically consistent. We solve the TOV equations, and the tidal Love number from Hinderer's equation (Hinderer 2008) with the surface-density correction for self-bound stars (Postnikov, Prakash and Lattimer 2010).

## 3. Results

Without pairing, the heaviest star is 1.54–1.80 M☉ across the B band, which fails PSR J0740+6620 (Cromartie et al. 2020; Fonseca et al. 2021). Unpaired strange matter is also unbound (E/A = 940–962 MeV).

With CFL pairing at m~s~ = 100 MeV and B = 60.1 MeV fm⁻³:

| Δ (MeV) | M~max~ (M☉) | Λ₁.₄ | E/A, CFL (MeV) |
|---|---|---|---|
| 82 | 2.004 | 459 | 864.1 |
| 97 | 2.125 | 575 | 832.1 |

These values come from an independent solver (script `cfl_gap_and_tides_2026-11-25.py`, run through `refinement_checks_2026-11-25.py`). They agree with the original scan (`pairing_gap_scan_2026-11-14.py`) to 0.01 M☉ and about 3% in Λ.

Across B = 57.1–63.7 MeV fm⁻³ the band moves to about 75–108 MeV. Because the stars fix the combination 4Δ² − m~s~² ≈ 1.7–2.8 × 10⁴ MeV², a running strange mass of 120–140 MeV corresponds to Δ ≈ 88–109 MeV.

**Ground state.** At B = 60.1, two-flavour matter with 2SC pairing at the same gap has E/A = 954–969 MeV, above the 930 MeV of iron, so ordinary nuclei are stable (`ns_ground_state_2026-11-18.py`). At the low edge of B the margin over 930 MeV shrinks to 5–7 MeV.

**Larger gaps fail.** At fixed a₄ (0.70–1.00) and m~s~ = 95 MeV, gaps of 90–125 MeV give heaviest stars of 2.20–2.61 M☉ but Λ₁.₄ = 652–1124. That is above the GW170817 limit of about 580 (Abbott et al. 2018) unless the 1.4 M☉ stars sit on a separate hadronic branch (Drago, Lavagno and Pagliara 2014).

## 4. What is and is not tested

The gap band was chosen so that the stars pass, so the star masses and Λ₁.₄ are not evidence for it. The testable statement is the band itself. A first-principles determination at μ ≈ 400–500 MeV that puts 4Δ² − m~s~² outside 1.7–2.8 × 10⁴ MeV² would rule out this branch. So would a firmly measured compact star well above 2.15 M☉: PSR J0952−0607 at 2.35 ± 0.17 M☉ (Romani et al. 2022) and the 2.6 M☉ object in GW190814 are the cases to watch. Functional-renormalisation-group, Dyson–Schwinger and effective-theory estimates currently span 50–150 MeV.

## Data and code

All scripts and saved outputs are in the TRPG code archive, doi:10.5281/zenodo.24329239 (all versions), in `code/` and `analysis/`. This result is one component of the framework described in Gelfman and Behrends-Gelfman, doi:10.5281/zenodo.24330763.

## Acknowledgements

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.

## References

- Abbott, B. P. et al. (2018). GW170817: measurements of neutron star radii and equation of state. *Phys. Rev. Lett.* 121, 161101.
- Alford, M., Braby, M., Paris, M. and Reddy, S. (2005). Hybrid stars that masquerade as neutron stars. *Astrophys. J.* 629, 969.
- Chodos, A., Jaffe, R. L., Johnson, K., Thorn, C. B. and Weisskopf, V. F. (1974). New extended model of hadrons. *Phys. Rev. D* 9, 3471.
- Cromartie, H. T. et al. (2020). Relativistic Shapiro delay measurements of an extremely massive millisecond pulsar. *Nature Astron.* 4, 72.
- DeGrand, T., Jaffe, R. L., Johnson, K. and Kiskis, J. (1975). Masses and other parameters of the light hadrons. *Phys. Rev. D* 12, 2060.
- Donoghue, J. F. and Johnson, K. (1980). Pion and nucleon form factors in the bag model. *Phys. Rev. D* 21, 1975.
- Drago, A., Lavagno, A. and Pagliara, G. (2014). Can very compact and very massive neutron stars both exist? *Phys. Rev. D* 89, 043014.
- Fonseca, E. et al. (2021). Refined mass and geometric measurements of the high-mass PSR J0740+6620. *Astrophys. J. Lett.* 915, L12.
- Hinderer, T. (2008). Tidal Love numbers of neutron stars. *Astrophys. J.* 677, 1216.
- Postnikov, S., Prakash, M. and Lattimer, J. M. (2010). Tidal Love numbers of neutron and self-bound quark stars. *Phys. Rev. D* 82, 024016.
- Romani, R. W. et al. (2022). PSR J0952−0607: the fastest and heaviest known Galactic neutron star. *Astrophys. J. Lett.* 934, L17.
- Théberge, S., Thomas, A. W. and Miller, G. A. (1980). Pionic corrections to the MIT bag model. *Phys. Rev. D* 22, 2838.
