---
title: "An exact right-angled unitarity triangle with λ = 1/√20 and A = √(2/3): a testable parametrisation of quark mixing"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "25 November 2026"
---

## Abstract

We take the right angle of the unitarity triangle at α as exact, as quark-mixing sum rules allow (Antusch, King, Malinsky and Spinrath 2010), and choose simple values for the other three quark-mixing parameters: λ = 1/√20, A = √(2/3), β = 22.5° and γ = 67.5°. The result is an exactly unitary CKM matrix (Kobayashi and Maskawa 1973) with Jarlskog invariant J = 2.94 × 10⁻⁵ (Jarlskog 1985) and sin 2β = 1/√2. All seven directly measured quantities lie within 1.6σ, with pulls from −1.2σ (|V~ub~|) to +1.6σ (α), giving χ² = 6.0 for seven quantities. The four values were chosen, not derived, so the χ² is not evidence. We present the form as a parametrisation with sharp, near-term tests. The matrix is exactly unitary and gives |V~ud~| = 0.97467. That is 3.1σ above the superallowed-decay value 0.97367 ± 0.00032, so the form bets that the first-row shortfall (the "Cabibbo anomaly") goes away. Belle II and LHCb measurements of sin 2β, α and γ at the 1° level will decide it.

## 1. Construction

In the exact form of the Wolfenstein parametrisation (Wolfenstein 1983; Charles et al. 2005), the inputs are λ = s₁₂ = 1/√20 = 0.22361 and s₂₃ = Aλ² with A = √(2/3). They are followed by ρ̄ + iη̄ = R~u~ e^iγ^, where R~u~ = sin β because the angle at α is 90°. The relation

s₁₃ e^iδ^ = (ρ̄ + iη̄) s₁₂ s₂₃ c₂₃ / [c₁₂ (1 − (ρ̄ + iη̄) s₂₃²)]

then gives an exactly unitary matrix in the standard form. Numerically, the unitarity error is 2 × 10⁻¹⁶.

## 2. Comparison with measurement

The measured values are PDG 2026 for |V~us~| = 0.22431, α = 84.1 (+3.7 −3.0)°, γ = 66.4 (+2.7 −2.8)° and the first row. For the other quantities we use PDG 2024 (Navas et al. 2024): |V~cb~| = 0.0411 ± 0.0012, |V~ub~| = (3.82 ± 0.20) × 10⁻³, sin 2β = 0.709 ± 0.011 and |V~td~/V~ts~| = 0.207 ± 0.0032 (0.001 statistical and 0.003 lattice, combined). We use direct measurements only, not global-fit values, which already assume unitarity.

| Quantity | Prediction | Measured | Pull |
|---|---|---|---|
| \|V~us~\| | 0.22361 | 0.22431 ± 0.00085 | −0.83 |
| \|V~cb~\| | 0.04082 | 0.0411 ± 0.0012 | −0.23 |
| \|V~ub~\| | 0.00358 | 0.00382 ± 0.00020 | −1.19 |
| sin 2β | 0.70711 | 0.709 ± 0.011 | −0.17 |
| α (°) | 90.0 | 84.1 (+3.7 −3.0) | +1.59 (upper error used) |
| γ (°) | 67.5 | 66.4 (+2.7 −2.8) | +0.39 (2.8 used) |
| \|V~td~/V~ts~\| | 0.21035 | 0.207 ± 0.0032 | +1.05 |

χ² = 5.98. If nothing had been chosen this would correspond to p = 0.54, but it has no such meaning here.

## 3. What was chosen, and what fails

The values 1/√20, √(2/3), 22.5° and 67.5° = 90° − 22.5° were selected because they are simple. A 3 × 3 down-quark mass matrix with fixed weights built to produce this triangle gives χ² = 11.5 at μ = 162.5 GeV (|V~us~| +2.2σ) and 26 at M~Z~, so the triangle is not yet tied to a mass matrix. The form is offered as a target for a mass-matrix origin, not as one.

## 4. Tests

- sin 2β firmly away from 1/√2 = 0.7071.
- α firmly away from 90°. Today the direct value is 1.6σ below.
- γ firmly outside about 64–68°.
- A confirmed first-row shortfall, which would also rule out any model that keeps the CKM matrix exactly unitary at low energy.

## Data and code

`analysis/quark_triangle_and_press_2026-11-25.py` in the TRPG code archive, doi:10.5281/zenodo.24329239. This result is one component of the framework in doi:10.5281/zenodo.24330763.

## Acknowledgements

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.

## References

- Antusch, S., King, S. F., Malinsky, M. and Spinrath, M. (2010). Quark mixing sum rules and the right unitarity triangle. *Phys. Rev. D* 81, 033008.
- Charles, J. et al. (CKMfitter Group) (2005). CP violation and the CKM matrix: assessing the impact of the asymmetric B factories. *Eur. Phys. J. C* 41, 1.
- Jarlskog, C. (1985). Commutator of the quark mass matrices in the standard electroweak model and a measure of maximal CP nonconservation. *Phys. Rev. Lett.* 55, 1039.
- Kobayashi, M. and Maskawa, T. (1973). CP-violation in the renormalizable theory of weak interaction. *Prog. Theor. Phys.* 49, 652.
- Navas, S. et al. (Particle Data Group) (2024). Review of particle physics. *Phys. Rev. D* 110, 030001; and the 2026 edition (summary tables and CKM review).
- Wolfenstein, L. (1983). Parametrization of the Kobayashi–Maskawa matrix. *Phys. Rev. Lett.* 51, 1945.
