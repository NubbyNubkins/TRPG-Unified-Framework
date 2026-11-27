---
title: "Extending a gauge-invariant lattice measure for chiral fermions to vector-like matter, gauge singlets and scalars"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "27 November 2026"
---

## Abstract

Craig (2026) has proposed a smooth, gauge-invariant overlap-fermion measure for the Standard Model's left-handed fermions, 10 ⊕ 5̄ per generation embedded in SU(5), on admissible lattice gauge fields. Many extensions of the Standard Model add three kinds of field to that content:

- vector-like pairs r ⊕ r̄, such as the 5 ⊕ 5̄ of the E₆ 27;
- gauge-singlet fermions, such as right-handed neutrinos with Majorana masses;
- scalars, such as the Higgs doublet or a majoron field.

We show that none of these adds a phase problem.

(i) In four dimensions the left-handed overlap projector in r̄ is the image of the one in r under one fixed antiunitary map, A·conj(·)·A⁻¹ with A = γ₅C⁻¹, for every gauge field. A vector-like pair therefore carries a canonical measure with no phase to fix. Its return phases cancel exactly, so Craig's integrability condition holds for the extended content.

(ii) A global U(1) acting on the fermions is an exact symmetry of the lattice measure in every topological sector if and only if its mixed anomalies with the gauge groups vanish. Thus B − L is exact with three singlets, while B and L each change only by 3 in one unit of weak topology.

(iii) Yukawa, Dirac and Majorana bilinears and scalar fields then give well-defined, gauge-invariant correlation functions at every finite lattice spacing, through a Pfaffian whose sign the measure fixes.

We check (i) numerically to 10⁻¹³ on a 3⁴ SU(3) lattice, with a two-dimensional control in which the map cannot exist and the cancellation indeed fails. We check (ii) on a two-dimensional U(1) lattice. The results are conditional on Craig's theorem, which we do not re-derive. They leave reflection positivity and the continuum limit open, as they are for the Standard Model. As an application, the three-generation 27 with the Standard Model gauge group, one Higgs doublet and a majoron field is defined at finite spacing.

## 1. Setting

On a finite periodic lattice we use the Wilson–Dirac operator D_W with r = 1, the Hermitian H = γ₅(D_W − 1), and the overlap operator D = 1 + γ₅ sign H (Narayanan and Neuberger, 1994; Neuberger, 1998). This operator obeys the Ginsparg–Wilson relation (Ginsparg and Wilson, 1982), so it avoids the doublers that the Nielsen–Ninomiya theorem forces on naive discretizations (Nielsen and Ninomiya, 1981). We follow Lüscher (1999):

- a left-handed field in representation R obeys P̂₋ψ = ψ, where P̂₋ = (1 + sign H)/2 depends on the gauge field U;
- its antifield obeys ψ̄P₊ = ψ̄, where P₊ = (1 + γ₅)/2 is fixed.

The fermion measure is fixed by a basis v_j(U) of the range of P̂₋. A change of basis by a unitary Q multiplies the measure by det Q⁻¹, a phase.

On fields that obey a plaquette bound (admissibility; Lüscher, 1999), Craig (2026) constructs a smooth, gauge-covariant measure current for R = 10 ⊕ 5̄ of SU(5), per generation, restricted to the Standard Model and its quotients. The current satisfies the integrability condition W[γ] = 𝒯[γ] on every closed curve of admissible fields, where 𝒯 is the return phase of the projector family. This fixes the measure up to one constant phase on each topological sector. Earlier steps toward this were made by Lüscher (1999) for abelian theories and by Kadoh and Kikukawa (2008) for the electroweak theory. Craig's paper was prepared with the aid of OpenAI's GPT-6 Astra. Reflection positivity at finite spacing is not established there.

## 2. Vector-like pairs

**Lemma 1.** Let C be a charge-conjugation matrix (Cγ_μC⁻¹ = −γ_μᵀ for Hermitian Euclidean γ_μ) and A = γ₅C⁻¹. For any compact gauge group, any unitary representation r with conjugate r̄, and any gauge field U,

  P̂₋\[r̄\](U) = A · conj(P̂₋\[r\](U)) · A⁻¹.

*Proof.* With Hermitian γ_μ, conj(γ_μ) = γ_μᵀ, and A·γ_μᵀ·A⁻¹ = γ₅C⁻¹γ_μᵀCγ₅ = −γ₅γ_μγ₅ = γ_μ. Since C⁻¹γ₅ᵀC = γ₅, also A·conj(γ₅)·A⁻¹ = γ₅, so A keeps handedness. Complex conjugation maps the links U_μ(x) of r to those of r̄, so A·conj(D_W[r])·A⁻¹ = D_W[r̄]. The same then holds for H, for sign H (a real function of a Hermitian matrix) and for P̂₋. ∎

**Corollary 1.** For the left-handed pair r ⊕ r̄, choose bases v_j for r and v′_j = A·conj(v_j) for r̄. A change v → vQ gives v′ → v′·conj(Q), so the pair's measure changes by |det Q|⁻² = 1. The pair's measure is therefore independent of the basis. A gauge transformation Λ acts on r̄ as conj(Λ) and commutes with the map, so the measure is gauge invariant and smooth on admissible fields. Its return phases obey 𝒯[r̄] = conj(𝒯[r]), and its measure current can be taken to be zero.

For R = R_C ⊕ ⊕ₖ(rₖ ⊕ r̄ₖ) ⊕ singlets, the product of Craig's measure for R_C, the canonical pair measures and the constant singlet measures satisfies W = 𝒯, since both sides multiply. A gauge singlet's projector does not depend on U. When a pair contains a copy of a representation already in R_C (the 5̄ in the 27), different assignments of the identical copies differ by a constant phase on each topological sector, which is the freedom Craig's measure already has.

In two dimensions no handedness-keeping antiunitary map exists. Every map A with A·conj(γ_μ)·A⁻¹ = γ_μ reverses γ₅, because there the conjugate of a left-mover is a left-mover. This is why L(q) ⊕ L(−q) is anomalous in two dimensions. It also gives a control: the cancellation must fail there.

**Numerical check.** We use SU(3) on a 3⁴ lattice with random links near the identity (smallest |H| eigenvalue 0.40), and compute the sign function from the exact spectrum.

- ‖P̂₋[r̄] − A·conj(P̂₋[r])·A⁻¹‖ = 9.6 × 10⁻¹⁴, against ‖P̂₋‖ = 22, and both projectors have rank 486.
- Along a closed loop of gauge fields, using 24 steps of parallel transport by polar decomposition, 𝒯[r] = exp(−0.00531 i) and 𝒯[r]𝒯[r̄] = exp(1.0 × 10⁻¹⁴ i).
- The same procedure in two dimensions (SU(3) on 6²) gives a product phase of 0.26 rad, as expected.

## 3. Global symmetries

**Lemma 2.** Rotate each left-handed field f by exp(i x_f α), with the antifields rotated oppositely. The measure is then multiplied by exp(−iα Σ_f x_f n_f(U)), where n_f is the index of the overlap operator in the representation of f. By the lattice index theorem (Hasenfratz, Laliena and Niedermayer, 1998; Lüscher, 1998), on admissible fields in a given topological sector,

  Σ_f x_f n_f = 2n₃·A(SU(3)²X) + 2n₂·A(SU(2)²X) + n₁·A(Y²X),

where n₃, n₂ and n₁ label the sector and A(·) are the mixed anomaly coefficients, with the usual normalisation T = ½ for a fundamental. The factor 2 on the non-abelian terms appears because one fundamental field has index 1, which is 2T, in one unit of topology; n₁ is normalised so that a field of unit hypercharge has index n₁. Cross terms vanish because the non-abelian generators are traceless. X is therefore an exact symmetry of the lattice measure in every sector if and only if the three mixed anomalies vanish. ∎

For the Standard Model with three gauge-singlet neutrinos, and for the full 27 (whose added states are vector-like), X = B − L has all three coefficients zero, so it is exact at every finite spacing. For X = B or X = L separately, A(SU(2)²X) = 3/2 for three generations. One unit of weak topology (n₂ = 1) then changes B and L by 2 × 3/2 = 3 each, the process of 't Hooft (1976), and no ΔB = 1 process arises from the measure. **Check:** on an 8 × 8 two-dimensional U(1) lattice the index is exactly q·m for charges q = 1, 2, 3 in flux sectors m = 0, 1, 2. Global charges with zero mixed sum leave the measure unchanged in every sector, and a charge with mixed sum 6 shifts it by 6 and 12 units.

## 4. Yukawa, Majorana and scalar terms

In the all-left-handed language, Yukawa couplings, Dirac masses of vector-like pairs and Majorana masses of singlets are field–field and antifield–antifield bilinears. The fermion integral of exp(−S) against the fixed top form defined by the measure is a well-defined number for any polynomial S: it is the Pfaffian of the bilinear form written in the measure's basis. It is gauge invariant whenever every term in S is. The measure fixes the Pfaffian's sign, so there is no square-root ambiguity. The sign is not positive in general, which matters for simulation but not for the definition.

The scalars are added with nearest-neighbour covariant hopping and quartic potentials bounded below, so on a finite lattice their integrals converge. A majoron field σ of B − L charge −2 couples as σ N^c N^c. By Lemma 2 its U(1) is then exact on the lattice, and any explicit breaking is a separate, deliberate term. The θ-term exp(iθ n) is exactly a constant phase per topological sector. This is the freedom that Craig's construction leaves, so the lattice neither fixes nor forbids θ.

## 5. Application: three 27s

Each generation contains:

- 10 ⊕ 5̄ ⊕ 1, the Standard Model fermions and N^c;
- 5 ⊕ 5̄, the vector-like D(3,1,−1/3) and L(1,2,1/2) with their partners;
- one further singlet.

The gauge group is SU(3) × SU(2) × U(1) only. The scalars are one Higgs doublet and a majoron field. All gauge and mixed-gravitational anomaly sums of the 27 vanish, and so do all six B − L anomaly sums. Craig's even-multiplicity condition counts SU(2) doublets on two colours of SU(5): 4 per generation in the 16 and 6 in the 27. Craig's theorem covers 3 × (10 ⊕ 5̄), Lemma 1 covers 3 × (5 ⊕ 5̄), and Section 4 covers the rest. Correlation functions of the model are therefore defined and gauge invariant on every finite admissible lattice. With the vector-like states included, SU(3) and SU(2) remain asymptotically free at one loop: b₃ = 5 and b₂ = 7/6, against 7 and 19/6 without them.

By contrast, gauging B − L, or using the Pati–Salam group (Pati and Salam, 1974), places the 16 in Spin(10). That falls outside Craig's covered families, and the proposed Spin(10) measure (Kikukawa, as cited by Craig) has open smoothness and locality. This model is one of the main-line readings of the framework of Gelfman and Behrends-Gelfman (2026), and its Appendix B states these results for it.

## 6. What remains open

- **Craig's theorem.** Everything here is conditional on it. It is a preprint with a Lean formalization, examined before release by the lattice field theorists M. Golterman, A. Kronfeld and Y. Shamir (as reported by OpenAI, 2026).
- **Reflection positivity.** It is open for the overlap chiral sector. Wilson fermions for the vector-like pairs and nearest-neighbour scalars are reflection positive on their own (Osterwalder and Seiler, 1978), but their couplings to the chiral sector have not been shown to keep this.
- **The continuum limit.** U(1) and the scalar self-couplings are expected to be trivial in four dimensions (Aizenman, 1981; Fröhlich, 1982; Aizenman and Duminil-Copin, 2021), as for the Standard Model.
- **The derivations.** Lemmas 1 and 2 are elementary, but they have not been examined by a lattice specialist.

## Data and code

The scripts `analysis/vectorlike_pair_measure_check_2026-11-27.py`, `analysis/bl_exact_on_lattice_2026-11-27.py`, `analysis/lattice_fit_check_2026-11-27.py` and `analysis/extras_theta_baryon_alternatives_2026-11-27.py`, with their saved outputs, are in the TRPG code archive, doi:10.5281/zenodo.24329239. The framework is described in doi:10.5281/zenodo.24330763.

## Acknowledgements

The construction extended here is due to N. Craig. We also acknowledge M. Golterman, A. Kronfeld and Y. Shamir, whose examination of it is reported by OpenAI (2026). Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.

## References

- Aizenman, M. (1981). Proof of the triviality of φ⁴_d field theory and some mean-field features of Ising models for d > 4. *Phys. Rev. Lett.* 47, 1.
- Aizenman, M. and Duminil-Copin, H. (2021). Marginal triviality of the scaling limits of critical 4D Ising and φ⁴₄ models. *Ann. Math.* 194, 163.
- Craig, N. (2026). A gauge-invariant measure for lattice chiral fermions. arXiv:2610.09093 [hep-lat].
- Fröhlich, J. (1982). On the triviality of λφ⁴_d theories and the approach to the critical point in d ≥ 4 dimensions. *Nucl. Phys. B* 200, 281.
- Gelfman, M. and Behrends-Gelfman, M. (2026). Dynamic Relational Transitioning Geometry. doi:10.5281/zenodo.24330763.
- Ginsparg, P. H. and Wilson, K. G. (1982). A remnant of chiral symmetry on the lattice. *Phys. Rev. D* 25, 2649.
- Hasenfratz, P., Laliena, V. and Niedermayer, F. (1998). The index theorem in QCD with a finite cut-off. *Phys. Lett. B* 427, 125.
- Kadoh, D. and Kikukawa, Y. (2008). A construction of the Glashow–Weinberg–Salam model on the lattice with exact gauge invariance. *JHEP* 05, 095.
- Lüscher, M. (1998). Exact chiral symmetry on the lattice and the Ginsparg–Wilson relation. *Phys. Lett. B* 428, 342.
- Lüscher, M. (1999). Abelian chiral gauge theories on the lattice with exact gauge invariance. *Nucl. Phys. B* 549, 295.
- Narayanan, R. and Neuberger, H. (1994). Chiral determinant as an overlap of two vacua. *Nucl. Phys. B* 412, 574.
- Neuberger, H. (1998). Exactly massless quarks on the lattice. *Phys. Lett. B* 417, 141.
- Nielsen, H. B. and Ninomiya, M. (1981). Absence of neutrinos on a lattice. *Nucl. Phys. B* 185, 20.
- OpenAI (2026). A step toward defining the Standard Model on a lattice. openai.com/index/step-toward-defining-the-standard-model (15 November 2026).
- Osterwalder, K. and Seiler, E. (1978). Gauge field theories on a lattice. *Ann. Phys.* 110, 440.
- Pati, J. C. and Salam, A. (1974). Lepton number as the fourth "color". *Phys. Rev. D* 10, 275.
- 't Hooft, G. (1976). Symmetry breaking through Bell–Jackiw anomalies. *Phys. Rev. Lett.* 37, 8.
