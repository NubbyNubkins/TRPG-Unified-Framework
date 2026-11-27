---
title: "No anomaly-free discrete symmetry splits the vector-like triplets and doublets of three 27s by one power of the breaking field"
author:
  - Matthew Gelfman (ORCID 0009-0000-6402-2258)
  - Melissa Behrends-Gelfman (ORCID 0009-0001-6281-1947)
date: "Version 1.2, 27 November 2026"
---

## Abstract

Models with three generations of vector-like colour triplets D(3,1,−1/3) and weak doublets L(1,2,−1/2), as in the 27 of E₆ (Gürsey, Ramond and Sikivie 1976), often need the doublets far lighter than the triplets. A natural mechanism is a discrete symmetry ℤ~N~ under which a scalar Φ gives D its mass at dimension four, while L can take mass only from Φ²/M, so that M~L~ ≈ M~D~²/M.

We show that no flavour-blind ℤ~N~ does this without a mixed ℤ~N~–SU(3)² anomaly. Requiring SM Yukawa couplings with one Higgs doublet makes the quarks' contribution to that anomaly vanish, which leaves 3·q~Φ~ from the three D states. Cancelling it forces Φ³ to be neutral. Then Φ² carries the charge of Φ\*, so L gets its mass at dimension four as well, and the hierarchy is lost. The argument holds for every N and every lepton, Higgs or quark charge assignment, and does not use the gravitational condition.

A brute-force scan of every flavour-blind assignment for N ≤ 18 confirms this. With the weaker "half-strength" form of the gauge conditions used in some of the literature, a single solution survives (ℤ₆, with the leptons at charge 3), and only if the gravitational condition is also taken at half strength. We correct a miscounted gravitational anomaly that circulated for this ℤ₆ case. Ways around the result need extra fields, a Green–Schwarz mechanism, or explicit breaking. The last is incompatible with using the symmetry for protection (Krauss and Wilczek 1989).

## 1. Setup

Take three generations with flavour-blind ℤ~N~ charges: q~Q~, q~u~, q~d~, q~l~, q~e~, q~N~ for the SM fermions (right-handed fields counted by their own charges), h for the Higgs doublet and q~Φ~ for a singlet scalar Φ. Yukawa invariance gives q~u~ = q~Q~ + h, q~d~ = q~Q~ − h, q~e~ = q~l~ − h and q~N~ = q~l~ + h, with the sign conventions absorbed into h.

Let Δ~D~ = q(D~R~) − q(D~L~) and Δ~L~ = q(L~R~) − q(L~L~). We require four things:

- **(a) D mass at dimension four:** Δ~D~ ≡ ∓q~Φ~ (mod N) and Δ~D~ ≢ 0.
- **(b) L mass only at dimension five:** Δ~L~ ≢ 0 and Δ~L~ ≢ ±q~Φ~, but Δ~L~ ≡ ∓2q~Φ~.
- **(c) No baryon-number-violating operators with a lepton:** QQQL, u^c^ u^c^ d^c^ e^c^ and QQu^c^ e^c^ are all forbidden.
- **(d) No mixed anomaly with SU(3)², SU(2)² or gravity.**

## 2. Anomaly conditions

For left-handed Weyl fields, an instanton of SU(M) has 2T(r) zero modes per field in representation r, so the 't Hooft vertex ('t Hooft 1976) carries ℤ~N~ charge A~M~ = Σ q~i~ · 2T(r~i~). Anomaly freedom requires A~M~ ≡ 0 (mod N) (Ibáñez and Ross 1991; Banks and Dine 1992). For gravity, the corresponding condition is A~grav~ = Σ q~i~ dim(r~i~) ≡ 0 modulo N, or modulo N/2 for even N.

The SM quarks contribute 3(2q~Q~ − q~u~ − q~d~) = 0 to A₃ by Yukawa invariance. The D states contribute 3(q(D~L~) − q(D~R~)) = −3Δ~D~ = ±3q~Φ~. Therefore

A₃ = ±3 q~Φ~.

## 3. Theorem

*No flavour-blind ℤ~N~ with three generations satisfies (a), (b) and the strict condition A₃ ≡ 0 (mod N).*

*Proof.* A₃ ≡ 0 requires 3q~Φ~ ≡ 0 (mod N), so 2q~Φ~ ≡ −q~Φ~. Condition (b) needs Δ~L~ ≡ ∓2q~Φ~ ≡ ±q~Φ~, which contradicts Δ~L~ ≢ ±q~Φ~. ∎

The proof uses neither (c) nor the SU(2)² and gravitational conditions. With N~g~ generations it generalises: anomaly freedom needs N~g~ q~Φ~ ≡ 0, and the obstruction is that 2q~Φ~ must not equal ±q~Φ~ or 0. For N~g~ = 3 the only possibility is the one excluded above.

## 4. Scan and the half-strength conventions

An exhaustive scan over all (q~Q~, q~l~, h, Δ~D~, Δ~L~) for N = 2 to 18 under (a) to (c) gives:

- **Strict gauge conditions (mod N):** no solution, with the gravitational condition at either strength or omitted.
- **Gauge conditions at half strength (mod N/2):** only ℤ₆ survives, and only if the gravitational condition is also taken mod N/2 or omitted. With the gravitational condition mod N, nothing survives.

For the ℤ₆ case (q~Q~ = 0, q~l~ = 3, q~Φ~ = 1, Δ~D~ = −1, Δ~L~ = −2), the sums are A₃ = 3, A₂ = 15 and A~grav~ = 21. The gravitational sum has been quoted as 30, which treats the lepton doublet and e^c^ as cancelling. Counted field by field, the lepton doublet, e^c^ and N^c^ give 2·3 − 3 − 3 = 0 per generation, so the Standard Model part is 0 and the vector-like part is 21. The Standard Model fields alone then pass every half-strength condition, so the vector-like states are not needed for anomaly freedom in that convention either.

Charge 3 for the leptons is lepton parity inside ℤ₆. That is what forbids the baryon-number-violating operators with a lepton; the "81 flavour combinations" sometimes quoted are one condition.

## 5. Ways around

- **Extra fields.** Extra fields charged under SU(3) can cancel A₃, but they add matter.
- **Green–Schwarz.** An axion-like field can shift to cancel the anomaly (Green and Schwarz 1984).
- **Explicit breaking.** The symmetry can be broken explicitly, for example by Planck-suppressed operators. But a symmetry broken by gravity gives no protection, and the "dimension-five gravitational bias" sometimes used to remove ℤ~N~ domain walls cannot be combined with anomaly freedom.

Without one of these, the relation M~L~ ≈ M~D~²/M must be taken as an assumption about couplings, not as a symmetry statement.

## Data and code

`analysis/farside_charge_rule_2026-11-25.py` and `analysis/refinement_checks_2026-11-25.py` in the TRPG code archive, doi:10.5281/zenodo.24329239. This result bears on one assumption of the framework in doi:10.5281/zenodo.24330763.

## Acknowledgements

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it.


This paper is licensed under CC BY 4.0; its scripts and data are under the MIT License.

## References

- Banks, T. and Dine, M. (1992). Note on discrete gauge anomalies. *Phys. Rev. D* 45, 1424.
- Green, M. B. and Schwarz, J. H. (1984). Anomaly cancellation in supersymmetric D = 10 gauge theory and superstring theory. *Phys. Lett. B* 149, 117.
- Gürsey, F., Ramond, P. and Sikivie, P. (1976). A universal gauge theory model based on E6. *Phys. Lett. B* 60, 177.
- Ibáñez, L. E. and Ross, G. G. (1991). Discrete gauge symmetry anomalies. *Phys. Lett. B* 260, 291.
- Krauss, L. M. and Wilczek, F. (1989). Discrete gauge symmetry in continuum theories. *Phys. Rev. Lett.* 62, 1221.
- 't Hooft, G. (1976). Symmetry breaking through Bell–Jackiw anomalies. *Phys. Rev. Lett.* 37, 8.
