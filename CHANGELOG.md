# Changelog

## 1.4.2 (27 November 2026)

New paper Appendix B and companion paper F: the particle side on a lattice. No result, rule, prediction, score or count changed.

- Craig (2026, arXiv:2610.09093) has proposed a smooth, gauge-invariant overlap-fermion measure for the Standard Model's left-handed fermions (10 ⊕ 5̄ of SU(5) per generation). Two lemmas extend it to the rest of the main line.
  - **Lemma A.** In four dimensions the left-handed projector in r̄ is the image of the one in r under the fixed antiunitary map A·conj(·)·A⁻¹ with A = γ₅C⁻¹. A vector-like pair (the far-side 5 ⊕ 5̄) therefore has a canonical, phase-free, gauge-invariant measure. Checked to 10⁻¹³ on a 3⁴ SU(3) lattice; in a two-dimensional control the cancellation fails, as it must. `analysis/vectorlike_pair_measure_check_2026-11-27.py`.
  - **Lemma B.** A global charge is exact on the lattice in every topological sector if and only if its mixed anomalies vanish. B − L, the majoron's symmetry, is exact on the 27. B and L change only by 3 in one unit of weak topology, so the lattice gives no proton decay. `analysis/bl_exact_on_lattice_2026-11-27.py`, `analysis/extras_theta_baryon_alternatives_2026-11-27.py`.
- With these lemmas the Higgs, the mass links, the Majorana mass and the majoron field are defined at every finite spacing. The fermion integral is a Pfaffian whose sign the measure fixes.
- Anomaly sums, Craig's even-count condition (4 doublets per generation in the 16, 6 in the 27) and the one-loop coefficients (b₃ = 5, b₂ = 7/6) are in `analysis/lattice_fit_check_2026-11-27.py`.
- The constant phase per topological sector that Craig's measure leaves free is exactly the strong-interaction angle, which stays a measured input.
- Only the main line is covered. The gauged B − L and 4 × 4 alternatives place the 16 in Spin(10), outside Craig's covered families.
- **Independent recheck of Craig's proof** (paper B.7, F Section 6): every finite, computable step was rechecked from his definitions, 18 of 18 checks pass and no error was found. That covers the colour-trace cancellation and where it stops, his two content tests, the free-field facts, his Lemma 5.1 and the algebra of his global step. His analytic estimates, the cited gap bound, and his Lean and two-dimensional supplements were not checked. `analysis/craig_proof_recheck_2026-11-27.py`.
- **Correction found by the recheck:** Craig's Proposition 11.1 and Table 3 already cover the full 27 (zero cubic sum, 6 doublets, and pairs Q ⊕ Q* are listed). The earlier text said the far-side pair lay outside his theorem. Lemma A is now stated as a second, explicit route for the pair; Lemma B and the corollary remain additions.
- **Still open:** Craig's theorem itself (a preprint), reflection positivity and the zero-spacing limit, as for the Standard Model. No number is derived.
- The PDFs of the paper and of companion paper F are now produced with XeLaTeX (the LibreOffice export dropped the bar in 5̄).
- Corrections after an outside check (27 November): the index-theorem sum in Lemma B now carries the factor 2 on the non-abelian terms (2n₃·A(SU(3)²X) + 2n₂·A(SU(2)²X) + n₁·A(Y²X), with T = ½ for a fundamental), matching the script; no result changes. The PDFs are rebuilt with XeLaTeX so barred symbols such as 5̄ print correctly, and the brackets in P̂₋[r̄] are restored in the Word files.

## 1.4.1 (25 November 2026)

The paper's Acknowledgements now state fully how AI assistants were used, as arXiv asks: "Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to carry out and check calculations, and to search the literature; every result was checked by rerunning the published code. The authors reviewed all content and take full responsibility for it." The README says the same. No result, rule, prediction or count changed.

## 1.4.0 (25 November 2026)

New tests, each rechecked with independent code. No particle-physics value of 1.3.3 changed. On the main line the scorecard (7 met, 2 partly met) and the input count (32 today, 28 if every pending measurement agrees, plus six stated choices) are unchanged. The framework is restated, and a start phase is kept as the alternative, as below.

- New lead: the exact right-angled quark triangle, λ = 1/√20, A = √(2/3), β = 22.5°, γ = 67.5°. All seven directly measured mixing quantities lie within 1.6σ (χ² = 6.0 for 7; the four numbers are chosen, so this is not evidence). It shares the first-row prediction: |V_ud| = 0.97467, 3.1σ above today's value. `analysis/quark_triangle_and_press_2026-11-25.py`.
- Tested and failed: a six-step (ℤ₆) charge that would turn the far-side weight rule into a symmetry. It gives 2.45 × 10⁸ GeV but is anomalous, and a search of every flavour-blind ℤ_N up to N = 18 finds no charge that works, so the rule stays fitted choice 2. `analysis/farside_charge_rule_2026-11-25.py`, `analysis/refinement_checks_2026-11-25.py`.
- Tested and failed (paper Appendix A.7, ledger section 6a): three generations from that charge; removing its walls by gravity; a very light majoron as all the dark matter without inflation; a pairing gap of 90–125 MeV (tidal deformability 652–1124 at 1.4 solar masses, above GW170817; `analysis/cfl_gap_and_tides_2026-11-25.py`); the matter excess with the Dirac link aligned to the neutrinos; a dark-energy formula from the containment pressure; a 3 × 3 down-quark matrix with set weights; simple fractions for the tilt of the seeds of galaxies; the Higgs mass as 246/2 GeV (self-coupling 1/8 is reached near 182 GeV, not at v; `analysis/refinement_checks_2026-11-25.py`); three generations in one octonionic matrix; integer loop rungs for the masses (matched by random masses 81% of the time).
- Rechecked with an independent solver: the pairing-gap band of Section 9 (82–97 MeV at m_s = 100 MeV) gives a heaviest star of 2.00–2.13 solar masses and Λ₁.₄ = 459–575 at B = 60.1 MeV/fm³, as stated.
- New paper Appendix A.8: where each open problem stands, with what holds, what failed and what is missing.
- Paper (`docs/TRPG_Paper.*`) is now version 1.4.0, with the status note dated 25 November 2026 and the input-count range stated after each "about 28", as on the website.
- **Five rules and one bet** (paper Section 2, README, ledger). Rule 1: one matrix of relationships read in three roles. Held states (matter) sit on the diagonal, the links between them are the twelve carriers, the link across one state is the Higgs, and space and time are the geometry of the relationships. Rule 2: balance of pull and resistance. Rule 3: the 45° rule for leptons. Rule 4: nothing extra. Rule 5: a smooth start. The bet is unchanged. Rules 1–4 restate the earlier three ideas and one rule; the earlier statement is kept in the paper.
- **Start phase, kept as the alternative start** (paper Sections 5 and 9; the main line follows CMB + DESI, which disfavour its tilt, and keeps the tilt as a measured input; `analysis/start_phase_starobinsky_2026-11-19.py`). The stiffness of space (an R² term; Starobinsky, 1980) stretches the universe for 54–55 e-folds.
  - It explains flatness and the uniform temperature of the sky.
  - It predicts r = 0.0040 ± 0.0003 and n_s = 0.9630–0.9636.
  - The tilt agrees with Planck (−0.1σ) and is 2.1σ low against Planck + ACT + SPT. It is 3.2–4.2σ low once DESI BAO is combined (McDonough and Ferreira, 2025), so it is under strain.
- Knock-on effects if the alternative start holds:
  - The coloured far-side states are never produced.
  - The far-side weak pairs are produced, so their decay coupling (above about 9 × 10⁻¹²) is now required.
  - The majoron's starting angle is averaged.
  - The majoron lead's r below 10⁻¹¹ moves to the cold-start alternative.
  - A seventh fitted choice is added: the Higgs coupling to curvature, ξ ≈ 0.06–0.15.
  - The matter excess gets a borderline non-thermal source and stays an input.
- **Count with the alternative start:** 33 today (one more than 1.3.3, for the decay coupling) and 28 if every pending measurement agrees, the tilt now included (`analysis/input_ledger_2026-11-19.py`; the 18 November script is kept).
- New test scripts, each with its saved output:
  - `fewer_particles_tests_2026-11-19` (counts from data, anomaly cancellation of one generation, force-strength meeting under removals, Koide, top condensation, preons);
  - `forces_resistances_matrix_2026-11-19` (lattice-style invariance of the matrix layout, pull and resistance, dark energy as tension, balance radii);
  - `four_left_out_parts_2026-11-19` (3 + 1 from a 2 × 2 complex matrix, quark-mixing relations at a consistent scale, uncertainty from waves, start conditions).
- New figure: `figures/relationship_matrix_diagonal_links.png`.
- Not changed: every particle-physics prediction, the scorecard (7 met, 2 partly met) and the derivation-tree image (the ledger's diagram gains rule 5).

## 1.3.3 (18 November 2026)

An audit of the whole framework against the newest measurements, and of every test for whether it was set after an earlier attempt failed. No rule or reading changed.

- New: `analysis/recalibrate_to_2026_values_2026-11-18.py` (with its output) recomputes every headline number from the Particle Data Group's 2026 values and JUNO 2026.
  - Higgs mass from the grain: 127.3 ± 1.5 GeV against 125.13 ± 0.11 (1.4σ; 2.1σ without the 0.5 GeV allowance for what the measured top mass means, so the pass rests on that allowance). It needs a top mass of 171.6 GeV against 172.60 ± 0.27, 3.7σ on the top mass alone (was 3.2σ).
  - Neutrino masses with JUNO's splittings: m₂ = 8.6 meV and sum 0.0591 eV; with the NuFIT 6.0 global fit 8.7 meV and 0.0594 eV. Both are now given as ranges (8.6–8.7 meV; 0.0591–0.0594 eV). Neutrinoless double beta decay: 1.3–4.1 meV.
  - Down quark: the Cabibbo relation gives m_s/m_d = 19.9 against 19.8 ± 0.3 from the measured masses (+0.3σ). Earlier versions compared it with an unsourced 20.4 ± 0.3 (−1.7σ).
  - Up quark: with the uncertainty of Q included, the check passes at 0.4σ with the decay value Q = 22.1 ± 0.7 and gives 2.0σ with the lattice value Q = 23.4 ± 0.5; this is now stated. It is a loose test, because Q is known to only 3%.
  - New: `analysis/higgs_three_loop_check_2026-11-18.py`. Estimated from the Standard Model's own two-to-three-loop shift (Buttazzo et al. 2013), the Higgs rule moves by −0.1 GeV, to 127.2 ± 1.5 GeV; no verdict changes. The two-loop code reproduces the published three-loop slopes to within 6%.
  - The direct measurement of the unitarity-triangle angle α (84.1°) sits 1.6σ below the 90° of the γ lead; now stated.
- The small neutrino splitting's move from 7.40 (v1.2.2) to 7.37 (v1.3.2) is explained: only the large-splitting input changed, from the global fit's 2.519 × 10⁻³ to JUNO's 2.509 × 10⁻³ eV². The fail band is now a 3σ test with the prediction's and JUNO's errors combined: outside about 7.13–7.62 at JUNO's full precision. The old band 7.32–7.48 was one standard deviation wide.
- New paper Appendix A.6 and ledger section 6: every test set after an earlier attempt failed (the grain, the far-side weight rule, equal strengths at the crossing, the compact-star pairing gap, the Higgs rule, TM1 and root signs, the 45° principle, the small-splitting ratio, the Cabibbo form and the up quark's Q). The compact-star masses moved from "passing strict tests" to "met by construction", and fission to "reproduced with fitted numbers".
- Added the four borrowed-theory tests that failed and were cited but missing: `lqg_grain`, `top_from_bottom`, `strong_cp`, `causal_set_and_e6` (all 2026-11-17).
- All 24 core scripts and the analysis scripts were rerun; every saved output was reproduced. `master_mass_table_2026-11-15.py` needs the `periodictable` package.
- The scorecard is unchanged: 7 met, 2 partly met. The input count is unchanged: 32 today, 28 if every pending measurement agrees, plus six stated choices.

## 1.3.2 (18 November 2026)

The paper was revised after three rounds of outside review (versions 1.3.0, 1.3.1 and 1.3.2 of the paper, all dated 18 November 2026). This release brings the repository to paper version 1.3.2.

- Paper (`docs/TRPG_Paper.*`) is now version 1.3.2. It cites the "all versions" DOIs, which always resolve to the latest version: code 10.5281/zenodo.24329239, paper 10.5281/zenodo.24330763.
- New checks in `analysis/`, each with its saved output:
  - `seesaw_full_matrix_2026-11-18.py`: the heavy neutrino partners from the full 3x3 seesaw with TM1 mixing (0.6–1.9 × 10⁵, 2.2–6.2 × 10⁸, 1.3–1.65 × 10¹² GeV). The earlier values (7.3 × 10⁵, 1.3 × 10⁹, 6.3 × 10¹⁰ GeV) used a shortcut that is exact only without mixing.
  - `ns_ground_state_2026-11-18.py`: the compact stars are self-bound strange quark stars; the Bodmer–Witten check passes (strange matter 820–880 MeV per baryon, below 930; up–down matter 941–1021, above 934, so ordinary nuclei are safe). Without pairing, strange matter is not the ground state.
  - `ns_running_ms_2026-11-18.py`: what the stars fix is 4Δ² − m_s² = 1.7–2.8 × 10⁴ MeV²; with a running strange mass of 120–140 MeV the gap is 88–109 MeV.
  - `crossing_vs_grain_2026-11-18.py`: the force strengths cross at 4.6 × 10¹³ GeV and separate again above it; no placement of the far-side states makes them equal at the grain (best spread 17). Equality where they cross is a stated boundary condition.
- Fitted choices: six (was five). The sixth is the coloured far-side states placed where the strengths cross, with the strengths taken as equal there. The input count is unchanged: 32 today, 28 if every pending measurement agrees.
- JUNO 2026 (arXiv:2609.37895): sin²θ₁₂ = 0.3036 ± 0.0064 puts the TM1 value 0.318 2.3σ high (3σ would fail it); the small splitting 7.388 ± 0.078 × 10⁻⁵ eV² is within 0.2σ of the lead; normal ordering is favoured at 2.1σ. `analysis/input_ledger_2026-11-18.py` now checks both.
- Neutrinoless double beta decay: 1.4–4.1 meV (was 1.3–3.9).
- No primordial gravitational waves now belongs to the majoron dark-matter lead (r below about 10⁻¹¹), not to the core.
- First-row quark-mixing unitarity is derived (the D states keep it to about 10⁻²³); a confirmed shortfall would falsify the far-side reading. A coloured spin-0 particle is added to the "nothing extra" falsifiers.
- The ledger, README and citation files are updated to match. The scorecard is unchanged: 7 met, 2 partly met.

## 1.2.2 (18 November 2026)

Bookkeeping only. No result, rule or prediction changed.

- New: `docs/ledger/TRPG_Ledger.md`, one ledger of every quantity (input, fitted, derived, predicted, lead or effective), with what each used and whether another published idea gives the same number; a "TRPG is false if…" list with the deciding experiments; and the neutrino masses as the headline prediction.
- New: `analysis/input_ledger_2026-11-18.py` (with its output), which recounts the inputs and reruns the break tests for the results that move the count (tau, down quark, lightest neutrino, Higgs) plus the strong-force strength and the JUNO lead.
- New: `figures/derivation_tree.png` and `figures/make_derivation_tree.py`, the chain from the rules to the outputs.
- The input count is now given as a range: 32 counting only savings already confirmed (tau and down quark), 28 if every pending measurement agrees (Higgs mass, lightest neutrino mass, two neutrino mixing numbers), plus five stated choices made with today's data in view.
- Accuracy: the lepton masses come from two numbers (scale and twist) with the 45° rule, not one; "Cabibbo angle from masses" and "down quark from the Cabibbo angle" are one relation, counted once; the neutrino masses are credited to the relation Brannen proposed in 2006.

## 1.2.0 (17 November 2026)

- Paper brought up to 17 November 2026: 27 measured inputs; the force-groups section rewritten (only the strengths remain); the strong coupling added as a dated prediction; results that wait on a measurement are marked "awaiting measurement" and are not counted as established until it returns.
- New folder `analysis/`: every test run from 14 to 17 November with its saved output, including the tests that failed.
- New folder `figures/`: the geometry figures and the scripts that draw them.
- New workbooks in `data/`: quark and gluon, energy, size and three-way carrier matrices.
- The website is no longer copied into this repository; it is live at https://trpgphysics.org.
- Working notes and conversation text removed from every script, output and workbook; the tests themselves are unchanged.
- ORCID added for Matthew Gelfman.
- Corrections: neutrino sum quoted as 0.0592 to 0.0594 eV, depending on which global fit of the splittings is used; Standard Model strength miss quoted as 8.9% at two loops (was 10.6% from a mixed method); the dark matter candidate's mass is not fixed (the earlier 1.4 keV and 3.4 keV values are withdrawn); quark pairing gap 82 to 97 MeV.
- Melissa Behrends-Gelfman added as co-author.
- The Higgs mass follows from a rule at the grain of space (Shaposhnikov and Wetterich, 2010): 127.2 ± 1.5 GeV with the far-side states against 125.20 ± 0.11 measured (128.8 ± 1.5 with the Standard Model's states alone). Awaiting measurement; the measured inputs go from 28 to 27. Script: `analysis/higgs_at_the_grain_2026-11-17.py`.
- No proton decay is proven for the framework's content from the tables of the 27. Script: `analysis/e6_tables_no_proton_decay_2026-11-17.py`.
- Dark energy: the main line sets its size by the grain of space and the horizon together, and it need not be constant.
- Paper title now ends "with testable predictions".
- 18 November 2026: the framework is now stated as three ideas, one rule (nothing extra) and one bet (the far side is real). The paper keeps the seven working principles under that statement.
- Inputs counted like for like: about 28, against about 34 for the Standard Model with gravity and cosmology. The earlier figure of 27 left out the neutrino mixing numbers and the two numbers for the seeds of galaxies.
- The strong coupling has one value, 0.1178: the grain of space is fixed as the Planck energy divided by √2. The choice was made with today's data in view, and the paper says so.
- Corrections: the 45° balance is stated for leptons only; flatness is consistent with the full medium, not derived from it; the Higgs rule is 3.2 standard deviations off on the top mass alone; the neutron-star prediction is the pairing gap (82 to 97 MeV), not the mass window; normal neutrino ordering holds given today's limit from cosmology (the rule's inverted solution sums to 0.103 eV); the majoron is a lead, outside the core.
- New open items: the origin of the seeds of galaxies, and where quantum probabilities come from.
- New tests: `analysis/independent_checks_2026-11-17.py` (the lepton and neutrino numbers recomputed from scratch) and `analysis/further_open_question_tests_2026-11-18.py` (five further tests, four of which fail).
- New: `docs/story/The_Full_Glass.md`, the framework told as a book, and `docs/theories/`, nine other theories of everything told the same way, with a comparison.
- ORCID added for Melissa Behrends-Gelfman.
- The matter sign (−1)^3(B−L) is stated as approximate for the far-side parts: an exact sign would leave a stable neutral far-side state, which dark-matter searches exclude (`analysis/further_open_question_tests_2026-11-18.py`, test 4).
- `analysis/leptogenesis_boltzmann_2026-11-14.py`: the closing line now carries the solver caution (results above about 1e-9 differ between solvers; settings that agree give about 2e-10, not reached). `solver_check_part30` renamed `leptogenesis_solver_check`.
- Wording leftovers removed from four scripts and one figure; no number changed.

## 1.1.1 (14 November 2026, night)

- Wording corrections in the documents and website so they say only what the theory says. Removed "the proton is the same object as a star". The website's theory page explains each idea four ways (simply, the picture, the math, for physicists).

## 1.1.0 (14 November 2026, evening)

Line-by-line check of all the math.

- Fixed: the tidal number Λ was read off a coarse grid by straight-line interpolation, which overstated it. It is now interpolated in log form on a finer grid. B = 54.1 gives Λ = 550 and passes. The 1.0.0 note saying it fails was itself wrong and is withdrawn.
- New: `proton_sets_B_2026-11-14.py` sets B = 60.1 (57.1 to 63.7) MeV/fm³ from the proton and Delta masses (cloudy bag model).
- New: `neutron_stars_no_knobs_2026-11-14.py`. With B from the proton and the measured strong force, the heaviest star is 1.54 to 1.80 solar masses (fails). With a quark pairing gap of 70 to 90 MeV it is 2.00 to 2.18, with Λ₁.₄ = 440 to 570 (passes). The gap is a prediction.
- New: `quark_link_pattern_2026-11-14.py` tests the 3/6, 4/6, 5/6 link-strength idea for quarks. It fails by about 7σ and is dropped.
- Fixed: the input count is about 25 against about 28 (the neutrino mixing and strong CP angle had been left out).
- Fixed: anomaly cancellation fixes the charges for any number of colours; 3 colours is measured, not forced.
- Fixed: antiparticle rules (colour cannot cancel against a colourless particle; different-flavour neutral quark pairs do not annihilate directly; nucleon-antinucleon pairs annihilate into pions).
- Relabelled as already known: the tau lifetime (lepton universality) and the Cabibbo relation (Gatto, Sartori and Tonin, 1968).
- Scorecard tally corrected and explained: 7 PASS, 2 PARTIAL.
- 24 scripts; `python run_all.py` gives 24 of 24 ok.

## 1.0.0 (14 November 2026)

First public release.

- 21 calculation scripts, all run by `python run_all.py` (21 of 21 pass).
- Every-particle-to-every-particle workbook, with antiparticle matrices and decays.
- Paper draft (PDF, Word and Markdown).
- Website source for trpgphysics.org.
- `neutron_star_tidal_fix_2026-11-13.py` now prints its results (its main loop had been left commented out).
- Purpose and method stated in the paper, the README and the book: the work starts from existing conclusions, discards any idea that disagreed with determined measurements, and labels choices made with today's data in view as fitted, not predicted.
- Every script was rerun (24 in `code/`, 128 in `analysis/`) and each saved output now matches its script. Fixes: leftover run markers and timing lines removed from six outputs; `ideas_for_last_four` and the leptogenesis solver check each given their own output file; the four `.json` data files that the leptogenesis scripts read are now included; `pairing_gap_scan` runs without arguments and covers every row of its saved output; one output brought up to date where a search moved in its last digits. No result changed.
- `requirements.txt` now lists Matplotlib, which the picture scripts need.

## 1.2.1 (18 November 2026)

- The lepton 45° balance is now derived for leptons: three placements of one shape make the family matrix circulant, and maximum uncertainty between self and relationship sets the self share to exactly one half, which is Koide's 2/3 and the 45° angle. It does not hold for quarks. Script: `analysis/derive_45_and_sqrt2_2026-11-18.py`.
- Which forces exist now follows from the perspective principle by a theorem: the part of Spin(9) that respects one chosen direction is exactly 8 + 3 + 1, whichever direction is chosen (Günaydin and Gürsey, 1973; Todorov and Dubois-Violette, 2018; Krasnov, 2019). The whole generation still uses one further stated reading. Script: `analysis/matrix_reading_from_rule4_2026-11-18.py`.
- Not derived, and said so: the grain factor √2 (five equally natural readings give five factors), and the rule for the far-side weak pairs' weight, which is now stated as one assumption (no direct link to the background value; Dimopoulos and Wilczek, 1981; Masiero and others, 1982).
- Paper: two reference entries that had slipped into the abstract are back in the reference list, and the list is in alphabetical order throughout.
- DOIs added to the paper, README and citation files: code, data and documents 10.5281/zenodo.24329240; the paper 10.5281/zenodo.24330764.
