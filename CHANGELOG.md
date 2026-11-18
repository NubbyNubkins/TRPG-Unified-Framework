# Changelog

## 1.2.0 (17 November 2026)

- Paper brought up to 17 November 2026: 27 measured inputs; the force-groups section rewritten (only the strengths remain); the strong coupling added as a dated prediction; results that wait on a measurement are marked "awaiting measurement" and are not counted as established until it returns.
- New folder `analysis/`: every test run from 14 to 17 November with its saved output, including the tests that failed.
- New folder `figures/`: the geometry figures and the scripts that draw them.
- New workbooks in `data/`: quark and gluon, energy, size and three-way carrier matrices.
- The website is no longer copied into this repository; it is live at https://trpgphysics.org.
- Working notes and conversation text removed from every script, output and workbook; the tests themselves are unchanged.
- ORCID added for Matthew Gelfman.
- Corrections: neutrino sum quoted as 0.0592 to 0.0594 eV, depending on which global fit of the splittings is used; Standard Model strength miss quoted as 8.9% at two loops (was 10.6% from a mixed method); the dark matter candidate's mass is not fixed (the earlier 1.4 keV and 3.4 keV values are withdrawn); quark pairing gap 82 to 97 MeV.
- Melissa Gelfman added as co-author.
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
- ORCID added for Melissa Gelfman.
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
- Purpose and method stated in the paper, the README and the book: the framework starts from what has been measured, discards any idea that disagreed with measurement, and labels choices made with today's data in view as fitted, not predicted.
- Every script was rerun (24 in `code/`, 128 in `analysis/`) and each saved output now matches its script. Fixes: leftover run markers and timing lines removed from six outputs; `ideas_for_last_four` and the leptogenesis solver check each given their own output file; the four `.json` data files that the leptogenesis scripts read are now included; `pairing_gap_scan` runs without arguments and covers every row of its saved output; one output brought up to date where a search moved in its last digits. No result changed.
- `requirements.txt` now lists Matplotlib, which the picture scripts need.
