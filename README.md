# Dynamic Relational Transitioning Geometry (TRPG)

**A relational, entanglement-based Theory of Everything by Matthew Gelfman and Melissa Behrends-Gelfman.**
Website: https://trpgphysics.org

TRPG is five rules and one bet (version 1.4.3; the content is that of 1.4.0, with a lattice definition of the particle side added as paper Appendix B and companion paper F).

- **Rule 1. One matrix of relationships**, seen one perspective at a time, with entanglement first. Matter is held states on the diagonal, forces are the links between them, the Higgs is the link across one state (its left side to its right side), and space and time are the geometry of the relationships ([`figures/relationship_matrix_diagonal_links.png`](figures/relationship_matrix_diagonal_links.png)).
- **Rule 2. Balance.** Whatever holds together is a pull balanced by a resistance (containment); mass is that balance.
- **Rule 3.** Within the lepton families, a particle's tie to itself balances its ties to the others (45°). For leptons it follows from the family symmetry and maximum uncertainty between self and relationship; it does not hold for quarks.
- **Rule 4. Nothing extra:** no proton decay, no fourth generation, no new forces, no mirror particles, one Higgs pair.
- **Rule 5. Smooth start.** The universe began smooth. Whether the stiffness of space then drove a short stretching (Starobinsky, 1980) is left to measurement: that alternative start would explain flatness and the uniform temperature of the sky and predicts r = 0.0040 and a ripple tilt of 0.963, which agrees with the CMB alone but is 3.2–4.2σ low once DESI is combined. The main line follows the combined data and keeps the tilt as a measured input.
- **The bet.** The far side of the relationship matrix is real. It gives the strong coupling and the Higgs-mass rule below. If the bet fails, the rules still stand.

Rules 1–4 restate the three ideas and one rule of earlier versions; rule 5 is new in 1.4.0. No particle-physics number changed in 1.4.0.

Combined with established physics, these reproduce known results and make dated predictions that experiments can check within a few years. The whole framework is told as a story anyone can read in [`docs/story/The_Full_Glass.md`](docs/story/The_Full_Glass.md), and set beside nine other theories of everything in [`docs/theories/`](docs/theories/00_Compare_the_Theories.md).

TRPG does not try to replace anyone. Its aim is to show that the great physicists were each correct, seeing the same thing from different perspectives. Every physicist whose work it uses is credited in the documents.

**Purpose and method.** This work starts from existing conclusions. Any idea that disagreed with determined measurements was discarded, and the record of each one is kept. Its purpose is to show that the Standard Model, gravity and quantum mechanics already fit together, building on the work of those who found them. Where a choice was made with today's data in view it is labelled as fitted, not predicted. Its predictions are only for measurements that have not yet been made.

## Start here: the ledger

[`docs/ledger/TRPG_Ledger.md`](docs/ledger/TRPG_Ledger.md) lists every quantity once, as input, fitted, derived or predicted, with what each one used and whether another published idea gives the same number. It also gives the input count as a range, a "TRPG is false if…" list with the deciding experiments, and the derivation tree ([`figures/derivation_tree.png`](figures/derivation_tree.png)). `analysis/input_ledger_2026-11-19.py` reruns the count and the break tests for the results that move it (tau, down quark, lightest neutrino, Higgs) plus the strong-force strength; the other numbers come from the scripts named in the paper or are standard results.

**The headline prediction is the neutrino one.** Normal ordering, masses 0.36, 8.6 to 8.7 and 50.1 to 50.3 meV, sum 0.0591 to 0.0594 eV (the range is JUNO 2026 against the NuFIT 6.0 global fit). It comes from the 45° rule that gives the tau mass, applied with one negative root, which is the relation Brannen proposed for neutrinos in 2006. JUNO, DESI and Euclid will decide it. Inverted ordering, or a sum measured clearly different from 0.059 eV, would falsify it. It is not 100% until those measurements return.

## The predictions

Each of these follows within the framework and is **awaiting measurement**. None is counted as established until its measurement returns.

| Prediction | Value | Decided by |
|---|---|---|
| Neutrino ordering and masses | normal (given today's limit from cosmology; JUNO 2026 favours normal at 2.1σ); 0.36, 8.6 to 8.7, 50.1 to 50.3 meV; sum **0.0591 to 0.0594 eV** | JUNO; DESI and Euclid |
| Strong coupling at the Z | **0.1178**, with the grain of space fixed as the Planck energy ÷ √2 (fitted, not predicted: chosen with today's data in view; the unreduced Planck energy would give 0.1196, which is not the framework's value). It also rests on the strengths being equal where they cross at 4.6 × 10¹³ GeV, a stated boundary condition: above that energy they separate again | lattice results at about 0.3% (today: 0.11873 ± 0.00056) |
| Small neutrino gap (a lead: the neutrino midpoint) | 7.37 ± 0.08 × 10⁻⁵ eV² with JUNO's large splitting (7.40 with the older global fit; only that input changed). JUNO 2026: 7.388 ± 0.078, within 0.2σ. The ratio 3 − 2√2 was proposed when the measured mass ratio was already known | JUNO at full precision |
| Neutrino mixing | sin²θ₁₂ = 0.318 (JUNO 2026: 0.3036 ± 0.0064, so 2.3σ high; 3σ would fail it) | JUNO |
| Neutrino timing angle | 253–292° | DUNE, Hyper-Kamiokande |
| Quark triangle angle γ | 64.1–67.6° | LHCb, Belle II |
| Quark pairing gap in compact (self-bound strange) stars | 4Δ² − m_s² = 1.7–2.8 × 10⁴ MeV²: about 82–97 MeV at m_s = 100 MeV, 88–109 MeV for a running m_s | theory and merger signals |
| Neutrinoless double beta decay | 1.3–4.1 meV | beyond planned reach |
| Primordial gravitational waves (alternative start, new in 1.4.0) | r = 0.0040 ± 0.0003 (limit today r < 0.036) if the alternative start phase holds; none measurable on the main line. Not seeing them at about 0.001 would rule out the alternative | Simons Observatory, LiteBIRD (early 2030s) |
| Tilt of the first ripples (new in 1.4.0) | n_s = 0.963. Under strain: −0.1σ against Planck, −2.1σ against Planck + ACT + SPT, −3.2 to −4.2σ once DESI is combined (McDonough and Ferreira, 2025) | CMB alone at higher precision; DESI–CMB agreement |
| Proton decay | none (proven for the framework's content; the grain alone gives about 10⁴⁴ years) | Hyper-Kamiokande |
| Higgs mass from the grain of space | 127.3 ± 1.5 GeV (measured 125.13 ± 0.11, PDG 2026: 1.4σ, or 2.1σ without the 0.5 GeV allowance for the top mass's meaning); needs a top mass of 171.6 GeV against 172.60 ± 0.27 (3.7σ on the top mass alone) | a sharper top-quark mass |

## Results you can rerun

| Result | TRPG | Measured | Script |
|---|---|---|---|
| Light bending at the Sun | 1.752″ | 1.75″ | `factor_two_gravity` |
| Age of the universe | 13.79 Gyr | 13.79 ± 0.02 | `expansion_and_entanglement` |
| Bag constant B from the proton and Delta masses | 60 (57–64) MeV/fm³ | 1975 bag fit ≈ 57 | `proton_sets_B` |
| Heaviest compact star (self-bound strange star; B from proton, pairing gap 82–97 MeV at m_s = 100 MeV; the gap is the prediction) | 2.0–2.15 M☉ | ≥ 2.0 | `neutron_stars_no_knobs`, `analysis/pairing_gap_scan` |
| Tidal Λ of a 1.4 M☉ star | 440–580 | < ~580 (GW170817) | `analysis/pairing_gap_scan` |
| Electron, muon, tau masses | 3 from 2 numbers (scale and twist) with the 45° rule | 0.003% | `master_action_M2_mass_couplings` |
| Cabibbo angle from masses (already known, GST 1968; the same relation gives the down quark, so it is counted once) | 0.2242 | 0.2243 | `build_relational_matrix` |
| Tau lifetime from the muon (already known, lepton universality) | 289.7 fs | 290.3 fs | `build_antiparticle_and_decays` |
| Quark charges, given 3 colours | +2/3, −1/3 | same | `area_law_and_charges` |
| Black-hole ringdown (GW150914) | 250 Hz, 4.0 ms | matches | `ringdown_check` |

**Status (25 November 2026, version 1.4.0):** on the main line the count is 32 measured numbers today and 28 if every pending measurement agrees, plus six choices made with today's data in view (33 and seven with the alternative start phase). New: the exact right-angled quark triangle as a lead; eleven further tests recorded as failed with their scripts; where each open problem stands (paper Appendix A.8). **Status (18 November 2026):** 7 of 9 Theory-of-Everything requirements met and 2 partly met. Counted like for like, the framework takes about 28 measured numbers, against about 34 for the Standard Model with gravity and cosmology. That is the figure if every pending measurement agrees; counting only savings already confirmed it is 32, plus six choices made with today's data in view (see the ledger). Saved and confirmed by measurement: the tau mass and the down quark mass. Saved only if pending measurements agree: the Higgs mass, the lightest neutrino mass and two neutrino mixing numbers. The Higgs mass follows from a rule at the grain of space (Shaposhnikov and Wetterich, 2010) and is awaiting measurement; other extra weak matter at high energy would help in the same way, so it does not single out TRPG.

- **Force groups: partly met, awaiting measurement.** Which forces exist, their charges and one-handed matter follow from the relationship matrix. Only the force strengths remain. They meet if three colourless far-side states of the matrix (named MMG) sit at 2.4 × 10⁸ GeV, which requires the strong coupling above. Not 100% until that measurement returns.
- **Quark link strengths: partly met, awaiting a rule.** Three numbers (the press on the down quarks and the two ladder openings) have no rule. No experiment closes this; it needs a new idea.
- Leads, not results: the matter excess and the dark matter candidate (a majoron wave), each with a named missing rule.

Every test carries an honest label, and the tests that failed are kept: their scripts and outputs are in `analysis/`.

## Quick start

```bash
git clone <this repository's URL>
cd TRPG-Unified-Framework
python -m pip install -r requirements.txt
python run_all.py
```

All 24 core scripts run in about 10 to 15 minutes on a normal computer (tested with Python 3.12 to 3.13, NumPy 2.5, SciPy 1.18, pandas 3.1); the three neutron-star scripts take 1 to 8 minutes each. On Windows, set `PYTHONIOENCODING=utf-8` (or run `python -X utf8`) so scripts that print Greek letters do not stop. Each script's output is saved in `results/`. See [`code/README.md`](code/README.md) for what each one shows.

## What is in this repository

| Folder | Contents |
|---|---|
| `docs/` | The paper (version 1.4.3, 27 November 2026) in PDF, Word and Markdown; `docs/story/` holds the framework told as a book; `docs/theories/` holds nine other theories of everything, each as a story, and a comparison |
| `code/` | The 24 core calculation scripts, run by `run_all.py` |
| `analysis/` | Every later test script (14 to 27 November 2026) with its saved output beside it |
| `data/` | The particle workbooks: every particle to every particle, quark and gluon, energy, size and three-way carrier matrices |
| `figures/` | The geometry figures and the scripts that draw them |
| `results/` | Saved output of the 24 core scripts |

## Use of AI tools

Due to medical necessity, AI language-model assistants were used, under the authors' direction, to draft and edit text from the authors' ideas, to do and check the calculations, and to search the literature. The authors reviewed all content and take full responsibility for it. The ideas are the authors'. Every result in this repository comes from a script that can be rerun.

## Cite and license

Please cite with the details in [`CITATION.cff`](CITATION.cff).

- **Code, data and documents (this repository), all versions:** https://doi.org/10.5281/zenodo.24329239
- **The paper, all versions:** https://doi.org/10.5281/zenodo.24330763

These DOIs always resolve to the latest version. Version 1.4.3 itself: code, data and documents [10.5281/zenodo.24578525](https://doi.org/10.5281/zenodo.24578525); paper [10.5281/zenodo.24578777](https://doi.org/10.5281/zenodo.24578777); companion paper F version 1.1 [10.5281/zenodo.24578833](https://doi.org/10.5281/zenodo.24578833); companion paper B version 1.2 [10.5281/zenodo.24578845](https://doi.org/10.5281/zenodo.24578845). Version 1.4.2: code, data and documents [10.5281/zenodo.24577279](https://doi.org/10.5281/zenodo.24577279); paper [10.5281/zenodo.24577226](https://doi.org/10.5281/zenodo.24577226). Version 1.4.1: code, data and documents [10.5281/zenodo.24552619](https://doi.org/10.5281/zenodo.24552619); paper [10.5281/zenodo.24553325](https://doi.org/10.5281/zenodo.24553325). Version 1.4.0: code [10.5281/zenodo.24551692](https://doi.org/10.5281/zenodo.24551692); paper [10.5281/zenodo.24552099](https://doi.org/10.5281/zenodo.24552099).

Companion papers (in `docs/companion_papers/`), each one component of the framework, all-versions DOIs:

- A. The colour-superconducting gap required in self-bound quark stars: [10.5281/zenodo.24553456](https://doi.org/10.5281/zenodo.24553456)
- B. Normal-ordered neutrino masses from the Koide relation: [10.5281/zenodo.24553498](https://doi.org/10.5281/zenodo.24553498)
- C. Vector-like states of the 27, α_s(M_Z) = 0.1178 and the Higgs mass: [10.5281/zenodo.24553500](https://doi.org/10.5281/zenodo.24553500)
- D. An exact right-angled unitarity triangle: [10.5281/zenodo.24553502](https://doi.org/10.5281/zenodo.24553502)
- E. No anomaly-free discrete symmetry splits the vector-like triplets and doublets: [10.5281/zenodo.24553504](https://doi.org/10.5281/zenodo.24553504)
- F. Extending a gauge-invariant lattice measure for chiral fermions to vector-like matter, gauge singlets and scalars (27 November 2026): [10.5281/zenodo.24577255](https://doi.org/10.5281/zenodo.24577255)

- **Everything in this repository (code, data, documents and website text):** © 2026 Matthew Gelfman and Melissa Behrends-Gelfman, released under the Creative Commons Attribution 4.0 International licence, CC BY 4.0 ([`LICENSE`](LICENSE); summary in [`LICENSE-DOCS.md`](LICENSE-DOCS.md)).
