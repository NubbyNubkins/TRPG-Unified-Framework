# Dynamic Relational Transitioning Geometry (TRPG)

**A relational, entanglement-based Theory of Everything by Matthew Gelfman and Melissa Gelfman.**
Website: https://trpgphysics.org

TRPG starts from a few plain principles: space is one completely full medium, outward push against inward pull (containment) decides what every structure becomes, no particle's scale is privileged, entanglement is fundamental, and local geometry keeps transitioning while averaging to flat. Combined with established physics, these reproduce known results and make dated predictions that experiments can check within a few years.

TRPG does not try to replace anyone. Its aim is to show that the great physicists were each correct, seeing the same thing from different perspectives. Every physicist whose work it uses is credited in the documents.

## The predictions

Each of these follows within the framework and is **awaiting measurement**. None is counted as established until its measurement returns.

| Prediction | Value | Decided by |
|---|---|---|
| Neutrino ordering and masses | normal; 0.36, 8.7, 50.1 to 50.3 meV; sum **0.0592 to 0.0594 eV** | JUNO; DESI and Euclid |
| Strong coupling at the Z | **0.1178** (0.1196 with the unreduced Planck grain) | lattice results at about 0.3% (today: 0.11873 ± 0.00056) |
| Small neutrino gap | 7.40 ± 0.06 × 10⁻⁵ eV² | JUNO at full precision |
| Neutrino mixing | sin²θ₁₂ = 0.318 | JUNO |
| Neutrino timing angle | 253–292° | DUNE, Hyper-Kamiokande |
| Quark triangle angle γ | 64.1–67.6° | LHCb, Belle II |
| Quark pairing gap in neutron stars | about 82–97 MeV | theory and merger signals |
| Neutrinoless double beta decay | 1.3–3.9 meV | beyond planned reach |
| Proton decay | none | Hyper-Kamiokande |

## Results you can rerun

| Result | TRPG | Measured | Script |
|---|---|---|---|
| Light bending at the Sun | 1.752″ | 1.75″ | `factor_two_gravity` |
| Age of the universe | 13.79 Gyr | 13.79 ± 0.02 | `expansion_and_entanglement` |
| Bag constant B from the proton and Delta masses | 60 (57–64) MeV/fm³ | 1975 bag fit ≈ 57 | `proton_sets_B` |
| Heaviest neutron star (B from proton, pairing gap 70–90 MeV) | 2.00–2.18 M☉ | ≥ 2.0 | `neutron_stars_no_knobs` |
| Tidal Λ of a 1.4 M☉ star | 440–570 | < ~580 (GW170817) | `neutron_stars_no_knobs` |
| Electron, muon, tau masses | 3 from 1 scale | 0.003% | `master_action_M2_mass_couplings` |
| Cabibbo angle from masses (already known, GST 1968) | 0.2242 | 0.2243 | `build_relational_matrix` |
| Tau lifetime from the muon (already known, lepton universality) | 289.7 fs | 290.3 fs | `build_antiparticle_and_decays` |
| Quark charges, given 3 colours | +2/3, −1/3 | same | `area_law_and_charges` |
| Black-hole ringdown (GW150914) | 250 Hz, 4.0 ms | matches | `ringdown_check` |

**Status (17 November 2026):** 7 of 9 Theory-of-Everything requirements met and 2 partly met, using 28 measured numbers.

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

All 24 core scripts run in a few minutes on a normal computer (tested with Python 3.12, NumPy 2.5, SciPy 1.18, pandas 3.1). Each script's output is saved in `results/`. See [`code/README.md`](code/README.md) for what each one shows.

## What is in this repository

| Folder | Contents |
|---|---|
| `docs/` | The paper (17 November 2026) in PDF, Word and Markdown |
| `code/` | The 24 core calculation scripts, run by `run_all.py` |
| `analysis/` | Every later test script (14 to 17 November 2026) with its saved output beside it |
| `data/` | The particle workbooks: every particle to every particle, quark and gluon, energy, size and three-way carrier matrices |
| `figures/` | The geometry figures and the scripts that draw them |
| `results/` | Saved output of the 24 core scripts |

## Use of AI tools

AI tools were used to write down the authors' ideas, to do and check the calculations, and to search the literature. The ideas are the authors'. Every result in this repository comes from a script that can be rerun.

## Cite and license

Please cite with the details in [`CITATION.cff`](CITATION.cff), or the Zenodo DOI once it exists.

- **Code and data:** MIT License ([`LICENSE`](LICENSE)).
- **Documents and website text:** © 2026 Matthew Gelfman, all rights reserved, quotation with credit welcome ([`LICENSE-DOCS.md`](LICENSE-DOCS.md)).
