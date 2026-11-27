# Later tests (14 to 18 November 2026)

Each `.py` file here is one test, and the `.txt` file with the same name is its saved output. Tests that failed are kept beside the ones that passed. Where a script says "an earlier test", it means another script in this folder.

These are not run by `run_all.py`: a few take several minutes, and some need SciPy's optimisers. To rerun one:

```bash
python analysis/<name>.py
```

The leptogenesis scripts (`leptogenesis_*`, `final_pass_*`, `ulysses_check_*`, `matter_excess_*`, `falsification_matter_excess_*`, `matter_and_dark_matter_final_*`) read one another and the saved `.json` files beside them, so run those from inside this folder; four of them also need the ULYSSES package (Granelli et al., 2021; `pip install ulysses`). The picture scripts need Matplotlib. `ideas_for_last_four` is in two parts, each with its own output. Two outputs, `force_pull_on_side_2026-11-16.txt` and `one_continuous_spiral_2026-11-17.txt`, are short calculations kept as tables without a separate script. Searches that start from random points can move in the last digit or two between library versions.

Scripts whose names start with `build_` or `add_` rebuild the workbooks in `data/` and expect to be run from that folder.

Four checks added for version 1.3.2 (18 November 2026): `seesaw_full_matrix` (heavy neutrino partners from the full 3x3 seesaw), `ns_ground_state` (the compact stars as self-bound strange stars, with the Bodmer–Witten check), `ns_running_ms` (the gap with a running strange-quark mass; it reads `ns_ground_state`) and `crossing_vs_grain` (whether the force strengths can be equal at the grain; it reads `force_strengths_final_2026-11-15.py`). `ns_ground_state` reads `neutron_stars_no_knobs_2026-11-14.py`; each of these finds its neighbours in this folder.

Added for version 1.3.3 (18 November 2026): `recalibrate_to_2026_values` recomputes every headline number from the Particle Data Group's 2026 values and JUNO 2026 (it reads `higgs_at_the_grain_2026-11-17.py`), `higgs_three_loop_check` (the Higgs rule with three-loop running; it reads the same script), and four borrowed-theory tests that failed: `lqg_grain` (loop quantum gravity's smallest area as the grain), `top_from_bottom` (the top mass from asymptotic safety), `strong_cp` (Nelson–Barr and the majoron as an axion) and `causal_set_and_e6` (causal-set dark energy and the E6 tables). `master_mass_table_2026-11-15.py` needs the `periodictable` package (`pip install periodictable`).

Added for version 1.4.2 (27 November 2026), the lattice checks of paper Appendix B and companion paper F: `lattice_fit_check` (anomaly sums, Craig's even-count condition, one-loop coefficients), `vectorlike_pair_measure_check` (Lemma A on a 3⁴ SU(3) lattice, with a two-dimensional control; it takes about a minute), `bl_exact_on_lattice` (Lemma B: the index on a two-dimensional U(1) lattice and the B − L sums) and `extras_theta_baryon_alternatives` (baryon number in steps of three; which pictures Craig's construction covers). They need only NumPy. `craig_proof_recheck` is an independent recheck of the computable steps of Craig's proof (colour-trace cancellation, his content tests, free-field facts, his Lemma 5.1, the global-step algebra; 18 checks, about five minutes); it checks steps, it does not re-prove the theorem.
