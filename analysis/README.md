# Later tests (14 to 18 November 2026)

Each `.py` file here is one test, and the `.txt` file with the same name is its saved output. Tests that failed are kept beside the ones that passed. Where a script says "an earlier test", it means another script in this folder.

These are not run by `run_all.py`: a few take several minutes, and some need SciPy's optimisers. To rerun one:

```bash
python analysis/<name>.py
```

Scripts whose names start with `build_` or `add_` rebuild the workbooks in `data/` and expect to be run from that folder.
