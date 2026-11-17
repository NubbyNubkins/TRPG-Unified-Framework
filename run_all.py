#!/usr/bin/env python3
"""Run every TRPG script and save its output to results/.

Usage:  python run_all.py            (run everything)
        python run_all.py ringdown   (run scripts whose name contains 'ringdown')

Scripts run with data/ as the working folder, so the ones that read or rebuild
the particle workbook find it there.
"""
import os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
CODE, DATA, RES = (os.path.join(ROOT, d) for d in ("code", "data", "results"))
FIRST = ["build_relational_matrix.py", "build_antiparticle_and_decays.py"]  # workbook builders run first

def main():
    os.makedirs(RES, exist_ok=True)
    scripts = FIRST + sorted(f for f in os.listdir(CODE) if f.endswith(".py") and f not in FIRST)
    if len(sys.argv) > 1:
        scripts = [s for s in scripts if sys.argv[1] in s]
    failed = []
    for s in scripts:
        t = time.time()
        r = subprocess.run([sys.executable, os.path.join(CODE, s)], cwd=DATA, capture_output=True, text=True)
        out = r.stdout + (("\n--- errors ---\n" + r.stderr) if r.returncode else "")
        with open(os.path.join(RES, s.replace(".py", ".txt")), "w", encoding="utf-8") as f:
            f.write(out)
        status = "ok" if r.returncode == 0 else "FAILED"
        if r.returncode:
            failed.append(s)
        print(f"{status:6} {time.time()-t:6.1f}s  {s}")
    print(f"\n{len(scripts)-len(failed)} of {len(scripts)} scripts ran. Outputs are in results/.")
    sys.exit(1 if failed else 0)

if __name__ == "__main__":
    main()
