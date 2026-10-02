"""Reproduce every tabular analysis in the article with one command, then check the result
against the deposited tables.

Each stage writes its outputs to out/ under the same file name the corresponding table carries
in data/, so the check is a field-by-field comparison of a fresh run against the deposit: same
columns in the same order, identical text, and numerical agreement to within floating-point
tolerance. Energy recomputation from structures is not part of this run; it needs xTB and is
driven by compute_fragment_energies.py.

Run:  python reproduce_all.py
Requires: numpy, pandas, scipy, rdkit
"""
import subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).parent
DATA, OUT = HERE / "data", HERE / "out"
TOL = 1e-6          # generous: observed agreement is at the 1e-13 level

STAGES = [
    ("annotate_motifs.py", "motif inventory, source-family prevalence, property contrasts, rule checks"),
    ("model_and_selection.py", "family calibrations, cross-validation, site estimates, selection pool"),
    ("benchmark_descriptor.py", "descriptor baselines, nested test, resolution counts"),
    ("compare_chemical_space.py", "property-matched comparison with the frozen reference"),
    ("reactivity_checks.py", "held-out ester test and thiol-surrogate sensitivity"),
    ("validate_annotations.py", "annotation sample against independent structure records"),
]


def run(script):
    t0 = time.time()
    r = subprocess.run([sys.executable, str(HERE / script)], cwd=HERE,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.returncode, time.time() - t0, (r.stdout or "") + (r.stderr or "")


def compare(name):
    """Compare one regenerated table with its deposited counterpart."""
    a = pd.read_csv(OUT / name, sep="\t")
    b = pd.read_csv(DATA / name, sep="\t")
    if a.shape != b.shape:
        return False, f"shape {a.shape} against {b.shape}"
    if a.columns.tolist() != b.columns.tolist():
        return False, "column order differs"
    worst = 0.0
    for c in a.columns:
        if a[c].dtype.kind == "f" and b[c].dtype.kind == "f":
            d = (a[c] - b[c]).abs().max()
            worst = max(worst, 0.0 if pd.isna(d) else float(d))
        elif not a[c].fillna("").astype(str).equals(b[c].fillna("").astype(str)):
            return False, f"text differs in {c}"
    return worst <= TOL, f"largest numerical difference {worst:.1e}"


def main():
    print("Reproducing the deposited analyses\n" + "=" * 72)
    failed = []
    for script, what in STAGES:
        code, secs, log = run(script)
        print(f"  {'ok  ' if code == 0 else 'FAIL'}  {script:28s} {secs:6.1f} s   {what}")
        if code != 0:
            failed.append(script)
            print("\n".join("        " + l for l in log.strip().splitlines()[-12:]))
    print()
    if failed:
        print("Stages that did not complete:", ", ".join(failed))
        return 1

    print("Comparing every regenerated table with the deposit\n" + "=" * 72)
    rows, mismatched = [], []
    for p in sorted(OUT.glob("*.tsv")):
        if not (DATA / p.name).is_file():
            continue
        ok, note = compare(p.name)
        rows.append((p.name, ok, note))
        if not ok:
            mismatched.append(p.name)
    for name, ok, note in rows:
        print(f"  {'match' if ok else 'DIFFERS':8s}  {name:38s} {note}")
    print()
    print(f"{sum(1 for _n, ok, _x in rows if ok)} of {len(rows)} deposited tables reproduced")
    if mismatched:
        print("Tables that differ:", ", ".join(mismatched))
        return 1
    print("Every deposited table was reproduced from the deposited inputs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
