# -*- coding: utf-8 -*-
"""Benchmark the addition energy against established electrophilicity descriptors.

The addition energy used in this work is one of several quantum-chemical quantities that have
been used to rank Michael acceptors. This script computes four of the conventional ones on the
same optimised fragment geometries, so the comparison is made on identical structures and
identical cross-validation splits rather than against literature values obtained some other way.

For each calibration electrophile the lowest-energy accepted conformer is taken from the
deposited calculation record and a single point is run with the same Hamiltonian and solvent
model, requesting condensed Fukui indices. From that single point:

    LUMO energy                 the frontier-orbital criterion in its simplest form
    global electrophilicity     omega = mu^2 / (2 eta), with mu = (E_HOMO + E_LUMO) / 2
                                and eta = E_LUMO - E_HOMO
    Fukui f(+) at C-beta        the condensed electrophilic-attack index at the reacting atom
    local electrophilicity      omega_k = omega * f(+) at C-beta

Each is then fitted as the single independent variable of an ordinary least-squares model
within each reaction family and cross-validated by withholding whole fragment groups, which is
the scheme the accompanying article reports, so the numbers sit beside the addition energy in
the same table.

Needs xTB and the deposited calculation records, so it is not part of reproduce_all.py. Its
output table is deposited, and the benchmark reported in the article is read from that table.

Descriptors are emitted for the calibration electrophiles and for every distinct
electrophile fragment of the released set, so a user who prefers one of these indices to
the addition energy can apply it to the released sites directly.

Inputs:   data/calibration_complete.tsv, data/selection_pool.tsv, data/quantum_records/
Output:   out/electrophilicity_comparison.tsv

Run:  python electrophilicity_indices.py
Requires: numpy, pandas, scipy, rdkit, xTB on PATH or at XTB_EXE
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import linregress

HERE = Path(__file__).parent
DATA, OUT = HERE / "data", HERE / "out"
# The calculation records are bulky; QUANTUM_RECORDS points at them wherever they are
# unpacked, and the deposited layout places them at data/quantum_records.
RECORDS = Path(os.environ.get("QUANTUM_RECORDS") or DATA / "quantum_records")
OUT.mkdir(exist_ok=True)
XTB = os.environ.get("XTB_EXE") or shutil.which("xtb") or "xtb"
HARTREE_EV = 27.211386245988

# C-beta is the alkene carbon remote from the carbonyl; atom 0 of this pattern.
ENONE = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=[OX1]")


def best_conformer(key):
    """The seed directory holding the lowest accepted energy for one calculation record."""
    r = json.loads((RECORDS / key / "result.json").read_text(encoding="utf-8"))
    done = [c for c in r["conformers"] if c["status"] == "converged"]
    if not done:
        return None, r
    return RECORDS / key / f"seed{min(done, key=lambda c: c['energy_hartree'])['seed']}", r


def beta_index(seed_dir):
    """Index of C-beta in the atom order the geometry file uses."""
    m = Chem.MolFromMolFile(str(seed_dir / "input.mol"), removeHs=False, sanitize=False)
    if m is None:
        return None
    Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    hits = m.GetSubstructMatches(ENONE)
    if not hits:
        return None
    elements = [a.GetSymbol() for a in m.GetAtoms()]
    xyz = (seed_dir / "xtbopt.xyz").read_text(encoding="utf-8").splitlines()
    if [l.split()[0] for l in xyz[2:2 + len(elements)]] != elements:
        return None                     # atom order must agree, or the index means nothing
    return hits[0][0]


def single_point(xyz, charge):
    """One GFN2-xTB/ALPB(water) single point with condensed Fukui indices."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        shutil.copy(xyz, work / "in.xyz")
        p = subprocess.run([XTB, "in.xyz", "--gfn", "2", "--alpb", "water",
                            "--chrg", str(charge), "--vfukui"],
                           cwd=work, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=900)
        return p.stdout + p.stderr


def parse(log):
    """Frontier-orbital energies in eV and the condensed f(+) column."""
    homo = lumo = None
    for line in log.splitlines():
        if "(HOMO)" in line:
            homo = float(line.split()[-2])
        elif "(LUMO)" in line and lumo is None:
            lumo = float(line.split()[-2])
    fplus, reading = [], False
    for line in log.splitlines():
        if re.search(r"#\s+f\(\+\)\s+f\(-\)\s+f\(0\)", line):
            reading = True
            continue
        if reading:
            m = re.match(r"\s*(\d+)([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
                         line)
            if m:
                fplus.append(float(m.group(3)))
            elif fplus:
                break
    return homo, lumo, fplus


def group_q2(x, y, groups):
    """Leave-one-fragment-group-out predictive Q2, the scheme the article reports."""
    pred = np.empty(len(y))
    for i in range(len(y)):
        k = groups != groups[i]
        if k.sum() < 3:
            pred[i] = y[k].mean() if k.any() else y.mean()
            continue
        f = linregress(x[k], y[k])
        pred[i] = f.intercept + f.slope * x[i]
    q2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return q2, float(np.sqrt(np.mean((y - pred) ** 2)))


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")
    pool = pd.read_csv(DATA / "selection_pool.tsv", sep=chr(9))
    cal_keys = set(cal.electrophile_job.dropna())
    pool_keys = set(pool.electrophile_job.dropna())
    rows, skipped = [], []
    for key in sorted(cal_keys | pool_keys):
        seed, rec = best_conformer(key)
        if seed is None:
            skipped.append((key, "no accepted conformer"))
            continue
        b = beta_index(seed)
        log = single_point(seed / "xtbopt.xyz", rec["charge"])
        homo, lumo, fplus = parse(log)
        if homo is None or lumo is None or not fplus:
            skipped.append((key, "single point did not report the expected properties"))
            continue
        mu, eta = (homo + lumo) / 2.0, lumo - homo
        omega = mu ** 2 / (2 * eta)
        fb = fplus[b] if (b is not None and b < len(fplus)) else np.nan
        rows.append(dict(electrophile_job=key, smiles=rec["smiles"], seed=seed.name,
                         scope="calibration" if key in cal_keys and key in pool_keys
                         else ("calibration" if key in cal_keys else "released set"),
                         homo_ev=round(homo, 4), lumo_ev=round(lumo, 4),
                         hardness_ev=round(eta, 4), omega_ev=round(omega, 4),
                         fukui_plus_beta=round(fb, 4) if fb == fb else np.nan,
                         omega_local_beta=round(omega * fb, 4) if fb == fb else np.nan))
        print(f"  {key}  LUMO {lumo:7.3f} eV   omega {omega:6.3f} eV   f+(Cb) "
              f"{fb:6.3f}" if fb == fb else f"  {key}  LUMO {lumo:7.3f} eV")
    idx = pd.DataFrame(rows)
    idx.to_csv(OUT / "electrophilicity_descriptors.tsv", sep="\t", index=False)
    for key, why in skipped:
        print(f"  skipped {key}: {why}")

    j = cal.merge(idx, on="electrophile_job", how="left")
    out = []
    for fam, g in j.groupby("family"):
        g = g.reset_index(drop=True)
        y, grp = g.log10_kGSH.values, g.fragment_smiles.values
        for label, col in [("addition energy", "dE_kcal"), ("LUMO energy", "lumo_ev"),
                           ("global electrophilicity", "omega_ev"),
                           ("Fukui f(+) at C-beta", "fukui_plus_beta"),
                           ("local electrophilicity at C-beta", "omega_local_beta")]:
            x = g[col].values.astype(float)
            ok = ~np.isnan(x)
            if ok.sum() < len(x):
                out.append(dict(family=fam, descriptor=label, n=int(ok.sum()),
                                fragments=int(pd.Series(grp[ok]).nunique()),
                                r2=np.nan, fragment_group_q2=np.nan, fragment_group_rmse=np.nan,
                                note="not defined for every calibrant"))
                continue
            fit = linregress(x, y)
            q2, rmse = group_q2(x, y, grp)
            out.append(dict(family=fam, descriptor=label, n=len(x),
                            fragments=int(pd.Series(grp).nunique()),
                            r2=round(fit.rvalue ** 2, 3), fragment_group_q2=round(q2, 3),
                            fragment_group_rmse=round(rmse, 3), note=""))
    cmp = pd.DataFrame(out)
    cmp.to_csv(OUT / "electrophilicity_comparison.tsv", sep="\t", index=False)
    print()
    print(cmp.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
