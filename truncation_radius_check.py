# -*- coding: utf-8 -*-
"""Test the fragment-fidelity diagnostics by widening the truncation radius.

The descriptor is evaluated on a radius-2 neighbourhood of the reacting core. Three
diagnostics flag the sites where that truncation is most likely to misrepresent the parent:
an alpha-prime heteroatom on the carbon flanking the carbonyl, a strained ring fused to the
enone ring, and a ring severed by the cut. The diagnostics are deterministic, but whether they
identify the sites that truncation actually perturbs is an empirical question, and this script
answers it.

For a seed-fixed, stratified subset of the released enone sites, the fragment and its enolate
adduct are rebuilt at radius 3 and recomputed with the same Hamiltonian, solvent model,
conformer protocol and acceptance criteria as the deposited radius-2 values. The quantity of
interest is the shift in the descriptor itself,

    ddE = dE(radius 3) - dE(radius 2)

in which any error common to the reactant and the adduct cancels. If the diagnostics do what
they are meant to do, flagged sites shift more than unflagged ones.

Inputs:   data/selection_pool.tsv, data/quantum_records/
Output:   out/truncation_radius_check.tsv, out/truncation_radius_summary.tsv

Run:  python truncation_radius_check.py [n_per_stratum]
Requires: numpy, pandas, rdkit, xTB on PATH or at XTB_EXE
"""
import hashlib
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.warning")
HERE = Path(__file__).parent
DATA, OUT = HERE / "data", HERE / "out"
RECORDS = Path(os.environ.get("QUANTUM_RECORDS") or DATA / "quantum_records")
JOBS = OUT / "radius3_records"
XTB = os.environ.get("XTB_EXE", "xtb")
HARTREE_KCAL = 627.5095
SEED = 42


def fragment(m, hit, radius):
    """The radius-n neighbourhood of the reacting core, aromatic rings completed."""
    keep = set(hit); front = set(hit)
    for _ in range(radius):
        front = {a.GetIdx() for i in front for a in m.GetAtomWithIdx(i).GetNeighbors()}
        keep |= front
    changed = True
    while changed:
        changed = False
        for ring in m.GetRingInfo().AtomRings():
            if (keep.intersection(ring)
                    and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                    and not set(ring) <= keep):
                keep.update(ring); changed = True
    rw = Chem.RWMol(m)
    for i in range(m.GetNumAtoms()):
        a = rw.GetAtomWithIdx(i)
        a.SetAtomMapNum(hit.index(i) + 1001 if i in hit else 0)
        if i in keep and any(n.GetIdx() not in keep for n in a.GetNeighbors()):
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            a.SetNumExplicitHs(0); a.SetNoImplicit(False)
    for i in sorted(set(range(m.GetNumAtoms())) - keep, reverse=True):
        rw.RemoveAtom(i)
    f = rw.GetMol(); Chem.SanitizeMol(f)
    return f, sorted(keep)


def product(f):
    """The enolate adduct: S-CH3 added at C-beta, the enolate written out explicitly."""
    ids = {a.GetAtomMapNum(): a.GetIdx() for a in f.GetAtoms() if a.GetAtomMapNum()}
    b, a, c, o = [ids[k] for k in range(1001, 1005)]
    rw = Chem.RWMol(f)
    for i, j, bt in [(b, a, Chem.BondType.SINGLE), (a, c, Chem.BondType.DOUBLE),
                     (c, o, Chem.BondType.SINGLE)]:
        bond = rw.GetBondBetweenAtoms(i, j)
        bond.SetBondType(bt); bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    rw.GetAtomWithIdx(o).SetFormalCharge(-1)
    rw.GetAtomWithIdx(b).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    s = rw.AddAtom(Chem.Atom("S")); me = rw.AddAtom(Chem.Atom("C"))
    rw.AddBond(b, s, Chem.BondType.SINGLE); rw.AddBond(s, me, Chem.BondType.SINGLE)
    p = rw.GetMol(); Chem.SanitizeMol(p)
    return p


def plain(m):
    m = Chem.Mol(m)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def key(smi):
    return hashlib.sha256(smi.encode()).hexdigest()[:20]


def job(smi):
    """Optimise one species from three ETKDGv3 starts; the deposited protocol, unchanged."""
    path = JOBS / key(smi); path.mkdir(parents=True, exist_ok=True)
    meta = path / "result.json"
    if meta.exists():
        return json.loads(meta.read_text(encoding="utf-8"))
    base = Chem.AddHs(Chem.MolFromSmiles(smi))
    charge = Chem.GetFormalCharge(base)
    result = dict(key=key(smi), smiles=smi, charge=charge,
                  method="GFN2-xTB/ALPB(water)", conformers=[])
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_STACKSIZE="256M")
    for seed in (42, 43, 44):
        wd = path / f"seed{seed}"; wd.mkdir(exist_ok=True)
        m = Chem.Mol(base)
        p = AllChem.ETKDGv3(); p.randomSeed = seed
        if AllChem.EmbedMolecule(m, p) != 0:
            p.useRandomCoords = True
            if AllChem.EmbedMolecule(m, p) != 0:
                result["conformers"].append(dict(seed=seed, status="embedding_failed")); continue
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:
            pass
        Chem.MolToXYZFile(m, str(wd / "input.xyz"))
        cmd = [XTB, "input.xyz", "--gfn", "2", "--opt", "tight", "--alpb", "water",
               "--chrg", str(charge), "--uhf", "0"]
        try:
            r = subprocess.run(cmd, cwd=wd, env=env, capture_output=True, text=True,
                               encoding="utf-8", errors="replace", timeout=900)
        except subprocess.TimeoutExpired:
            result["conformers"].append(dict(seed=seed, status="timeout")); continue
        (wd / "stdout.log").write_text(r.stdout, encoding="utf-8")
        energies = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
        ok = (r.returncode == 0 and "GEOMETRY OPTIMIZATION CONVERGED" in r.stdout
              and bool(energies) and (wd / "xtbopt.xyz").exists())
        status = "converged" if ok else "optimization_failed"
        if ok:
            lines = (wd / "xtbopt.xyz").read_text().splitlines()[2:]
            xyz = np.array([[float(x) for x in l.split()[1:4]]
                            for l in lines if len(l.split()) == 4])
            pt = Chem.GetPeriodicTable()
            for bd in m.GetBonds():
                i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
                rad = (pt.GetRcovalent(m.GetAtomWithIdx(i).GetAtomicNum())
                       + pt.GetRcovalent(m.GetAtomWithIdx(j).GetAtomicNum()))
                if np.linalg.norm(xyz[i] - xyz[j]) > 1.35 * rad:
                    status = "connectivity_changed"; break
        result["conformers"].append(dict(seed=seed, status=status,
                                         energy_hartree=float(energies[-1]) if energies else None))
    good = [c for c in result["conformers"] if c["status"] == "converged"]
    result["status"] = "ok" if good else "failed"
    result["energy_hartree"] = min([c["energy_hartree"] for c in good], default=None)
    meta.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def one_site(row):
    """Rebuild one site at radius 3 and return its descriptor, or the reason it could not be."""
    m = Chem.MolFromSmiles(row["parent_smiles"])
    if m is None:
        return dict(site_id=row["site_id"], status="parent unreadable")
    hit = [int(x) for x in str(row["core_atom_indices"]).split(";")]
    try:
        f3, keep3 = fragment(m, hit, 3)
        p3 = product(f3)
    except Exception as exc:
        return dict(site_id=row["site_id"], status=f"construction failed: {type(exc).__name__}")
    fs, ps = plain(f3), plain(p3)
    rf, rp = job(fs), job(ps)
    if rf["status"] != "ok" or rp["status"] != "ok":
        return dict(site_id=row["site_id"], status="no accepted energy at radius 3",
                    fragment_smiles_r3=fs)
    dE3 = (rp["energy_hartree"] - rf["energy_hartree"] - row["thiolate_hartree"]) * HARTREE_KCAL
    return dict(site_id=row["site_id"], status="ok", fragment_smiles_r3=fs,
                heavy_atoms_r2=row["heavy_atoms_r2"],
                heavy_atoms_r3=Chem.MolFromSmiles(fs).GetNumHeavyAtoms(),
                dE_r2=row["dE_kcal"], dE_r3=round(dE3, 4),
                ddE=round(dE3 - row["dE_kcal"], 4))


def main():
    n_per = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    pool = pd.read_csv(DATA / "selection_pool.tsv", sep="\t")
    ep = pool[pool.family == "enone"].copy()
    ep["flagged"] = ep.fragment_fidelity_flags != "none"
    thio = json.loads((RECORDS / ep.thiolate_job.iloc[0] / "result.json")
                      .read_text(encoding="utf-8"))["energy_hartree"]
    ep["thiolate_hartree"] = thio
    ep["heavy_atoms_r2"] = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in ep.fragment_smiles]

    rng = np.random.default_rng(SEED)
    pick = pd.concat([g.iloc[rng.permutation(len(g))[:n_per]]
                      for _k, g in ep.groupby("flagged")], ignore_index=True)
    print(f"Rebuilding {len(pick)} released enone sites at radius 3 "
          f"({int(pick.flagged.sum())} flagged, {int((~pick.flagged).sum())} unflagged)\n")

    rows = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures = {ex.submit(one_site, r._asdict()): r.site_id for r in pick.itertuples()}
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result(); rows.append(res)
            tag = (f"ddE {res['ddE']:+8.3f} kcal/mol" if res["status"] == "ok"
                   else res["status"])
            print(f"  {i:3d}/{len(pick)}  {res['site_id']:28s} {tag}")

    out = pick[["site_id", "record_id", "flagged", "fragment_fidelity_flags"]].merge(
        pd.DataFrame(rows), on="site_id", how="left")
    out.to_csv(OUT / "truncation_radius_check.tsv", sep="\t", index=False)

    ok = out[out.status == "ok"].copy()
    ok["abs_ddE"] = ok.ddE.abs()
    summary = ok.groupby("flagged").agg(
        sites=("site_id", "size"),
        median_abs_shift=("abs_ddE", "median"),
        mean_abs_shift=("abs_ddE", "mean"),
        p90_abs_shift=("abs_ddE", lambda s: float(np.quantile(s, 0.9))),
        max_abs_shift=("abs_ddE", "max"),
        median_heavy_atoms_r2=("heavy_atoms_r2", "median"),
        median_heavy_atoms_r3=("heavy_atoms_r3", "median")).round(3).reset_index()
    summary["stratum"] = np.where(summary.flagged, "carries a fidelity flag",
                                  "fidelity-clean")
    summary = summary.drop(columns="flagged")
    summary.to_csv(OUT / "truncation_radius_summary.tsv", sep="\t", index=False)
    print()
    print(summary.to_string(index=False))
    if ok.flagged.nunique() == 2:
        a = ok.loc[ok.flagged, "abs_ddE"].values
        b = ok.loc[~ok.flagged, "abs_ddE"].values
        # rank-based effect size, no distributional assumption
        delta = (np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])) / (len(a) * len(b))
        print(f"\nCliff's delta, flagged against fidelity-clean: {delta:+.3f}")


if __name__ == "__main__":
    main()
