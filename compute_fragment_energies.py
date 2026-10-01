"""Compute the local-fragment methanethiolate-addition energy descriptor for every
mapped alpha,beta-unsaturated carbonyl site, with persistent per-species xTB records.

For each site a radius-2 fragment is extracted around the Cbeta=Calpha-Ccarbonyl=O core
(aromatic rings completed, cut valences hydrogen-capped). The enolate adduct is built by
adding S-CH3 at Cbeta. Reactant and adduct are optimized with GFN2-xTB and ALPB water from
three ETKDGv3 seeds (42, 43, 44); the lowest accepted energy defines each species. The
descriptor is dE = E(adduct) - E(fragment) - E(CH3S-).

Usage:
    python compute_fragment_energies.py calibration          # 36 calibrants -> calibration_complete.tsv
    python compute_fragment_energies.py all_natural_products  # screened sites -> all_natural_products_energies.tsv

Inputs (data/): calibration_primary_source_manifest.tsv, anpdb_annotations.tsv
Environment: set XTB_EXE to the xtb executable (xtb 6.7.1). Requires numpy, pandas, rdkit, xtb.
"""
from pathlib import Path
import sys, os, json, hashlib, subprocess, re
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor, as_completed

RDLogger.DisableLog("rdApp.warning")
DATA = Path("data"); OUT = Path("out"); JOBS = OUT / "quantum_records"; JOBS.mkdir(parents=True, exist_ok=True)
XTB = os.environ.get("XTB_EXE", "xtb"); CORE = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=[OX1]")
QUINONE = ["O=[CX3]1[#6]=,:[#6][CX3](=O)[#6]=,:[#6]1", "O=[CX3]1[CX3](=O)[#6]=,:[#6][#6]=,:[#6]1"]
# aromatic-stabilized flavonol / 3-hydroxy-gamma-pyranone: excluded (not a Michael acceptor)
FLAVONOL = Chem.MolFromSmarts("[OX2;R]C=C([OX2H1])C(=O)")


def fragment(m, hit):
    keep = set(hit); front = set(hit)
    for _ in range(2):
        front = {a.GetIdx() for i in front for a in m.GetAtomWithIdx(i).GetNeighbors()}; keep |= front
    changed = True
    while changed:
        changed = False
        for ring in m.GetRingInfo().AtomRings():
            if keep.intersection(ring) and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring) and not set(ring) <= keep:
                keep.update(ring); changed = True
    rw = Chem.RWMol(m)
    for i in range(m.GetNumAtoms()):
        a = rw.GetAtomWithIdx(i); a.SetAtomMapNum(hit.index(i) + 1001 if i in hit else 0)
        if i in keep and any(n.GetIdx() not in keep for n in a.GetNeighbors()):
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED); a.SetNumExplicitHs(0); a.SetNoImplicit(False)
    for i in sorted(set(range(m.GetNumAtoms())) - keep, reverse=True):
        rw.RemoveAtom(i)
    f = rw.GetMol(); Chem.SanitizeMol(f)
    return Chem.MolFromSmiles(Chem.MolToSmiles(f)), sorted(keep)


def species(m):
    m = Chem.Mol(m)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def product(f):
    ids = {a.GetAtomMapNum(): a.GetIdx() for a in f.GetAtoms() if a.GetAtomMapNum()}
    b, a, c, o = [ids[k] for k in range(1001, 1005)]
    rw = Chem.RWMol(f)
    for i, j, bt in [(b, a, Chem.BondType.SINGLE), (a, c, Chem.BondType.DOUBLE), (c, o, Chem.BondType.SINGLE)]:
        bond = rw.GetBondBetweenAtoms(i, j); bond.SetBondType(bt)
        bond.SetStereo(Chem.BondStereo.STEREONONE); bond.SetBondDir(Chem.BondDir.NONE)
    rw.GetAtomWithIdx(o).SetFormalCharge(-1); rw.GetAtomWithIdx(b).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    s = rw.AddAtom(Chem.Atom("S")); me = rw.AddAtom(Chem.Atom("C"))
    rw.AddBond(b, s, Chem.BondType.SINGLE); rw.AddBond(s, me, Chem.BondType.SINGLE)
    p = rw.GetMol(); Chem.SanitizeMol(p)
    return p


def family(m, hit):
    b, a, c, o = hit
    ns = [n for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() not in (a, o)]
    if not ns:
        return "enal"
    if len(ns) != 1:
        return None
    if ns[0].GetAtomicNum() == 8 and ns[0].GetTotalNumHs() > 0:
        return None
    return {6: "enone", 8: "ester"}.get(ns[0].GetAtomicNum())


def build_maps(which):
    rows, fails = [], []
    if which == "calibration":
        inp = pd.read_csv(DATA / "calibration_primary_source_manifest.tsv", sep="\t")
        items = [(r.calibrant_id, r.canonical_smiles, r.family) for r in inp.itertuples()]
    else:
        ann = pd.read_csv(DATA / "anpdb_annotations.tsv", sep="\t").set_index("molecule_id")
        items = [(i, ann.loc[i, "parent_smiles"], None) for i in ann.index if i != "Mol_00042"]
    for mid, smi, expected in items:
        m = Chem.MolFromSmiles(smi); seen = set(); quinone = set()
        flavonol = {a for fh in m.GetSubstructMatches(FLAVONOL) for a in fh}
        if not expected:
            for patt in QUINONE:
                for qh in m.GetSubstructMatches(Chem.MolFromSmarts(patt)):
                    quinone.update(qh)
        for hit in m.GetSubstructMatches(CORE):
            if not expected and (hit[2] in quinone or set(hit) & flavonol):
                continue
            fam = family(m, hit)
            if fam is None or (expected and fam != expected) or hit[0] in seen:
                continue
            seen.add(hit[0])
            try:
                f, keep = fragment(m, hit); p = product(f)
                rows.append(dict(record_id=mid, site_id=f"{mid}:atom{hit[0]}", family=fam, parent_smiles=smi,
                                 core_atom_indices=";".join(map(str, hit)), retained_parent_atoms=";".join(map(str, keep)),
                                 mapped_fragment_smiles=Chem.MolToSmiles(f), fragment_smiles=species(f), adduct_smiles=species(p)))
            except Exception as e:
                fails.append(dict(record_id=mid, core=str(hit), reason=str(e)))
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / f"{which}_site_maps.tsv", sep="\t", index=False)
    pd.DataFrame(fails, columns=["record_id", "core", "reason"]).to_csv(OUT / f"{which}_fragment_exclusions.tsv", sep="\t", index=False)
    return frame


def key(smi):
    return hashlib.sha256(smi.encode()).hexdigest()[:20]


def job(smi):
    path = JOBS / key(smi); path.mkdir(exist_ok=True); meta = path / "result.json"
    if meta.exists():
        return json.loads(meta.read_text(encoding="utf-8"))
    base = Chem.AddHs(Chem.MolFromSmiles(smi)); charge = Chem.GetFormalCharge(base)
    result = dict(key=key(smi), smiles=smi, charge=charge, uhf=0, method="GFN2-xTB/ALPB(water)", conformers=[])
    env = os.environ.copy(); env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_STACKSIZE="256M")
    for seed in (42, 43, 44):
        wd = path / f"seed{seed}"; wd.mkdir(exist_ok=True); m = Chem.Mol(base)
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
        cmd = [XTB, "input.xyz", "--gfn", "2", "--opt", "tight", "--alpb", "water", "--chrg", str(charge), "--uhf", "0"]
        (wd / "command.json").write_text(json.dumps(cmd))
        try:
            r = subprocess.run(cmd, cwd=wd, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180)
        except subprocess.TimeoutExpired:
            result["conformers"].append(dict(seed=seed, status="timeout")); continue
        (wd / "stdout.log").write_text(r.stdout, encoding="utf-8")
        energies = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
        ok = r.returncode == 0 and "GEOMETRY OPTIMIZATION CONVERGED" in r.stdout and bool(energies) and (wd / "xtbopt.xyz").exists()
        status = "converged" if ok else "optimization_failed"
        if ok:
            lines = (wd / "xtbopt.xyz").read_text().splitlines()[2:]
            xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines if len(l.split()) == 4])
            pt = Chem.GetPeriodicTable()
            for b in m.GetBonds():
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                radii = pt.GetRcovalent(m.GetAtomWithIdx(i).GetAtomicNum()) + pt.GetRcovalent(m.GetAtomWithIdx(j).GetAtomicNum())
                if np.linalg.norm(xyz[i] - xyz[j]) > 1.35 * radii:
                    status = "connectivity_changed"; break
        result["conformers"].append(dict(seed=seed, status=status, energy_hartree=float(energies[-1]) if energies else None))
    good = [c for c in result["conformers"] if c["status"] == "converged"]
    result["status"] = "ok" if good else "failed"
    result["energy_hartree"] = min([c["energy_hartree"] for c in good], default=None)
    meta.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "calibration"
    frame = build_maps(which)
    smis = sorted(set(frame.fragment_smiles) | set(frame.adduct_smiles) | {"C[S-]"})
    print(which, "site records", len(frame), "unique species", len(smis), flush=True)
    with ProcessPoolExecutor(max_workers=6) as pool:
        futs = {pool.submit(job, s): s for s in smis}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if i % 25 == 0 or r["status"] != "ok":
                print(i, "/", len(smis), r["status"], flush=True)
    energies = {s: json.loads((JOBS / key(s) / "result.json").read_text()) for s in smis}
    thiolate = energies["C[S-]"]["energy_hartree"]
    rows = []
    for r in frame.to_dict("records"):
        e, p = energies[r["fragment_smiles"]], energies[r["adduct_smiles"]]
        ok = e["status"] == p["status"] == "ok" and thiolate is not None
        r.update(electrophile_job=key(r["fragment_smiles"]), adduct_job=key(r["adduct_smiles"]),
                 thiolate_job=key("C[S-]"), status="ok" if ok else "failed",
                 dE_kcal=627.5095 * (p["energy_hartree"] - e["energy_hartree"] - thiolate) if ok else None)
        rows.append(r)
    energies_frame = pd.DataFrame(rows)
    energies_frame.to_csv(OUT / f"{which}_energies.tsv", sep="\t", index=False)
    if which == "calibration":
        manifest = pd.read_csv(DATA / "calibration_primary_source_manifest.tsv", sep="\t")
        merged = manifest.merge(energies_frame, left_on="calibrant_id", right_on="record_id", suffixes=("", "_site"))
        merged.to_csv(OUT / "calibration_complete.tsv", sep="\t", index=False)
    print("completed", which, flush=True)


if __name__ == "__main__":
    main()
