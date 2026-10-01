"""Benchmark the local-fragment addition energy against simpler descriptors.

Three comparators are fitted on identical fragment-group splits, so that every model is
asked the same question: predict a compound whose local fragment was withheld from training.

    family mean        no structural information
    structural counts  carbon substituents at Calpha and Cbeta, plus a ring indicator
    fingerprint        256-bit Morgan fingerprint of the local fragment, ridge regression
    addition energy    the dE descriptor used in the article

The script also measures how finely each descriptor separates the released selection pool,
which is the property that governs site-level ranking.

Inputs (data/): calibration_complete.tsv, site_reactivity_estimates.tsv
Output (out/):  descriptor_baseline_comparison.tsv, descriptor_resolution.tsv

Run:  python benchmark_descriptor.py
Requires: numpy, pandas, rdkit
"""
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)
CORE = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=[OX1]")


def core_roles(mol, core):
    """Return (beta, alpha, carbonyl, oxygen) from the substructure match."""
    target = set(core)
    for match in mol.GetSubstructMatches(CORE):
        if set(match) == target:
            return match
    return None


def structural_features(parent_smiles, core_atom_indices):
    """Carbon substituents at Calpha and Cbeta outside the core, and a ring indicator."""
    m = Chem.MolFromSmiles(parent_smiles)
    if m is None:
        return None
    core = [int(i) for i in str(core_atom_indices).split(";")]
    roles = core_roles(m, core)
    if roles is None:
        return None
    beta, alpha, carbonyl, _ = roles
    count = lambda i: sum(1 for n in m.GetAtomWithIdx(i).GetNeighbors()
                          if n.GetAtomicNum() == 6 and n.GetIdx() not in core)
    cyclic = int(any(beta in r and alpha in r and carbonyl in r for r in m.GetRingInfo().AtomRings()))
    return count(alpha), count(beta), cyclic


def grouped_predictions(X, y, groups, ridge=0.0):
    """Leave-one-fragment-group-out predictions for a linear model."""
    pred = np.full(len(y), np.nan)
    for i in range(len(y)):
        train = groups != groups[i]
        Xt, yt = X[train], y[train]
        if ridge:
            beta = np.linalg.solve(Xt.T @ Xt + ridge * np.eye(X.shape[1]), Xt.T @ (yt - yt.mean()))
            pred[i] = X[i] @ beta + yt.mean()
        elif train.sum() >= X.shape[1]:
            beta, *_ = np.linalg.lstsq(Xt, yt, rcond=None)
            pred[i] = X[i] @ beta
        else:
            pred[i] = yt.mean()
    return pred


def q2_rmse(y, pred):
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2), np.sqrt(np.mean((y - pred) ** 2))


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")
    rows = []
    for family, g in cal.groupby("family"):
        g = g.reset_index(drop=True)
        y = g.log10_kGSH.values; groups = g.fragment_smiles.values; ones = np.ones((len(g), 1))
        feats = np.array([structural_features(r.parent_smiles, r.core_atom_indices) for r in g.itertuples()])
        fp = np.array([list(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=256))
                       for s in g.fragment_smiles], dtype=float)
        mean_pred = np.array([y[groups != groups[i]].mean() for i in range(len(y))])
        models = {
            "family mean": (mean_pred, 0),
            "structural counts": (grouped_predictions(np.column_stack([feats, ones]), y, groups), 4),
            "fingerprint ridge": (grouped_predictions(fp, y, groups, ridge=1.0), 256),
            "addition energy": (grouped_predictions(np.column_stack([g.dE_kcal.values, ones]), y, groups), 2),
        }
        for name, (pred, ncoef) in models.items():
            q2, rmse = q2_rmse(y, pred)
            rows.append(dict(family=family, n=len(g), fragments=len(set(groups)), descriptor=name,
                             coefficients=ncoef, fragment_group_q2=round(q2, 3), fragment_group_rmse=round(rmse, 3)))
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "descriptor_baseline_comparison.tsv", sep="\t", index=False)
    print(comp.pivot(index="family", columns="descriptor", values="fragment_group_q2").to_string())

    # How finely does each descriptor separate the sites it is deployed on?
    est = pd.read_csv(DATA / "site_reactivity_estimates.tsv", sep="\t")
    res = []
    for label, sub in [("all enone sites", est[est.family == "enone"]),
                       ("enone selection pool", est[(est.family == "enone") &
                        (est.prediction_use_status.str.contains("in applicability domain"))])]:
        feats = [structural_features(r.parent_smiles, r.core_atom_indices) for r in sub.itertuples()]
        feats = [f for f in feats if f]
        counts = pd.Series(feats).value_counts()
        res.append(dict(site_set=label, sites=len(sub), distinct_structural_states=len(counts),
                        largest_structural_group=int(counts.iloc[0]),
                        distinct_addition_energies=int(sub.dE_kcal.round(3).nunique())))
    resolution = pd.DataFrame(res)
    resolution.to_csv(OUT / "descriptor_resolution.tsv", sep="\t", index=False)
    print()
    print(resolution.to_string(index=False))


if __name__ == "__main__":
    main()
