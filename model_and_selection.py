"""Fit the family glutathione-reactivity calibrations and build the site-level
selection pool from the deposited fragment addition energies.

Inputs (data/):
    calibration_complete.tsv        36 calibrants with dE_kcal and log10_kGSH
    all_natural_products_energies.tsv   site-level addition energies
    anpdb_annotations.tsv           compound names
Outputs (out/):
    calibration_diagnostics.tsv, crossvalidation.tsv, enone_diagnostics.tsv,
    site_reactivity_estimates.tsv, selection_pool.tsv, site_estimate_summary.tsv

Run:  python model_and_selection.py
Requires: numpy, pandas, scipy, rdkit
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import linregress, t
from rdkit import Chem

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)
ENONE = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=[OX1]")   # Cbeta=Calpha-Ccarbonyl=O


def core_roles(mol, core):
    """Return (beta, alpha, carbonyl, oxygen) atom indices for the enone/ester core,
    identified from the substructure match rather than the stored index order."""
    target = set(core)
    for match in mol.GetSubstructMatches(ENONE):
        if set(match) == target:
            return match          # (beta, alpha, carbonyl, oxygen)
    return core[0], core[1], core[2], core[3]


def fit_families(cal):
    diag, cv = [], []
    for fam, g in cal.groupby("family"):
        x, y = g.dE_kcal.values, g.log10_kGSH.values
        frag = g.fragment_smiles.values
        fit = linregress(x, y)
        loo, grouped = [], []
        for i in range(len(g)):
            p_loo = linregress(x[np.arange(len(g)) != i], y[np.arange(len(g)) != i])
            p_grp = linregress(x[frag != frag[i]], y[frag != frag[i]])
            loo.append(p_loo.intercept + p_loo.slope * x[i])
            grouped.append(p_grp.intercept + p_grp.slope * x[i])
            cv.append(dict(calibrant_id=g.iloc[i].calibrant_id, family=fam, observed=y[i],
                           loo_prediction=loo[-1], fragment_group_prediction=grouped[-1]))
        resid = y - (fit.intercept + fit.slope * x); tss = np.sum((y - y.mean()) ** 2)
        diag.append(dict(family=fam, n=len(g), unique_fragments=len(set(frag)),
                         r2=fit.rvalue ** 2, slope=fit.slope, intercept=fit.intercept, p_value=fit.pvalue,
                         loo_q2=1 - np.sum((y - loo) ** 2) / tss, loo_rmse=np.sqrt(np.mean((y - np.array(loo)) ** 2)),
                         fragment_group_q2=1 - np.sum((y - grouped) ** 2) / tss,
                         fragment_group_rmse=np.sqrt(np.mean((y - np.array(grouped)) ** 2)),
                         residual_sd=np.sqrt(np.sum(resid ** 2) / (len(g) - 2)),
                         energy_min=x.min(), energy_max=x.max(), energy_mean=x.mean(), sxx=np.sum((x - x.mean()) ** 2)))
    return pd.DataFrame(diag), pd.DataFrame(cv)


def enone_influence(cal, d):
    g = cal[cal.family == "enone"]; x, y = g.dE_kcal.values, g.log10_kGSH.values
    n = len(x); xbar = x.mean(); sxx = np.sum((x - xbar) ** 2); s = d.loc["enone", "residual_sd"]
    yhat = d.loc["enone", "intercept"] + d.loc["enone", "slope"] * x
    h = 1 / n + (x - xbar) ** 2 / sxx
    resid = y - yhat; stud = resid / (s * np.sqrt(1 - h)); cook = (resid ** 2 / (2 * s ** 2)) * (h / (1 - h) ** 2)
    loo = {r.calibrant_id: r.loo_prediction for r in CV[CV.family == "enone"].itertuples()}
    return pd.DataFrame(dict(calibrant_id=g.calibrant_id, compound=g.compound, dE_kcal=x.round(3),
                             log10_kGSH_obs=y.round(3), log10_kGSH_fit=yhat.round(3), leverage=h.round(3),
                             studentized_residual=stud.round(3), cooks_distance=cook.round(4),
                             loo_prediction=[round(loo[c], 3) for c in g.calibrant_id]))


def domain_flags(row, d):
    m = Chem.MolFromSmiles(row.parent_smiles); core = [int(i) for i in str(row.core_atom_indices).split(";")]
    beta, alpha, carbonyl, _ = core_roles(m, core)
    flags = []
    if any(a.GetIsAromatic() for a in Chem.MolFromSmiles(row.fragment_smiles).GetAtoms()):
        flags.append("aromatic environment absent from calibration")
    if any(nb.GetAtomicNum() not in (6, 1) for i in (beta, alpha)
           for nb in m.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() not in core):
        flags.append("heteroatom substitution at alpha or beta carbon")
    if row.family == "ester" and m.GetAtomWithIdx(carbonyl).IsInRing():
        flags.append("cyclic ester absent from calibration")
    if row.family == "enone":
        rings = [len(r) for r in m.GetRingInfo().AtomRings() if beta in r and alpha in r and carbonyl in r]
        if any(sz != 5 for sz in rings):
            flags.append("ring size absent from calibration")
    lo, hi = d.loc[row.family, "energy_min"], d.loc[row.family, "energy_max"]
    inside = lo <= row.dE_kcal <= hi
    return ("within" if inside else "outside"), ("; ".join(flags) if flags else "no listed structural alert")


def fidelity_flags(row):
    """Record how faithfully the radius-2 fragment represents its parent. Reported per
    site; these describe local-model fidelity, not position within the calibration."""
    m = Chem.MolFromSmiles(row.parent_smiles)
    core = [int(i) for i in str(row.core_atom_indices).split(";")]
    keep = set(int(i) for i in str(row.retained_parent_atoms).split(";")) if isinstance(row.retained_parent_atoms, str) else set(core)
    beta, alpha, carbonyl, oxy = core_roles(m, core)
    flags = []
    alpha_prime = [n.GetIdx() for n in m.GetAtomWithIdx(carbonyl).GetNeighbors() if n.GetIdx() not in (alpha, oxy)]
    if any(nb.GetAtomicNum() not in (1, 6) for i in alpha_prime
           for nb in m.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() != carbonyl):
        flags.append("alpha-prime heteroatom")
    rings = [set(x) for x in m.GetRingInfo().AtomRings()]
    enone_rings = [x for x in rings if beta in x and alpha in x and carbonyl in x]
    if any(len(x) <= 4 and (x & e) and x is not e for e in enone_rings for x in rings):
        flags.append("strained ring fusion")
    if any((x & keep) and not x.issubset(keep) for x in rings):
        flags.append("ring severed by truncation")
    return "; ".join(flags) if flags else "none"


FAMILY_EVIDENCE = {"enone": "cross-validated calibration",
                   "ester": "weak calibration; ranked hypothesis only",
                   "enal": "no quantitative model; non-quantitative annotation"}


def use_status(family, erange, domain):
    if family == "enal":
        return "non-quantitative annotation"
    if erange == "within" and domain == "no listed structural alert":
        tag = "ranked hypothesis" if family == "ester" else "internally cross-validated candidate"
        return tag + " (in applicability domain)"
    why = []
    if erange != "within":
        why.append("energy outside calibration range")
    if domain != "no listed structural alert":
        why.append("structural alert: " + domain)
    return "extrapolation / exploratory (" + "; ".join(why) + ")"


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")
    assert len(cal) == 36 and cal.status.eq("ok").all()
    global CV
    diag, CV = fit_families(cal)
    di = diag.set_index("family")
    diag.to_csv(OUT / "calibration_diagnostics.tsv", sep="\t", index=False)
    CV.to_csv(OUT / "crossvalidation.tsv", sep="\t", index=False)
    enone_influence(cal, di).to_csv(OUT / "enone_diagnostics.tsv", sep="\t", index=False)
    print("calibration:", diag[["family", "n", "r2", "loo_q2", "fragment_group_q2"]].round(3).to_dict("records"))

    site = pd.read_csv(DATA / "all_natural_products_energies.tsv", sep="\t")
    names = pd.read_csv(DATA / "anpdb_annotations.tsv", sep="\t").set_index("molecule_id").mol_name
    rows = []
    for r in site[site.status == "ok"].itertuples():
        f = di.loc[r.family]; pred = f.intercept + f.slope * r.dE_kcal
        sem = f.residual_sd * np.sqrt(1 / f.n + (r.dE_kcal - f.energy_mean) ** 2 / f.sxx)
        pse = np.sqrt(f.residual_sd ** 2 + sem ** 2); crit = t.ppf(0.975, f.n - 2)
        erange, domain = domain_flags(r, di)
        keep = erange == "within" and domain == "no listed structural alert" and r.family != "enal"
        status = use_status(r.family, erange, domain)
        rows.append(dict(record_id=r.record_id, site_id=r.site_id, family=r.family, model_status=status,
                         family_model_evidence=FAMILY_EVIDENCE[r.family], prediction_use_status=status,
                         energy_domain_status="within calibration range" if erange == "within" else "outside calibration range",
                         structural_domain_status=domain, fragment_fidelity_flags=fidelity_flags(r),
                         compound=names.get(r.record_id, ""), parent_smiles=r.parent_smiles,
                         core_atom_indices=r.core_atom_indices, retained_parent_atoms=r.retained_parent_atoms,
                         mapped_fragment_smiles=r.mapped_fragment_smiles, fragment_smiles=r.fragment_smiles,
                         adduct_smiles=r.adduct_smiles, electrophile_job=r.electrophile_job,
                         adduct_job=r.adduct_job, thiolate_job=r.thiolate_job, status=r.status, dE_kcal=r.dE_kcal,
                         pred_log10_kGSH=pred, mean_ci_low=pred - crit * sem, mean_ci_high=pred + crit * sem,
                         prediction_interval_low=pred - crit * pse, prediction_interval_high=pred + crit * pse,
                         energy_range=erange, structural_domain=domain,
                         selection_scope="enone or ester; within energy range; no listed structural alert" if keep
                         else "chemical comparison with applicability flags"))
    res = pd.DataFrame(rows)
    groups = {f: i + 1 for i, f in enumerate(sorted(res.fragment_smiles.dropna().unique()))}
    res["fragment_group_id"] = res.fragment_smiles.map(groups)
    res.to_csv(OUT / "site_reactivity_estimates.tsv", sep="\t", index=False)
    pool = res[res.selection_scope.str.startswith("enone or ester")]
    pool.to_csv(OUT / "selection_pool.tsv", sep="\t", index=False)
    summ = res.groupby("family").agg(sites=("site_id", "size"), structures=("record_id", "nunique")).reset_index()
    summ.to_csv(OUT / "site_estimate_summary.tsv", sep="\t", index=False)
    for fam in ("enone", "ester"):
        g = pool[pool.family == fam]
        print(f"pool {fam}: {len(g)} sites / {g.record_id.nunique()} structures")
    print(f"pool total: {len(pool)} sites / {pool.record_id.nunique()} structures")


if __name__ == "__main__":
    main()
