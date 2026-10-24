"""External-validation metrics, cross-validation optimism, and an applicability-domain check.

Three analyses, each written to its own table.

1. Cross-validation optimism. Fragment-based descriptors repeat: several compounds in a
   congeneric series can share one local fragment. Leave-one-out then predicts a compound
   partly from a near-duplicate of itself. Comparing leave-one-out with
   leave-one-fragment-group-out quantifies how far the naive statistic overstates
   predictivity, and the fragment-group size distribution shows why.

2. External-validation metrics for the independent same-assay test set, in the forms
   conventional for QSAR reporting: Q2F1, Q2F2, Q2F3, the concordance correlation
   coefficient, and the r2m family. Because tied compounds must receive one prediction,
   the table also carries the lowest root mean square error any model on this descriptor
   can reach on this set, which separates descriptor resolution from calibration error.

3. Held-out error by leverage stratum. Each family is split at its median leverage and the
   leave-one-out error of each stratum is reported, so the per-site leverage released with
   every prediction can be read against measured held-out accuracy.

Inputs (data/):  calibration_complete.tsv, external_reactivity_predictions.tsv
Outputs (out/):  crossvalidation_optimism.tsv, validation_metrics.tsv,
                 applicability_performance.tsv

Run:  python validation_metrics.py
Requires: numpy, pandas, scipy
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import linregress

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)


def loo_predictions(x, y):
    out = np.empty(len(y))
    for i in range(len(y)):
        k = np.arange(len(y)) != i
        f = linregress(x[k], y[k])
        out[i] = f.intercept + f.slope * x[i]
    return out


def group_predictions(x, y, groups):
    out = np.empty(len(y))
    for i in range(len(y)):
        k = groups != groups[i]
        if k.sum() < 3:
            out[i] = y[k].mean() if k.any() else y.mean()
            continue
        f = linregress(x[k], y[k])
        out[i] = f.intercept + f.slope * x[i]
    return out


def q2(y, pred, ref=None):
    ref = y.mean() if ref is None else ref
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - ref) ** 2)


def ccc(a, b):
    """Lin's concordance correlation coefficient."""
    va, vb = a.var(ddof=0), b.var(ddof=0)
    return 2 * np.cov(a, b, ddof=0)[0, 1] / (va + vb + (a.mean() - b.mean()) ** 2)


def rm2(obs, prd):
    """Roy's r2m and its companion, from the fitted and through-origin correlations."""
    r2 = np.corrcoef(obs, prd)[0, 1] ** 2
    k = np.sum(obs * prd) / np.sum(prd ** 2)          # observed regressed on predicted
    r02 = 1 - np.sum((obs - k * prd) ** 2) / np.sum((obs - obs.mean()) ** 2)
    kp = np.sum(obs * prd) / np.sum(obs ** 2)
    r02p = 1 - np.sum((prd - kp * obs) ** 2) / np.sum((prd - prd.mean()) ** 2)
    m = r2 * (1 - np.sqrt(abs(r2 - r02)))
    mp = r2 * (1 - np.sqrt(abs(r2 - r02p)))
    return r2, r02, m, mp, (m + mp) / 2, abs(m - mp)


def cluster_bootstrap(x, y, groups, draws=2000, seed=42):
    """Percentile intervals for the fit statistics, resampling whole fragment groups.

    Compounds sharing a local fragment are not independent observations, so the resampling
    unit is the fragment group rather than the compound. Draws that leave fewer than three
    distinct descriptor values cannot support a fit and are not counted.
    """
    rng = np.random.default_rng(seed)
    names = np.array(sorted(set(groups)))
    idx = {g: np.flatnonzero(groups == g) for g in names}
    r2s, q2s = [], []
    for _ in range(draws):
        take = np.concatenate([idx[g] for g in rng.choice(names, len(names), replace=True)])
        xb, yb, gb = x[take], y[take], groups[take]
        if len(np.unique(np.round(xb, 6))) < 3:
            continue
        r2s.append(linregress(xb, yb).rvalue ** 2)
        q2s.append(q2(yb, group_predictions(xb, yb, gb)))
    lo_r2, hi_r2 = np.percentile(r2s, [2.5, 97.5])
    lo_q2, hi_q2 = np.percentile(q2s, [2.5, 97.5])
    return dict(usable_draws=len(r2s), r2_ci_low=round(float(lo_r2), 3),
                r2_ci_high=round(float(hi_r2), 3),
                fragment_group_q2_ci_low=round(float(lo_q2), 3),
                fragment_group_q2_ci_high=round(float(hi_q2), 3))


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")

    # --- 1. cross-validation optimism -------------------------------------------------
    rows = []
    for fam, g in cal.groupby("family"):
        g = g.reset_index(drop=True)
        x, y, grp = g.dE_kcal.values, g.log10_kGSH.values, g.fragment_smiles.values
        sizes = pd.Series(grp).value_counts()
        q_loo = q2(y, loo_predictions(x, y))
        q_grp = q2(y, group_predictions(x, y, grp))
        rows.append(dict(
            family=fam, n=len(g), distinct_fragments=int(sizes.size),
            largest_fragment_group=int(sizes.iloc[0]),
            compounds_sharing_a_fragment=int(sizes[sizes > 1].sum()),
            fragment_redundancy=round(1 - sizes.size / len(g), 3),
            loo_q2=round(q_loo, 3), fragment_group_q2=round(q_grp, 3),
            optimism=round(q_loo - q_grp, 3),
            loo_rmse=round(float(np.sqrt(np.mean((y - loo_predictions(x, y)) ** 2))), 3),
            fragment_group_rmse=round(float(np.sqrt(np.mean((y - group_predictions(x, y, grp)) ** 2))), 3)))
    opt = pd.DataFrame(rows).sort_values("optimism", ascending=False)
    opt.to_csv(OUT / "crossvalidation_optimism.tsv", sep="\t", index=False)
    print("Cross-validation optimism")
    print(opt.to_string(index=False))

    # --- 2. external-validation metrics -----------------------------------------------
    ext = pd.read_csv(DATA / "external_reactivity_predictions.tsv", sep="\t")
    obs = ext.measured_log10_kGSH.values.astype(float)
    prd = ext.predicted_log10_kGSH.values.astype(float)
    tr = cal[cal.family == "ester"].log10_kGSH.values
    r2, r02, m, mp, mbar, mdelta = rm2(obs, prd)
    # Tied compounds share one prediction, so their within-tie scatter cannot be fitted away.
    # The best constant predictor for a tie is its mean, which fixes the attainable minimum.
    ties = ext.dE_kcal.round(3)
    floor = float(np.sqrt(sum(float(((g - g.mean()) ** 2).sum())
                              for _k, g in ext.measured_log10_kGSH.groupby(ties)) / len(obs)))
    mets = pd.DataFrame([dict(
        set="independent same-assay esters", n=len(obs),
        distinct_descriptor_values=int(ext.dE_kcal.round(3).nunique()),
        largest_tied_group=int(ext.dE_kcal.round(3).value_counts().iloc[0]),
        response_spread_within_largest_tie=round(float(
            obs[ext.dE_kcal.round(3) == ext.dE_kcal.round(3).value_counts().index[0]].max() -
            obs[ext.dE_kcal.round(3) == ext.dE_kcal.round(3).value_counts().index[0]].min()), 3),
        rmse=round(float(np.sqrt(np.mean((obs - prd) ** 2))), 3),
        rmse_floor_from_ties=round(floor, 3),
        mae=round(float(np.mean(np.abs(obs - prd))), 3),
        q2_f1=round(float(1 - np.sum((obs - prd) ** 2) / np.sum((obs - tr.mean()) ** 2)), 3),
        q2_f2=round(float(1 - np.sum((obs - prd) ** 2) / np.sum((obs - obs.mean()) ** 2)), 3),
        q2_f3=round(float(1 - (np.sum((obs - prd) ** 2) / len(obs)) /
                          (np.sum((tr - tr.mean()) ** 2) / len(tr))), 3),
        ccc=round(float(ccc(obs, prd)), 3),
        r2=round(float(r2), 3), r2_through_origin=round(float(r02), 3),
        rm2=round(float(m), 3), rm2_reverse=round(float(mp), 3),
        rm2_mean=round(float(mbar), 3), rm2_delta=round(float(mdelta), 3))])
    mets.to_csv(OUT / "validation_metrics.tsv", sep="\t", index=False)
    print()
    print("External-validation metrics")
    print(mets.T.to_string(header=False))

    # --- 3. applicability-domain check ------------------------------------------------
    rows = []
    for fam, g in cal.groupby("family"):
        g = g.reset_index(drop=True)
        x, y = g.dE_kcal.values, g.log10_kGSH.values
        n = len(x); h = 1 / n + (x - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2)
        err = np.abs(y - loo_predictions(x, y))
        cut = np.median(h)
        for label, mask in [("leverage at or below the family median", h <= cut),
                            ("leverage above the family median", h > cut)]:
            rows.append(dict(family=fam, stratum=label, n=int(mask.sum()),
                             median_leverage=round(float(np.median(h[mask])), 3),
                             mean_absolute_loo_error=round(float(err[mask].mean()), 3),
                             max_absolute_loo_error=round(float(err[mask].max()), 3)))
    adp = pd.DataFrame(rows)
    adp.to_csv(OUT / "applicability_performance.tsv", sep="\t", index=False)

    # --- 4. interval estimates for the fit statistics ---------------------------------
    rows = []
    for fam, g in cal.groupby("family"):
        g = g.reset_index(drop=True)
        x, y, grp = g.dE_kcal.values, g.log10_kGSH.values, g.fragment_smiles.values
        fit = linregress(x, y)
        rows.append(dict(family=fam, n=len(g), fragments=int(pd.Series(grp).nunique()),
                         r2=round(fit.rvalue ** 2, 3),
                         fragment_group_q2=round(q2(y, group_predictions(x, y, grp)), 3),
                         **cluster_bootstrap(x, y, grp)))
    ints = pd.DataFrame(rows)
    ints.to_csv(OUT / "calibration_intervals.tsv", sep="\t", index=False)
    print()
    print("Cluster-bootstrap intervals over fragment groups")
    print(ints.to_string(index=False))
    print()
    print("Held-out error by leverage stratum")
    print(adp.to_string(index=False))


if __name__ == "__main__":
    main()
