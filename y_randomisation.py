"""Chance-correlation check for the family calibrations by response permutation.

The response is permuted within a family and the model refitted, many times, so that the
fitted and cross-validated performance of the real model can be compared with the
distribution obtained when no true relationship remains. A model that owes its performance
to chance produces permuted statistics close to the real ones.

Both the ordinary fit and the fragment-group cross-validation are permuted, because the
fragment-group split is the stricter test and the one the article reports.

Inputs (data/):  calibration_complete.tsv
Output (out/):   y_randomisation.tsv

Run:  python y_randomisation.py
Requires: numpy, pandas, scipy
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import linregress

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)
N_PERM = 2000
SEED = 42


def group_q2(x, y, groups):
    """Leave-one-fragment-group-out Q squared, as used for the reported models."""
    pred = np.empty(len(y))
    for i in range(len(y)):
        keep = groups != groups[i]
        if keep.sum() < 3:
            pred[i] = y[keep].mean() if keep.any() else y.mean()
            continue
        f = linregress(x[keep], y[keep])
        pred[i] = f.intercept + f.slope * x[i]
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")
    rng = np.random.default_rng(SEED)
    rows = []
    for fam, g in cal.groupby("family"):
        g = g.reset_index(drop=True)
        x = g.dE_kcal.values; y = g.log10_kGSH.values; groups = g.fragment_smiles.values
        real_r2 = linregress(x, y).rvalue ** 2
        real_q2 = group_q2(x, y, groups)
        r2p = np.empty(N_PERM); q2p = np.empty(N_PERM)
        for k in range(N_PERM):
            yp = rng.permutation(y)
            r2p[k] = linregress(x, yp).rvalue ** 2
            q2p[k] = group_q2(x, yp, groups)
        rows.append(dict(
            family=fam, n=len(g), permutations=N_PERM,
            r2=round(real_r2, 3), r2_permuted_mean=round(float(r2p.mean()), 3),
            r2_permuted_p95=round(float(np.percentile(r2p, 95)), 3),
            r2_permuted_max=round(float(r2p.max()), 3),
            r2_p_value=round(float((r2p >= real_r2).mean()), 4),
            fragment_group_q2=round(real_q2, 3),
            q2_permuted_mean=round(float(q2p.mean()), 3),
            q2_permuted_p95=round(float(np.percentile(q2p, 95)), 3),
            q2_permuted_max=round(float(q2p.max()), 3),
            q2_p_value=round(float((q2p >= real_q2).mean()), 4)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "y_randomisation.tsv", sep="\t", index=False)
    print(out[["family", "n", "r2", "r2_permuted_mean", "r2_p_value",
               "fragment_group_q2", "q2_permuted_mean", "q2_p_value"]].to_string(index=False))
    print()
    for r in out.itertuples():
        print(f"{r.family}: real r2 {r.r2:.3f} against permuted mean {r.r2_permuted_mean:.3f} "
              f"(max {r.r2_permuted_max:.3f}), p = {r.r2_p_value:.4f}; "
              f"group Q2 {r.fragment_group_q2:.3f} against permuted mean {r.q2_permuted_mean:.3f}, "
              f"p = {r.q2_p_value:.4f}")


if __name__ == "__main__":
    main()
