"""Two robustness checks on the reactivity model, reproduced from the deposited
fragment addition energies.

1. Held-out ester test. Five alpha,beta-unsaturated esters from the same chemoassay,
   excluded from the 36-compound calibration, are predicted from the ester regression
   and compared with their measured glutathione rate constants.
2. Thiol-surrogate sensitivity. Methanethiolate and ethanethiolate addition energies over
   a family-spanning subset are compared, and the offset between them is expressed on the
   predicted log scale through the enone slope, giving the size of the descriptor shift
   that the choice of minimal thiol carries.

Inputs (data/):
    calibration_complete.tsv, external_reactivity_predictions.tsv, thiol_surrogate_check.tsv
Output (out/):
    external_reactivity_predictions.tsv (predictions reproduced from the deposited energies)

Run:  python reactivity_checks.py
Requires: numpy, pandas, scipy
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import linregress

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)


def main():
    cal = pd.read_csv(DATA / "calibration_complete.tsv", sep="\t")
    e = cal[cal.family == "ester"]
    fit = linregress(e.dE_kcal.values, e.log10_kGSH.values)

    held = pd.read_csv(DATA / "external_reactivity_predictions.tsv", sep="\t")
    held["predicted_log10_kGSH"] = fit.intercept + fit.slope * held.dE_kcal
    rmsd = np.sqrt(np.mean((held.measured_log10_kGSH - held.predicted_log10_kGSH) ** 2))
    held.to_csv(OUT / "external_reactivity_predictions.tsv", sep="\t", index=False)
    print(f"held-out ester test: n={len(held)}  RMSD={rmsd:.2f} log units  "
          f"(all within energy range: {bool(held.within_energy_range.all())})")

    # Surrogate sensitivity. The 12 rows resolve to 6 distinct energy pairs, because compounds
    # sharing a local fragment share a calculation, so the comparison is reported as the size and
    # the spread of the offset rather than as a correlation over repeated values.
    surr = pd.read_csv(DATA / "thiol_surrogate_check.tsv", sep="\t")
    surr["offset"] = surr.dE_ethyl - surr.dE_methyl
    pairs = surr.drop_duplicates(subset=["dE_methyl", "dE_ethyl"])
    by_family = ", ".join(f"{k} {v}" for k, v in pairs.family.value_counts().sort_index().items())
    shift = pairs.offset.abs()
    print(f"thiol surrogate: {len(surr)} compounds resolving to {len(pairs)} distinct energy pairs ({by_family})")
    print(f"  ethanethiolate is more exothermic by {shift.min():.2f} to {shift.max():.2f} kcal/mol")

    en = cal[cal.family == "enone"]; x, y = en.dE_kcal.values, en.log10_kGSH.values
    slope = abs(linregress(x, y).slope)
    loo = [linregress(np.delete(x, i), np.delete(y, i)) for i in range(len(x))]
    loo_rmse = np.sqrt(np.mean([(y[i] - (f.intercept + f.slope * x[i])) ** 2 for i, f in enumerate(loo)]))
    spread = shift[pairs.family == "enone"].max() - shift[pairs.family == "enone"].min()
    print(f"  enone offsets span {spread:.2f} kcal/mol, which the enone slope of {slope:.4f} log units "
          f"per kcal/mol maps to {spread * slope:.2f} log units")
    print(f"  for comparison, the enone leave-one-out error is {loo_rmse:.3f} log units")


if __name__ == "__main__":
    main()
