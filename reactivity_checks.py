"""Two robustness checks on the reactivity model, reproduced from the deposited
fragment addition energies.

1. Held-out ester test. Five alpha,beta-unsaturated esters from the same chemoassay,
   excluded from the 36-compound calibration, are predicted from the ester regression
   and compared with their measured glutathione rate constants.
2. Thiol-surrogate robustness. Methanethiolate and ethanethiolate addition energies
   over a family-spanning subset are compared to confirm that the minimal-thiol
   surrogate preserves the reactivity ordering the model uses.

Inputs (data/):
    calibration_complete.tsv, external_reactivity_predictions.tsv, thiol_surrogate_check.tsv
Output (out/):
    external_reactivity_predictions.tsv (predictions reproduced from the deposited energies)

Run:  python reactivity_checks.py
Requires: numpy, pandas, scipy
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import linregress, pearsonr, spearmanr

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

    surr = pd.read_csv(DATA / "thiol_surrogate_check.tsv", sep="\t")
    pr = pearsonr(surr.dE_methyl, surr.dE_ethyl)[0]; sp = spearmanr(surr.dE_methyl, surr.dE_ethyl)[0]
    offset = (surr.dE_ethyl - surr.dE_methyl).mean()
    print(f"thiol surrogate: n={len(surr)}  Pearson={pr:.2f}  Spearman={sp:.2f}  mean offset={offset:.2f} kcal/mol")


if __name__ == "__main__":
    main()
