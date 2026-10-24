# -*- coding: utf-8 -*-
"""Write the released prediction set as a table a reader can use without running anything.

The working table that the modelling produces is keyed by record and site identifiers, which
is what the pipeline needs but not what a reader needs. This script joins it to the provenance
record so that every released site carries the compound name, the source organism, its family
and kingdom beside the predicted rate constant, its interval and the diagnostics that qualify
it, and sorts it so the fastest predicted sites come first within each family.

Nothing is recomputed here. Every number is copied from the released table, so this is a
presentation of the deposited result rather than a second derivation of it.

Inputs:   data/selection_pool.tsv, data/anpdb_record_provenance.tsv
Output:   out/released_sites_annotated.tsv

Run:  python released_sites_table.py
Requires: pandas
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
DATA, OUT = HERE / "data", HERE / "out"
OUT.mkdir(exist_ok=True)

COLUMNS = ["compound", "source_species", "genus", "source_family", "kingdom", "family",
           "model_status", "dE_kcal", "pred_log10_kGSH", "prediction_interval_low",
           "prediction_interval_high", "fragment_fidelity_flags", "fidelity_clean",
           "record_id", "site_id", "parent_smiles", "fragment_smiles"]


def main():
    pool = pd.read_csv(DATA / "selection_pool.tsv", sep=chr(9))
    prov = pd.read_csv(DATA / "anpdb_record_provenance.tsv", sep=chr(9))

    keep = prov[["molecule_id", "source_species", "family", "kingdom"]].rename(
        columns={"molecule_id": "record_id", "family": "source_family"})
    out = pool.merge(keep, on="record_id", how="left")
    # the genus is the first token of a source string that names one organism
    single = ~out.source_species.fillna("").str.contains(r"\|")
    out["genus"] = out.source_species.where(single).str.split().str[0]
    out["fidelity_clean"] = out.fragment_fidelity_flags.eq("none").map({True: "yes", False: "no"})
    out = out[COLUMNS].sort_values(["family", "pred_log10_kGSH"],
                                   ascending=[True, False]).reset_index(drop=True)
    out.to_csv(OUT / "released_sites_annotated.tsv", sep=chr(9), index=False)

    named = int(out.source_species.notna().sum())
    print(f"released_sites_annotated.tsv: {len(out)} sites, {out.record_id.nunique()} structures")
    print(f"  source organism recorded for {named} of {len(out)} rows")
    print(f"  source families represented: {out.source_family.dropna().nunique()}")
    print(f"  kingdoms: {out.kingdom.dropna().value_counts().head(5).to_dict()}")
    enone = out[out.family == "enone"]
    print(f"  enone sites with a predicted rate: {len(enone)}, "
          f"{int(enone.fidelity_clean.eq('yes').sum())} of them fidelity-clean")
    print()
    print("Fastest ten predicted enone sites:")
    cols = ["compound", "source_species", "source_family", "pred_log10_kGSH", "fidelity_clean"]
    print(enone.head(10)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
