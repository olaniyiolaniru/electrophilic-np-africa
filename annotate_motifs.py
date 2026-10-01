"""Reapply the electrophilic-motif library to the standardized structures and
reproduce the motif inventory, source-family prevalence, and property contrasts.

Inputs (data/):
    motif_library.tsv          22 SMARTS definitions with candidate-atom indices
    anpdb_annotations.tsv      11,424 standardized structures (parent_smiles, family)
    internal_challenge.tsv     43-record boundary/challenge set
    exclusion_rules.tsv        declarative chemical exclusions applied with the library
Outputs (out/):
    motif_matches.tsv, motif_prevalence.tsv, family_prevalence.tsv,
    family_prevalence_normalized.tsv, molecular_descriptors.tsv,
    property_effects.tsv, boundary_tests.tsv, challenge_result.tsv

Run:  python annotate_motifs.py
Requires: numpy, pandas, scipy, rdkit
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, QED

DATA = Path("data"); OUT = Path("out"); OUT.mkdir(exist_ok=True)

# accepted-family synonyms for the source-family ranking
FAM_SYN = {"Labiatae": "Lamiaceae", "Umbelliferae": "Apiaceae", "Compositae": "Asteraceae",
           "Asteraceae-Compositae": "Asteraceae", "Leguminosae": "Fabaceae", "Leguminosae/Fabaceae": "Fabaceae",
           "Gramineae": "Poaceae", "Gramineae-Poaceae": "Poaceae", "Cruciferae": "Brassicaceae",
           "Cruciferae-Brassicaceea": "Brassicaceae", "Guttiferae": "Clusiaceae", "Clusiaceae-Guttiferae": "Clusiaceae",
           "Palmae": "Arecaceae"}

BOUNDARY = [("cyclic ketone", "O=C1C=CCCC1", "enone_cyclic", True),
            ("lactone exclusion", "O=C1OC=CC1", "enone_cyclic", False),
            ("conjugated lactone", "O=C1C=CCO1", "unsat_lactone", True),
            ("lactam exclusion", "O=C1C=CCN1", "enone_cyclic", False),
            ("iodoacetamide", "NC(=O)CI", "alpha_halo_carbonyl", True),
            ("fluoroacetamide", "NC(=O)CF", "alpha_halo_carbonyl", False),
            ("naphthoquinone", "O=C1C=CC(=O)c2ccccc12", "para_quinone", True),
            ("saturated dione", "O=C1CCC(=O)CC1", "para_quinone", False),
            ("epoxide", "CC1CO1", "epoxide", True),
            ("nitrile", "CC#N", "nitrile", True)]


def load_library():
    lib = pd.read_csv(DATA / "motif_library.tsv", sep="\t")
    return [(r, Chem.MolFromSmarts(r["smarts"])) for r in lib.to_dict("records")]


def load_exclusions():
    """Read the declarative chemical exclusions from data/exclusion_rules.tsv.

    The table is the single definition of each rule: its SMARTS and the motifs it applies to
    are taken from there, so the shipped rule and the code that applies it cannot diverge.
    EX01 removes the 3-hydroxy-4-oxo-pyranone (flavonol) core, an aromatic-stabilized
    vinylogous ester rather than a ketone enone, whose perception as an aliphatic enone is
    sensitive to the RDKit aromaticity model.
    """
    out = []
    for r in pd.read_csv(DATA / "exclusion_rules.tsv", sep="\t").to_dict("records"):
        patt = Chem.MolFromSmarts(r["smarts"])
        assert patt is not None, f"unparsable exclusion SMARTS in rule {r['rule_id']}"
        out.append((r["rule_id"], patt, {s.strip() for s in str(r["applies_to_motif_ids"]).split(";")}))
    assert {rid for rid, _p, _m in out} == {"EX01"}, "exclusion_rules.tsv does not match the deposited rule set"
    return out


def match(mol, lib, excl):
    masked = [(motifs, {a for h in mol.GetSubstructMatches(patt) for a in h}) for _rid, patt, motifs in excl]
    out = []
    for r, patt in lib:
        for hit in mol.GetSubstructMatches(patt, uniquify=True):
            if any(r["motif_id"] in motifs and set(hit) & atoms for motifs, atoms in masked):
                continue
            centres = sorted({hit[int(i)] for i in str(r["candidate_query_indices"]).split(";")})
            out.append((r, hit, centres))
    return out


def main():
    lib = load_library(); excl = load_exclusions()
    for name, smi, motif, expected in BOUNDARY:
        obs = motif in {r["motif_id"] for r, h, c in match(Chem.MolFromSmiles(smi), lib, excl)}
        assert obs == expected, name
    pd.DataFrame([dict(case=n, smiles=s, motif=m, expected=e,
                       observed=(m in {r["motif_id"] for r, h, c in match(Chem.MolFromSmiles(s), lib, excl)}))
                 for n, s, m, e in BOUNDARY]).to_csv(OUT / "boundary_tests.tsv", sep="\t", index=False)

    ann = pd.read_csv(DATA / "anpdb_annotations.tsv", sep="\t")
    matches, desc, rows_out = [], [], []
    for row in ann.to_dict("records"):
        m = Chem.MolFromSmiles(row["parent_smiles"]); m = Chem.MolFromSmiles(Chem.MolToSmiles(m))
        hits = match(m, lib, excl)
        row["is_electrophilic"] = int(bool(hits))
        row["n_motif_matches"] = len(hits)
        row["motif_ids"] = ";".join(sorted({r["motif_id"] for r, h, c in hits}))
        row["mechanism_groups"] = ";".join(sorted({r["mechanism_group"] for r, h, c in hits}))
        rows_out.append(row)
        for j, (r, h, c) in enumerate(hits, 1):
            matches.append(dict(molecule_id=row["molecule_id"], motif_match_id=f"{row['molecule_id']}::M{j:02d}",
                                motif_id=r["motif_id"], mechanism_group=r["mechanism_group"],
                                matched_atom_indices=";".join(map(str, h)),
                                candidate_reactive_atom_indices=";".join(map(str, c)),
                                parent_smiles=Chem.MolToSmiles(m)))
        desc.append(dict(molecule_id=row["molecule_id"], is_electrophilic=row["is_electrophilic"],
                         MW=Descriptors.MolWt(m), logP=Descriptors.MolLogP(m), TPSA=Descriptors.TPSA(m),
                         Fsp3=rdMolDescriptors.CalcFractionCSP3(m),
                         n_stereocentres=len(Chem.FindMolChiralCenters(m, includeUnassigned=True)),
                         n_rings=rdMolDescriptors.CalcNumRings(m),
                         n_arom_rings=rdMolDescriptors.CalcNumAromaticRings(m), HBD=Descriptors.NumHDonors(m),
                         HBA=Descriptors.NumHAcceptors(m), QED=QED.qed(m), BertzCT=Descriptors.BertzCT(m)))
    annotated = pd.DataFrame(rows_out); st = pd.DataFrame(matches); ds = pd.DataFrame(desc)
    st.to_csv(OUT / "motif_matches.tsv", sep="\t", index=False)
    ds.to_csv(OUT / "molecular_descriptors.tsv", sep="\t", index=False)

    counts = st.groupby("motif_id").molecule_id.nunique().reset_index(name="n_unique_compounds")
    counts.to_csv(OUT / "motif_prevalence.tsv", sep="\t", index=False)
    fam = annotated[annotated.molecule_id != "Mol_00042"].groupby("family").is_electrophilic.agg(["size", "sum"]).reset_index()
    fam.columns = ["family", "n_unique", "n_electrophilic"]; fam["pct_electrophilic"] = 100 * fam.n_electrophilic / fam.n_unique
    fam.to_csv(OUT / "family_prevalence.tsv", sep="\t", index=False)
    sg = fam[~fam.family.str.contains(r"\|")].copy(); sg["family"] = sg.family.map(lambda x: FAM_SYN.get(x, x))
    norm = sg.groupby("family", as_index=False).agg(n_unique=("n_unique", "sum"), n_electrophilic=("n_electrophilic", "sum"))
    norm["pct_electrophilic"] = 100 * norm.n_electrophilic / norm.n_unique
    norm.sort_values("pct_electrophilic", ascending=False).to_csv(OUT / "family_prevalence_normalized.tsv", sep="\t", index=False)

    eff = []
    for c in ds.columns[2:]:
        a = ds.loc[ds.is_electrophilic == 1, c]; b = ds.loc[ds.is_electrophilic == 0, c]
        u = mannwhitneyu(a, b, alternative="two-sided")
        eff.append(dict(descriptor=c, electrophilic_median=a.median(), nonelectrophilic_median=b.median(),
                        p_value=u.pvalue, cliff_delta=2 * u.statistic / (len(a) * len(b)) - 1))
    pd.DataFrame(eff).to_csv(OUT / "property_effects.tsv", sep="\t", index=False)

    ch = pd.read_csv(DATA / "internal_challenge.tsv", sep="\t")
    ch["flagged"] = [int(bool(match(Chem.MolFromSmiles(s), lib, excl))) for s in ch.smiles]
    ch.to_csv(OUT / "challenge_result.tsv", sep="\t", index=False)
    tp = int(((ch.is_covalent == 1) & (ch.flagged == 1)).sum()); tn = int(((ch.is_covalent == 0) & (ch.flagged == 0)).sum())
    fp = int(((ch.is_covalent == 0) & (ch.flagged == 1)).sum()); fn = int(((ch.is_covalent == 1) & (ch.flagged == 0)).sum())

    michael = int(annotated.mechanism_groups.str.contains("michael_acceptor").sum())
    cand = len({(r["molecule_id"], int(i)) for r in matches for i in r["candidate_reactive_atom_indices"].split(";")})
    ast = norm[norm.family == "Asteraceae"].iloc[0]
    print(f"structures {len(annotated)}  electrophilic {int(annotated.is_electrophilic.sum())}  "
          f"matches {len(st)}  candidate atoms {cand}  Michael {michael}")
    print(f"challenge TP {tp} TN {tn} FP {fp} FN {fn} | Asteraceae {int(ast.n_electrophilic)}/{int(ast.n_unique)}")

    # Machine-checkable inventory assertions: fail loudly if the environment (e.g. a change
    # in RDKit aromaticity handling) alters the deposited motif inventory.
    per_motif = counts.set_index("motif_id").n_unique_compounds.to_dict()
    expected = {"structures": (len(annotated), 11424), "electrophilic": (int(annotated.is_electrophilic.sum()), 3414),
                "matches": (len(st), 5671), "candidate_atoms": (cand, 6222), "michael": (michael, 2944),
                "acrylate_ester": (per_motif.get("acrylate_ester"), 1545), "enone_cyclic": (per_motif.get("enone_cyclic"), 812),
                "chalcone": (per_motif.get("chalcone"), 79), "challenge": ((tp, tn, fp, fn), (22, 19, 2, 0))}
    for name, (got, exp) in expected.items():
        assert got == exp, f"inventory drift in {name}: got {got}, expected {exp} (check RDKit version, pinned to 2026.03.6)"
    print("inventory assertions passed")


if __name__ == "__main__":
    main()
