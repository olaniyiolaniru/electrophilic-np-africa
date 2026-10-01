# Electrophilic motifs and glutathione reactivity in African natural products

Code and data for an atom-resolved map of electrophilic functionality across the African
Natural Products Database (ANPDB), coupled to a source-calibrated glutathione (GSH)
reactivity model. A 22-pattern screen identifies 3,414 motif-bearing structures among
11,424 records; local-fragment methanethiolate-addition energies are calibrated against 36
literature GSH rate constants; and site-resolved estimates define a focused pool of 519
structures for thiol-reactivity study. Every released estimate carries its prediction
interval, applicability flags, and fragment-fidelity diagnostics, which together define a
205-structure enone applicability domain and a 41-structure high-fidelity tier within it.

## Scripts
Each script reads the tables in `data/` and writes its outputs to `out/`.

- `annotate_motifs.py` - reapply the 22-pattern motif library to the standardized
  structures and reproduce the motif inventory, source-family prevalence, property
  contrasts, and boundary/challenge tests. Chemical exclusions are read from
  `data/exclusion_rules.tsv`, which is the single declaration of each rule.
- `compute_fragment_energies.py` - build the radius-2 fragment and enolate adduct for each
  mapped site and compute the three-start GFN2-xTB/ALPB-water addition energy (needs xTB).
- `model_and_selection.py` - fit the per-family calibrations, cross-validate, and build the
  site-level predictions, selection pool, and applicability and fragment-fidelity diagnostics.
- `benchmark_descriptor.py` - compare the addition energy with a family mean, three-feature
  structural substitution counts, and a 256-bit fragment fingerprint on identical
  fragment-group splits, and measure how finely each descriptor separates the released sites.
- `compare_chemical_space.py` - property-matched comparison with the frozen COCONUT
  reference at four Tanimoto thresholds with paired bootstrap intervals.
- `reactivity_checks.py` - held-out ester test and the methanethiolate/ethanethiolate
  surrogate-robustness check.
- `validate_annotations.py` - validate a stratified annotation sample against independent
  PubChem structure records and primary-literature identifiers.

## Requirements
Python 3.14 with RDKit, pandas, numpy, and scipy. Exact versions used to generate the
deposited tables are pinned in `requirements.txt` and `environment.yml` (RDKit installs most
easily with conda). `annotate_motifs.py` asserts the headline inventory counts and fails
loudly if a different environment changes them, because SMARTS matching and aromaticity
perception can shift the inventory between RDKit versions. Energy calculations use xTB 6.7.1,
available separately from https://github.com/grimme-lab/xtb; set `XTB_EXE` to its executable.

## Reproduction
The modeling results reproduce from the deposited energies without xTB:

    python model_and_selection.py     # calibration, cross-validation, pool (598 sites / 519 structures)
    python annotate_motifs.py         # 3,414 motif-bearing structures; 5,671 matches (asserted)
    python benchmark_descriptor.py    # descriptor comparison and resolution on fragment-group splits
    python compare_chemical_space.py  # property-matched COCONUT comparison
    python reactivity_checks.py       # held-out ester RMSD and thiol-surrogate correlation

Recomputing the addition energies from structures requires xTB and is driven by
`compute_fragment_energies.py`. The optimized geometries and per-start xTB calculation
records are large and are supplied as a separate archive accompanying the article.

## Data sources
ANPDB: https://african-compounds.org (SHA-256 147e46ba64deacdda12d0fc4ebe8cb99c27e5527339225e47b79f46ff065c28c). COCONUT September 2024: https://coconut.naturalproducts.net (SHA-256 1900a3824324dec786c4b41b040cc39c02dcb985987fb8eb238c4857a7c1e49a). All derived tables in `data/` are traceable to these sources by identifier and
checksum. Compound names, species strings and literature fields are carried exactly as
retrieved, so a few of them contain source typography; where a source string needed an
accepted reading for presentation, both forms are kept side by side.

## Citation
If you use this resource, please cite the accompanying article in the Journal of Natural
Products and this repository (see `CITATION.cff`).

## License
Code is released under the MIT License (`LICENSE`).
