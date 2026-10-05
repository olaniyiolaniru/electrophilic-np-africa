# Electrophilic motifs and glutathione reactivity in African natural products

Archived release: https://doi.org/10.5281/zenodo.23168309 (always the newest version)

Code and data for an atom-resolved map of electrophilic functionality across the African
Natural Products Database (ANPDB), coupled to a source-calibrated glutathione (GSH)
reactivity model. A 22-pattern screen identifies 3,414 motif-bearing structures among
11,424 records; local-fragment methanethiolate-addition energies are calibrated against 36
literature GSH rate constants; and site-resolved estimates define a focused pool of 519
structures for thiol-reactivity study. Every released estimate carries its prediction
interval, applicability flags, and fragment-fidelity diagnostics, which together define a
205-structure enone applicability domain and a 41-structure high-fidelity tier within it.

## Reproducing everything with one command

    python reproduce_all.py

That runs the full workflow in order and then compares every regenerated table with its
deposited counterpart, reporting the column order, the text fields and the largest numerical
difference for each. It needs no xTB and no network, and it reports the agreement for every table it checks.

## Scripts
Each script reads the tables in `data/` and writes its outputs to `out/`. They can be run
individually in any order; `reproduce_all.py` simply runs them all and checks the result.

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
  fragment-group splits; fit the nested structural, energy and combined models on those same
  splits; and measure how finely each descriptor separates the released sites. Writes
  `descriptor_baseline_comparison.tsv`, `descriptor_nested_comparison.tsv` and
  `descriptor_resolution.tsv`.
- `compare_chemical_space.py` - property-matched comparison with the frozen COCONUT
  reference at four Tanimoto thresholds with paired bootstrap intervals.
- `reactivity_checks.py` - held-out ester test, and the methanethiolate/ethanethiolate
  surrogate sensitivity expressed on the predicted log scale.
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
The modeling results reproduce from the deposited energies without xTB. `model_and_selection.py`
asserts its output schema and rewrites the released tables column for column, so a rerun is
comparable to the deposit field by field, not only in its predictions:

    python model_and_selection.py     # calibration, cross-validation, pool (598 sites / 519 structures)
    python annotate_motifs.py         # 3,414 motif-bearing structures; 5,671 matches (asserted)
    python benchmark_descriptor.py    # descriptor baselines, nested test, resolution counts
    python compare_chemical_space.py  # property-matched COCONUT comparison
    python reactivity_checks.py       # held-out ester RMSD and thiol-surrogate sensitivity

Each script writes into `out/` under the same file name the corresponding table carries in
`data/`, so a fresh run can be compared with the deposited table field by field. On the
versions pinned here every regenerated analysis table matches in shape, column order and text,
with numerical agreement at the limit of floating-point arithmetic.

Recomputing the addition energies from structures requires xTB and is driven by
`compute_fragment_energies.py`. The optimized geometries, program output, optimization
trajectories and per-start energies for all 1,486 calculated species are deposited alongside this
package in the archived record at https://doi.org/10.5281/zenodo.23168310, with a manifest
and per-file checksums. Each site row names
the three calculation keys behind its predicted rate constant, so any estimate in the article
can be followed to the calculations it came from.

## Data sources
ANPDB: https://african-compounds.org (SHA-256 147e46ba64deacdda12d0fc4ebe8cb99c27e5527339225e47b79f46ff065c28c). COCONUT September 2024: https://coconut.naturalproducts.net (SHA-256 1900a3824324dec786c4b41b040cc39c02dcb985987fb8eb238c4857a7c1e49a). All derived tables in `data/` are traceable to these sources by identifier and
checksum. Compound names, species strings and literature fields are carried exactly as
retrieved, so a few of them contain source typography; where a source string needed an
accepted reading for presentation, both forms are kept side by side.

## Citation
Please cite the accompanying article in the Journal of Natural Products together with the
archived record. Cite https://doi.org/10.5281/zenodo.23168310 for the exact version you used,
or https://doi.org/10.5281/zenodo.23168309 to refer to the resource across versions.
`CITATION.cff` carries the machine-readable form.

## Licence and third-party content
The analysis code in this repository is released under the MIT Licence (`LICENSE`).

The tables in `data/` are not original code. `anpdb_annotations.tsv` and
`anpdb_record_provenance.tsv` derive from the African Natural Products Database export named
under Data sources, and `coconut_frozen_reference.tsv` with its matched-pair and
nearest-neighbour tables derives from the COCONUT export named there. Both are redistributed in
derived form so that the reported results can be verified, and both remain subject to the terms
of their originating databases. The kinetic records are transcribed from the two cited primary
papers and remain subject to the publisher's terms. Each source is identified by URL and SHA-256
checksum, so the provenance of any row can be established. Anyone reusing the underlying
collections themselves, rather than checking the results reported here, should obtain them from
the sources named above under those sources' own terms.
