# External GSH reactivity set - held-out alpha,beta-unsaturated esters

Purpose: literature-measured GSH second-order rate constants (kGSH; log10 kGSH,
kGSH in L mol-1 min-1) for alpha,beta-unsaturated carbonyls that are HELD OUT of the
36-compound calibration set, for use as an external test set.

Assay (all values, same chemoassay): phosphate buffer pH 7.4, 25 C, 0.14 mM GSH,
DTNB detection at 412 nm; apparent second-order rate on total GSH. No cross-assay
conversion was applied. log10 kGSH is dimensionless with kGSH in L mol-1 min-1.

## Primary target paper (confirmed identity)

Bohme, A.; Thaens, D.; Schramm, F.; Paschke, A.; Schuurmann, G.
"Thiol Reactivity and Its Impact on the Ciliate Toxicity of alpha,beta-Unsaturated
Aldehydes, Ketones, and Esters." Chem. Res. Toxicol. 2010, 23 (12), 1905-1912.
DOI 10.1021/tx100226n. PMID 20923215.
- Title/authors/volume/issue/pages verified via Crossref (api.crossref.org/works/10.1021/tx100226n)
  and PubMed/EuropePMC (PMID 20923215).

### Access attempted for the PRIMARY full text / SI
- https://pubs.acs.org/doi/10.1021/tx100226n  -> HTTP 403 (paywall).
- ACS Supporting Information https://pubs.acs.org/doi/suppl/10.1021/tx100226n/suppl_file/tx100226n_si_001.pdf
  -> HTTP 403 (returns a Cloudflare HTML challenge, not the PDF).
- Unpaywall (api.unpaywall.org/v2/10.1021/tx100226n): is_oa=false, oa_status=closed,
  best_oa_location=null, oa_locations=[].
- Semantic Scholar Graph API (DOI:10.1021/tx100226n): isOpenAccess=false, openAccessPdf=None.
- PubMed / EuropePMC: no PMCID, no open full text.
Conclusion: no open-access copy of the 2010 primary paper exists; its Table could not be
opened directly. Values below were therefore taken from a citable secondary compilation
that faithfully reprints the Bohme data (see cross-validation), with primary attribution
to Bohme 2010.

## Secondary source actually used for the numeric log kGSH values

Mulliner, D. "Quantenchemische Modellierung der Thiol-Addition an Michael-Akzeptoren zur
quantitativen Vorhersage ihrer elektrophilen Reaktivitaet und aquatischen Toxizitaet."
Dissertation, TU Bergakademie Freiberg, degree conferred 31 Jan 2014 (Schuurmann group).
- Table 3.1 ("Alle 74 Verbindungen und Informationen zu deren Toxizitaet, Hydrophobie und
  Reaktivitaet") lists, per compound: CAS, log EC50, log Kow, and log kGSH (L mol-1 min-1).
- Table 3.1 footnote: "Alle Daten zur Toxizitaet ... und Reaktivitaet[13,40,41] wurden der
  Literatur entnommen." The three reactivity references are, per the thesis bibliography:
    [40] Bohme, Thaens, Paschke, Schuurmann 2009, Chem. Res. Toxicol. 22:742-750 (= tx800492x)
    [41] Bohme, Thaens, Schramm, Paschke, Schuurmann 2010, Chem. Res. Toxicol. 23:1905-1912 (= tx100226n)  <-- the target paper
    [13] Wondrousch, Bohme, Thaens, Ost, Schuurmann 2010, "Local Electrophilicity" (QM paper reusing Bohme kGSH)
- Access route: PDF (mulliner_thesis.pdf) already present in the project scratchpad from
  earlier work in this build; text extracted with pymupdf (PyMuPDF 1.28.2), pages 43-46.

### Why the numbers are trustworthy (cross-validation)
Every log kGSH in Mulliner Table 3.1 that overlaps the project's 36-compound calibration
set matches the calibration value (calibration_primary_source_manifest.tsv) to the printed
2 decimal places, e.g.: methyl acrylate 1.06 (CAL15 1.057), ethyl acrylate 1.03 (CAL16),
n-propyl acrylate 1.01 (CAL17), n-butyl acrylate 0.93 (CAL18), methyl crotonate -0.79
(CAL19), ethyl crotonate -0.79 (CAL20), methyl methacrylate -1.14 (CAL21), ethyl
methacrylate -1.24 (CAL22), methyl tiglate -2.15 (CAL23); ketones 1-penten-3-one 3.10
(CAL01), 1-octen-3-one 3.03 (CAL03), etc. The thesis is thus a faithful reprint of the
Bohme data.

## The five held-out esters (exact text taken)

The 2010 paper's ester set was independently enumerated from the paper's OWN dataset
registry: QsarDB entry for tx100226n (handle 10967/15, QDB DOI 10.15152/QDB.15,
file tx100226n.qdb.zip -> compounds/compounds.xml). That registry lists 16 esters
(compound IDs 31-46). Nine alkene-ester rows are the Bohme-2009 esters already in the
36-compound calibration; two are propiolates (excluded, alkyne conjugated to C=O). The
remaining FIVE alkene esters are unique to the 2010 paper and are the held-out set below.
(These five are absent from Bohme 2009 Table 2 per the project's kinetic_source_notes.md;
they are present in the tx100226n registry -> primary source = Bohme 2010, ref [41].)
Note: in propargyl/allyl acrylate and propargyl methacrylate the C=C / C#C in the ALCOHOL
(OR) group is NOT conjugated to the carbonyl; the reacting Michael system is the
acrylate/methacrylate C=C, so these are genuine alpha,beta-unsaturated ESTERS (not propiolates).

Exact rows quoted from Mulliner 2014 Table 3.1 (format: Name | Nr. | CAS | logEC50 | logKow | log kGSH):

1. "Propargylacrylat | C6 | 10477-47-1 | -4.05 | 0.94 | 1.71"
   -> propargyl acrylate; log10 kGSH = 1.71.  PubChem CID 82655, SMILES C=CC(=O)OCC#C
   (name/CAS 10477-47-1 -> CID 82655, IUPAC prop-2-ynyl prop-2-enoate).

2. "Allylacrylat | E7 | 999-55-3 | -3.68 | 1.57 | 1.29"
   -> allyl acrylate; log10 kGSH = 1.29.  PubChem CID 13835, SMILES C=CCOC(=O)C=C
   (CAS 999-55-3 -> CID 13835, IUPAC prop-2-enyl prop-2-enoate).

3. "Propargylmethacrylat | B7 | 13861-22-8 | -2.63 | 1.49 | -0.66"
   -> propargyl methacrylate; log10 kGSH = -0.66.  PubChem CID 83778,
   SMILES CC(=C)C(=O)OCC#C (CAS 13861-22-8 -> CID 83778, IUPAC prop-2-ynyl 2-methylprop-2-enoate).

4. "tert-Butylacrylat | B9 | 1663-39-4 | -3.27 | 2.09 | 0.40"
   -> tert-butyl acrylate; log10 kGSH = 0.40.  PubChem CID 15458,
   SMILES CC(C)(C)OC(=O)C=C (CAS 1663-39-4 -> CID 15458, IUPAC tert-butyl prop-2-enoate).

5. "Methyl-2-octenoat | B12 | 7367-81-9 | -3.76 | 3.10 | -0.11"
   -> methyl trans-2-octenoate; log10 kGSH = -0.11.  PubChem CID 5364532,
   SMILES CCCCC/C=C/C(=O)OC (CAS 7367-81-9 -> CID 5364532, IUPAC methyl (E)-oct-2-enoate).
   (CAS 7367-81-9 is the E/trans isomer; the QDB registry and thesis both use this CAS.)

All five PubChem lookups done via PUG-REST name/<CAS>/property/IUPACName,SMILES/JSON
(pubchem.ncbi.nlm.nih.gov), 2026-09-17.

### kGSH (linear, L mol-1 min-1) column
Left BLANK in the TSV: the linear kGSH values are printed only in the 2010 primary Table
(not openable) and in the Bohme 2016 compilation; the secondary source used here prints
log10 kGSH only. For convenience, kGSH = 10^(log10 kGSH) gives approx: propargyl acrylate
51, allyl acrylate 19.5, tert-butyl acrylate 2.5, propargyl methacrylate 0.22,
methyl trans-2-octenoate 0.78 (these are COMPUTED back-transforms, not source-printed).

### Confidence rationale (all five = high)
- Compound identity: independently confirmed in the 2010 paper's own QsarDB registry.
- log kGSH value: exact 2-dp figure from a citable dissertation whose overlapping values
  reproduce the calibration set to 2 dp (airtight reprint).
- Residual caveat (documented, not a guess): the thesis footnote attributes reactivity to
  refs 13/40/41 collectively, not per row; assignment of these five specifically to ref [41]
  (Bohme 2010) rests on (a) their absence from Bohme 2009 Table 2 and (b) their presence in
  the tx100226n registry. This is a strong, documented inference rather than a per-row
  citation in the thesis. Values are printed to 2 dp only.

## Ketones / enones in the 2010 paper - none held out

The 2010 paper's ketone set (QDB registry IDs 16-30) is: 1-penten-3-one, 1-hexen-3-one,
1-octen-3-one, 3-hexyn-2-one, 3-penten-2-one, 2-octen-4-one, 2-cyclopenten-1-one,
4-hexen-3-one, 3-hepten-2-one, 3-octen-2-one, 3-nonen-2-one, 3-methyl-3-penten-2-one,
4-methyl-3-penten-2-one, 2-methyl-2-cyclopenten-1-one, 3-methyl-2-cyclopenten-1-one.
- 3-hexyn-2-one is an ynone (C#C conjugated to C=O) -> excluded.
- The remaining 14 enones are ALL already in the 36-compound calibration set (CAL01-CAL14,
  sourced to Bohme 2009 tx800492x Table 2, and reprinted identically in Mulliner Table 3.1).
Therefore NO alpha,beta-unsaturated ketone from the 2010 paper is held out; there are no
external enone test points to add from this paper.

## Excluded triple-bond compounds noted (per request, not added)
From the 2010 paper's registry and the 2009 paper: methyl propiolate, ethyl propiolate
(propiolate esters), 3-hexyn-2-one (ynone), 2-octynal, phenyl propiolaldehyde (ynals).
Excluded because the reacting conjugated system contains a C#C triple bond.

## Other compilations checked (not needed as source)
- Schwobel et al., "Measurement and Estimation of Electrophilic Reactivity for Predictive
  Toxicology," Chem. Rev. 2011, 111, 2562-2596, DOI 10.1021/cr100098n (reprints Bohme
  values with citation) - ACS SI download returned a Cloudflare HTML page, not usable.
- Townsend, Farrar, Grayson, ACS Omega 2022, 7, DOI 10.1021/acsomega.2c03739 (PMC9352231):
  open, but its log(kGSH) set = only 23 MAs (nine esters = the 2009 nine, seven aldehydes,
  seven ketones); does NOT contain any of the five held-out esters. Attributes its Bohme
  ester data to the Bohme 2016 compilation (Chem. Res. Toxicol. 2016, 29, 952-962).
- QsarDB tx100226n archive (10.15152/QDB.15) stores only logIGC50 toxicity + logKow, not
  kGSH; used here solely to enumerate the exact 2010 compound list.
