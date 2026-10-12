# SDRF snapshot

The exact SDRF files used by every quantmsdiann run in the manuscript, copied
verbatim on 2026-10-12. They are kept here because the curated copies in
[bigbio/sdrf-annotated-datasets](https://github.com/bigbio/sdrf-annotated-datasets)
keep evolving; these are the versions the results were produced from.

| Folder | Runs |
|---|---|
| `ProteoBench/` | ProteoBench modules 5 (PXD062685), 7, 9 (PXD049412) and 10 (PXD070049); DIA-NN 1.8.1, 2.5.1, 2.5.1-enterprise |
| `PXD046357/`, `PXD049412/`, `MSV000093870/` | Single-cell cohorts (PXD049412 here is the A549/H460 cohort, not Module 9) |
| `PXD003539/`, `PXD004701/`, `PXD030304/`, `PXD017199/`, `PXD041421/` | Bulk cell-line cohorts |
| `PXD034128/`, `PXD034623/`, `PXD049692/` | Phosphoproteomics |
| `PXD064049/` | Spatial Deep Visual Proteomics |
| `PXD071075/` | Node-count sweep (Fig. 2a) |
| `Paired_standalone_vs_quantmsdiann/` | Paired standalone DIA-NN vs quantmsdiann runs (Supplementary Note 8) |

Which run used which file: `run_metadata/runs.tsv` (column `sdrf`).

Sources: the `sdrf/` folder each run published on the PRIDE FTP; for the
phosphoproteomics runs and the PXD071075 sweep, whose FTP folders have no
`sdrf/`, the input SDRF named in the run's SLURM log on the EMBL-EBI cluster.

## Known errors in these files (kept as used)

Four instrument accessions are wrong (checked against PSI-MS in OLS). The copies
uploaded to PRIDE are corrected (`scripts/build_submission_px.py`):

| File(s) | Written | Correct |
|---|---|---|
| ProteoBench/PXD049412, ProteoBench/ProteoBench_Module_7 | Orbitrap Astral, MS:1003356 (= Orbitrap Ascend) | MS:1003378 |
| ProteoBench/PXD070049 | ZenoTOF 8600, MS:1003387 (= MetaboAnnotation) | MS:1003443 |
| ProteoBench/PXD062685 | timsTOF SCP, MS:1003005 (= timsTOF Pro) | MS:1003231 |
| PXD017199 | Q Exactive Plus, MS:1001911 (= Q Exactive) | MS:1002634 |

The ProteoBench SDRFs list only *Homo sapiens* as organism, although modules
5, 7 and 10 are human/yeast/E. coli mixtures.
