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

## Supplementary Table 1: datasets, search settings and SDRFs

One row per dataset; when a dataset was run with several DIA-NN versions, all used the same SDRF, FASTA
and peptide length. Values from `run_metadata/runs.tsv`; full parameters per run in Supplementary Table S6.

| Accession | DIA-NN version | FASTA database | Peptide length | SDRF |
|---|---|---|---|---|
| **ProteoBench benchmarks** | | | | |
| PXD049412 (Module 9) | 1.8.1, 2.5.1, 2.5.1 (Enterprise) | ProteoBenchFASTA_DDAQuantification_noecoli.fasta | 7-30 | [PXD049412.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/ProteoBench/PXD049412.sdrf.tsv) |
| PXD062685 (Module 5) | 1.8.1, 2.5.1, 2.5.1 (Enterprise) | ProteoBenchFASTA_MixedSpecies_HYE.fasta | 7-30 | [PXD062685.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/ProteoBench/PXD062685.sdrf.tsv) |
| PXD070049 (Module 10) | 1.8.1, 2.5.1, 2.5.1 (Enterprise) | ProteoBenchFASTA_MixedSpecies_HYE.fasta | 7-30 | [PXD070049.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/ProteoBench/PXD070049.sdrf.tsv) |
| Module 7 | 1.8.1, 2.5.1, 2.5.1 (Enterprise) | ProteoBenchFASTA_MixedSpecies_HYE.fasta | 7-30 | [ProteoBench_Module_7.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/ProteoBench/ProteoBench_Module_7.sdrf.tsv) |
| **Bulk cell-line panel** | | | | |
| PXD003539 | 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-30 | [PXD003539.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD003539/PXD003539.sdrf.tsv) |
| PXD030304 | 2.5.1 | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD030304.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD030304/PXD030304.sdrf.tsv) |
| PXD004701 | 2.5.1 | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD004701.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD004701/PXD004701.sdrf.tsv) |
| PXD041421 | 2.5.0 | Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta | 7-40 | [PXD041421.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD041421/PXD041421.sdrf.tsv) |
| PXD017199 | 2.5.0 | Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta | 7-40 | [PXD017199.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD017199/PXD017199.sdrf.tsv) |
| **Single-cell proteomics** | | | | |
| PXD046357 | 1.8.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD046357.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD046357/PXD046357.sdrf.tsv) |
| PXD049412 | 1.8.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD049412.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD049412/PXD049412.sdrf.tsv) |
| MSV000093870 | 2.5.1 | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [MSV000093870-plexDIA.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/MSV000093870/MSV000093870-plexDIA.sdrf.tsv) |
| **Phosphoproteomics** | | | | |
| PXD034128 (biological study) | 2.5.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD034128-phospho-biological-study.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD034128/PXD034128-phospho-biological-study.sdrf.tsv) |
| PXD034128 (high-speed DIA) | 2.5.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD034128-phospho-highspeed-DIA.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD034128/PXD034128-phospho-highspeed-DIA.sdrf.tsv) |
| PXD034623 | 2.5.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD034623-M2-Galectin1-phospho-DIA.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD034623/PXD034623-M2-Galectin1-phospho-DIA.sdrf.tsv) |
| PXD049692 | 2.5.1, 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-40 | [PXD049692-NK-fibrin-IL15-phospho-diaPASEF.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD049692/PXD049692-NK-fibrin-IL15-phospho-diaPASEF.sdrf.tsv) |
| **Spatial proteomics** | | | | |
| PXD064049 | 2.5.1 (Enterprise) | Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta | 7-30 | [PXD064049-CHP212-MYCN-DVP-diaPASEF-carbamidomethyl.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD064049/PXD064049-CHP212-MYCN-DVP-diaPASEF-carbamidomethyl.sdrf.tsv) |
| **Scalability (Fig. 2a node-count sweep)** | | | | |
| PXD071075 | 2.5.0 | UP000005640_9606.fasta | 6-40 (pipeline default) | [PXD071075.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/PXD071075/PXD071075.sdrf.tsv) |
| **Standalone DIA-NN vs quantmsdiann (Supplementary Note 8)** | | | | |
| PXD049412 | 2.5.1 (standalone and quantmsdiann) | ProteoBenchFASTA_DDAQuantification_noecoli.fasta | 7-30 | [PXD049412_local.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/Paired_standalone_vs_quantmsdiann/PXD049412_local.sdrf.tsv) |
| PXD070049 | 2.5.1 (standalone and quantmsdiann) | ProteoBenchFASTA_MixedSpecies_HYE.fasta | 7-30 | [PXD070049_local.sdrf.tsv](https://github.com/ypriverol/quantmsdiann-paper/blob/main/sdrf/Paired_standalone_vs_quantmsdiann/PXD070049_local.sdrf.tsv) |

Sources: the `sdrf/` folder each run published on the PRIDE FTP; for the
phosphoproteomics runs and the PXD071075 sweep, whose FTP folders have no
`sdrf/`, the input SDRF named in the run's SLURM log on the EMBL-EBI cluster.

## Known errors in these files (kept as used)

Four instrument accessions are wrong (checked against PSI-MS in OLS). The SDRFs
uploaded to PRIDE are these files unchanged; only the instrument list in the PRIDE
submission summary uses the correct accessions (`scripts/build_submission_px.py`):

| File(s) | Written | Correct |
|---|---|---|
| ProteoBench/PXD049412, ProteoBench/ProteoBench_Module_7 | Orbitrap Astral, MS:1003356 (= Orbitrap Ascend) | MS:1003378 |
| ProteoBench/PXD070049 | ZenoTOF 8600, MS:1003387 (= MetaboAnnotation) | MS:1003443 |
| ProteoBench/PXD062685 | timsTOF SCP, MS:1003005 (= timsTOF Pro) | MS:1003231 |
| PXD017199 | Q Exactive Plus, MS:1001911 (= Q Exactive) | MS:1002634 |

The ProteoBench SDRFs list only *Homo sapiens* as organism, although modules
5, 7 and 10 are human/yeast/E. coli mixtures.
