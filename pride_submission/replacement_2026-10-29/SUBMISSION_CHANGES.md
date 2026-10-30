# PXD085732: file changes for the revision

Staging on Codon: `/hps/nobackup/juan/pride/reanalysis/submission-quantmsdiann-r1/` (the current
`submission-quantmsdiann/` stays untouched). Every unchanged file keeps its name and content.

## Replace (same file name, new content)

| File | New content | Status |
|---|---|---|
| `PXD017199_atlas.zip` | DIA-NN 2.5.1 run (cell-lines-proteomes, June 2026), contaminants-only FASTA, no entrapment | to stage |
| `PXD017199_atlas_PXD017199.sdrf.tsv` | as-run SDRF of that run | to stage |
| `PXD017199.qpx.zip` | QPX of that run | to stage |
| `PXD041421_atlas.zip` | DIA-NN 2.5.1 run (cell-lines-proteomes, June 2026), contaminants-only FASTA, no entrapment | to stage |
| `PXD041421_atlas_PXD041421.sdrf.tsv` | as-run SDRF of that run | to stage |
| `PXD041421.qpx.zip` | QPX of that run | to stage |
| `Paired-*_v2_5_0*` (6 files in this folder) | paired standalone DIA-NN 2.5.0 vs quantmsdiann v2.3.0 rerun | staged here |

## Remove, then add under a new name

| Remove | Add | Why |
|---|---|---|
| `PXD003539_v2_5_1_enterprise.zip` | `PXD003539_v2_5_1.zip` | NCI-60 rerun with standard DIA-NN 2.5.1, same settings as the other cell lines (done 2026-10-29) |
| `PXD003539_v2_5_1_enterprise_PXD003539.sdrf.tsv` | `PXD003539_v2_5_1_PXD003539.sdrf.tsv` | as above |
| (none) | `PXD003539.qpx.zip` | QPX of the new run (the run produced `qpx/`) |

## Sources of the new files

The harmonised runs are on the PRIDE FTP (2026-10-30), in the layout of the other cell lines:
`quantmsdiann-benchmarks/cell-lines/{PXD003539,PXD017199,PXD041421}/v2_5_1/` (copied from the
Codon run folders; `quant_tables/diann_report.parquet` sizes match the run outputs byte for byte).
The zips above are built from these folders, with the same layout as `PXD004701_v2_5_1.zip`.

## Unchanged

All other files, including `PXD004701_v2_5_1.zip` and `PXD030304_v2_5_1.zip` (already the harmonised
2.5.1 runs), and their SDRFs and QPX zips.

## Also update in the submission form

- Data processing protocol: cell-line cohorts all DIA-NN 2.5.1 (quantmsdiann de98e98, v2.2.0dev),
  contaminants-only FASTA; no entrapment searches.
- `submission.px`: regenerate after staging (file list and checksums).
