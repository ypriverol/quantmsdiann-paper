# Reproducibility and data provenance

Every figure, table and number in the manuscript is derived from public data by a script in
this repository. DIA-NN reports are not committed (tens to hundreds of MB each): each stage
downloads what it needs from the PRIDE FTP into `data/` (git-ignored) and caches it.

## Data sources

| Source | Location | What |
|---|---|---|
| Our reanalyses | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/` | One folder per dataset and DIA-NN version (`v1_8_1`, `v2_5_1`, `v2_5_1_enterprise`, ...): `quant_tables/diann_report.{parquet,tsv}`, matrices, `pipeline_info/`, `sdrf/`. Subtrees `proteobench/`, `single-cell/`, `cell-lines/`, `phospho/`, `spatial/`, `PXD071075_cluster_sizes/`. |
| Cell-line cohorts | `.../quantmsdiann-benchmarks/cell-lines/<PXD>/v2_5_1/` | PXD003539, PXD004701, PXD030304, PXD017199, PXD041421: one harmonised run each (DIA-NN 2.5.1, reviewed human + contaminants FASTA, no entrapment, same parameters). Read through `scripts/cellline_reports.py` (local copy in `data/cell-lines-v2_5_1/<PXD>/` first, never purged by the rebuild). |
| Original deposits | `https://ftp.pride.ebi.ac.uk/pride/data/archive/<YYYY>/<MM>/<PXD>/`, MassIVE, figshare | The originally deposited results compared against in Fig. 4, Fig. S7 and Table S4. |
| Run inputs | [`sdrf/`](sdrf/) (this repository) | SDRF, final Nextflow parameters and FASTA of every run; `sdrf/runs.tsv` indexes all 38 runs. |
| Pipeline and containers | <https://github.com/bigbio/quantmsdiann> (v2.2.0; the cell-line cohorts used commit de98e98, the paired comparison v2.3.0), <https://github.com/bigbio/quantms-containers> | Workflow and DIA-NN container recipes. |

The same results are deposited as a PRIDE reanalysis dataset and on Zenodo (built by the
deposition scripts below; accession and DOI in the paper's Data Availability).

## Counting rule

[methods.md](methods.md) §1, implemented by `count_report` in `scripts/rebuild.py`: global
counts at `Lib.Q.Value` / `Lib.PG.Q.Value` <= 0.01, per-run counts at `Q.Value` /
`PG.Q.Value` <= 0.01, decoys dropped, no contaminant or quantity filter. Counts always come
from the report, never from the `*_matrix.tsv` files.

## Figures and tables

`R` = `python -m scripts.rebuild --only`. Output folder is `analysis/figures/`.

| Item | Output | Command |
|---|---|---|
| Fig. 1 workflow | `manuscript/fig1_workflow.svg` (drawn, no data) | none |
| Fig. 2 runtime scaling + ProteoBench accuracy | `manuscript/fig2_validation.svg` | `R queue_sweep benchmarks proteobench_accuracy performance_trace fig2_validation` |
| Fig. 3 single-cell | `manuscript/fig3_single_cell_combined.svg` | `R single_cell_tables single_cell_combined` |
| Fig. 4 reanalysis vs original | `manuscript/fig_reanalysis_improvement.svg` | `R reanalysis_pg_counts reanalysis_improvement_data reanalysis_improvement` |
| Fig. S1 per-step runtime | `performance/FigS1_runtime_per_step.png` (panels `runtime_per_step.svg`, `mdc_cluster_runtime.svg`) | `R performance_trace mdc_cluster_runtime && python -m scripts.compose_figs1` |
| Fig. S2 memory / CPU per step | `performance/resources_per_step.svg` | `R performance_trace` |
| Fig. S3, Table S3 wall-clock vs number of files | `review/supplementary_revision/FigS7_runtime_vs_files.svg`, `TableS4_runtime_vs_scale.tsv` | see below |
| Fig. S4 peptides per protein (PXD003539) | `PXD003539/supp_peptides_per_protein.svg` | `R pxd003539` |
| Fig. S5 accession overlap (PXD003539) | `PXD003539/supp_venn_protein_accessions.svg` | `R pxd003539` |
| Fig. S6 pan-cohort of cell-line reanalyses | `manuscript/fig4_cellline_atlas.svg` | `R atlas` (needs numpy<2) |
| Fig. S7, Table S4 proteins by peptide evidence | `review/supplementary_revision/FigS8_proteins_by_peptide_evidence.svg`, `TableS5_proteins_by_peptide_evidence.tsv` | see below |
| Fig. S8, Table S6 standalone DIA-NN vs quantmsdiann | `review/supplementary_revision/FigS6_standalone_vs_quantmsdiann.svg`, `TableS3_standalone_vs_quantmsdiann.tsv` | see below |
| Table S1 datasets, search settings, SDRFs | [`sdrf/README.md`](sdrf/README.md) | `python -m scripts.build_run_metadata` (writes `sdrf/runs.tsv`) |
| Table S2 search parameters of every run | `review/supplementary_revision/TableS6_run_parameters.tsv` | `python -m scripts.build_run_metadata` |
| Tables S5 (software) and S7 (feature comparison) | supplementary text | no data |
| Table S8 precursors and peptides per dataset and version | `review/supplementary_revision/TableS8_identifications.tsv` | `python -m scripts.review_table_identifications` |
| Response Table 2, integration of the cell-line cohorts | `review/integration/` | `python -m scripts.review_integration_fdr && python -m scripts.review_integration_figure` |
| Manuscript numbers | `data/paper_numbers.tsv` | `R paper_numbers` |

File names under `review/supplementary_revision/` predate the final supplement numbering; the
table above maps each file to its number in the supplement.

Revision items (Figs. S3, S7, S8; Tables S3, S4, S6, S8), in this order:

```bash
python -m scripts.review_desktop_vs_quantmsdiann   # paired runs -> Table S6 statistics
python -m scripts.review_desktop_overview          # Fig. S8
python -m scripts.review_runtime_vs_scale          # wall-clock vs files / GB (Fig. S3, Table S3)
python -m scripts.review_threshold_regimes         # streams the Fig. 4 reports
python -m scripts.review_peptide_overlap
python -m scripts.review_supplementary             # Figs. S3, S7, S8; Tables S3, S4, S6
python -m scripts.review_table_identifications     # Table S8
```

The paired standalone/quantmsdiann runs (Fig. S8) are read from
`data/proteobench-review-experiments/codon/` (DIA-NN 2.5.0 standalone and quantmsdiann v2.3.0,
same cluster, same predicted library, files in name order); their reports, matrices and launch
scripts (`sdrf/Paired_standalone_vs_quantmsdiann/`) are part of the PRIDE / Zenodo deposit.

## Original (baseline) counts in Fig. 4

Recomputed from the public deposits, not transcribed from the papers:

- NCI-60, Guo 2019 (PXD003539): deposited OpenSWATH `feature_alignment_requant_matrix.tsv`.
- ProCan, Goncalves 2022 (PXD030304): figshare peptide-count matrix (doi:10.6084/m9.figshare.19345397).
- Sun 2023 (PXD004701): deposited OpenSWATH `1R_all_id.txt` and library.
- Spatial DVP (PXD064049): deposited `DIANN_results.zip` (DIA-NN 1.8.1).

## Deposition of the results (PRIDE reanalysis dataset and Zenodo)

1. `python -m scripts.build_pride_raw_list` writes `pride_submission/raw_files.{tsv,txt}`.
2. `python -m scripts.build_submission_manifest` crawls the run folders on the FTP.
3. `python -m scripts.build_submission_px` writes the SDRF copies (identical to `sdrf/`), `raw_files.raw.zip` (the raw list, uploaded as the RAW-type file),
   `submission_files_final.tsv`, `zip_members.tsv` (contents of each `<dataset>_<version>.zip` run
   bundle), `zenodo_files.tsv` and `submission.px`.
4. PRIDE: on the EBI cluster, `python -m scripts.stage_pride_submission --dest <dir>`
   (88 files, ~89 GB: one zip per run, the SDRFs, the `<PXD>.qpx.zip` files and `raw_files.raw.zip` with `raw_files.txt` and `raw_files.tsv`; `runs.tsv`, `extra.config` and Table S2 (`TableS6_run_parameters.tsv`) are not uploaded because they are in the code repository's Zenodo record, linked from the dataset description), then upload `<dir>` with `submission.px`.
5. Zenodo: `ZENODO_TOKEN=... python -m scripts.zenodo_upload` creates a draft (~20 GB, streamed
   from the FTP; everything but the 44 GB ProCan report), with metadata from
   `docs/zenodo_data_metadata.json`.
6. Code: `.zenodo.json` and `CITATION.cff` describe this repository for the GitHub-Zenodo release.
