# Reproducibility & data provenance

Every number and figure in the manuscript is derived from **public** data by a
script in this repository. The large DIA-NN reports themselves are *not* checked
in (they are 10s–100s of MB each); instead each figure script either downloads
them on first run and caches them, or reads a small derived table produced by a
generator that does. This file maps every figure to its source data and the
command that regenerates it.

## Data sources

| Source | Base URL / location | What |
|---|---|---|
| **Our reanalyses** (this pipeline) | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/` | DIA-NN reports + matrices for every dataset we reprocessed, one tree per `DIA-NN` version (`v1_8_1`, `v2_5_1`, `v2_5_1_enterprise`), under `quant_tables/`. Subtrees: `single-cell/`, `proteobench/`, `PXD071075_cluster_sizes/`, plus per-PXD dirs (spatial, phospho). |
| **Original deposits** | `https://ftp.pride.ebi.ac.uk/pride/data/archive/<YYYY>/<MM>/<PXD>/` | The authors' originally deposited results we compare against (e.g. PXD003539 OpenSWATH matrix, 2020/06; PXD064049 `DIANN_results.zip`, 2025/07). |
| **quantms collections** | `https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/` | Earlier curated reanalysis collections reused for some cohort figures (e.g. ProCan `proteo-search-data/ProCan/`). |
| **Pipeline** | `https://github.com/bigbio/quantmsdiann` | The Nextflow/DIA-NN workflow. |
| **SDRF tooling** | `https://github.com/bigbio/sdrf-pipelines` | `convert-diann` (SDRF → per-run DIA-NN config). |
| **Containers** | `https://github.com/bigbio/quantms-containers` | Per-`DIA-NN`-version Singularity/Docker recipes. |

## Counting conventions (the one place the rules live)

The canonical rule is [methods.md §1](methods.md); the primitive that implements it
is `count_report` in [`scripts/rebuild.py`](scripts/rebuild.py). Exactly one
admissible q-value filter per quantity, and nothing else:

- **protein groups, per run**: `PG.Q.Value ≤ 0.01`; **precursors, per run**: `Q.Value ≤ 0.01`.
- **protein groups, global** (unique per dataset): `Lib.PG.Q.Value ≤ 0.01`;
  **precursors, global**: `Lib.Q.Value ≤ 0.01`.
- **No** contaminant/target filter, **no** positive-quantity filter (zeros counted);
  decoys dropped. A §4 guard flags datasets where >1.2% of protein groups are
  multi-accession (possible inference inflation).
- counts come from the **report** (`diann_report.parquet`/`.tsv`), never the `*_matrix.tsv`
  files (those bake in a version-dependent `--matrix-spec-q` and are not comparable across versions).
- plexDIA and phosphoproteomics are processed for workflow support only; their
  identification counts are not benchmarked across DIA-NN versions (runtime figures only).

## Public availability (verified 2026-06-17)

Reproducibility from **public** data depends on the reanalysis reports being
deposited under `quantmsdiann-benchmarks/`. Verified by HTTP HEAD:

| Dataset / figure | Public on FTP? |
|---|---|
| single-cell PXD046357 (HeLa Astral) / PXD049412 (A549/H460) (Fig 3) | ✅ live |
| ProteoBench modules — Module_7, PXD062685, PXD049412, PXD070049 (Fig 2) | ✅ live |
| PXD071075 cluster-size sweep (Fig 1) | ✅ live |
| NCI-60 PXD003539 (Fig 4, via `quantms-collections`) | ✅ live |
| ProCan PXD030304 (Fig 4) — `quantms-collections/absolute-expression-2.0/cell-lines/PXD030304/` (figure source; also mirrored at `quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1/`) | ✅ live |
| plexDIA MSV000093870 (support only, not benchmarked) — `single-cell/MSV000093870/v2_5_1/` | ✅ live |
| spatial PXD064049 (Fig 4) — `spatial/PXD064049/v2_5_1_enterprise/` | ✅ live |
| phospho PXD049692 / PXD034128 / PXD034623 (Supp) — `phospho/<ds>/v{2_5_1,2_5_1_enterprise}/` | ✅ live |

All reanalysis reports are now deposited under `quantmsdiann-benchmarks/`
(verified by HTTP HEAD). Each figure's generator pulls its report from the FTP:
single-cell → the `single_cell_tables` stage; spatial → the figure stage
downloads directly. ProCan's per-tissue
generator is the one remaining item.

## Run inputs and parameters (committed)

- `sdrf/`: the exact SDRF of every run, copied verbatim (the curated copies in
  sdrf-annotated-datasets keep changing). Known instrument-accession errors are
  listed in `sdrf/README.md`.
- `sdrf/<folder>/<run>/`: final Nextflow parameters (`params.json`),
  launch command (`launch.txt`) and DIA-NN options (`diann_config.cfg`) per run,
  next to the run's SDRF; `sdrf/runs.tsv` indexes all 38 runs. Rebuild with
  `python -m scripts.build_run_metadata`.

## Reviewer-response analyses

| Script | Output (`analysis/figures/review/`) |
|---|---|
| `review_desktop_vs_quantmsdiann` | paired standalone vs quantmsdiann statistics |
| `review_desktop_overview` | `supplementary_revision/FigS6_*` |
| `review_runtime_vs_scale` | wall-clock vs files / GB |
| `review_threshold_regimes`, `review_peptide_overlap` | peptide-evidence counts and peptide overlap (streams the Fig. 4 reports) |
| `review_supplementary` | `supplementary_revision/FigS7-S8`, `TableS3-S5` |
| `build_run_metadata` | `supplementary_revision/TableS6_run_parameters.tsv` |

The paired runs read `data/proteobench-review-experiments/x/` (extracted from
the review-experiment archives; git-ignored).

## Deposition (PRIDE reanalysis dataset + Zenodo)

1. `python -m scripts.build_pride_raw_list` writes `pride_submission/raw_files.{tsv,txt}`.
2. `python -m scripts.build_submission_manifest` crawls the run folders on the FTP.
3. `python -m scripts.build_submission_px` writes the corrected SDRF copies,
   `submission_files_final.tsv`, `zenodo_files.tsv` and `submission.px`.
4. PRIDE: on the EBI cluster, `python -m scripts.stage_pride_submission --dest <dir>`
   (~64 GB), then upload `<dir>` with `submission.px`.
5. Zenodo: `ZENODO_TOKEN=... python -m scripts.zenodo_upload` creates a draft
   (~20 GB, streamed from the FTP; everything but the 44 GB ProCan report),
   with metadata from `docs/zenodo_data_metadata.json`. Publish it on zenodo.org.
6. Code: `.zenodo.json` describes this repository for the GitHub-Zenodo release.

## Per-figure provenance

| Figure | Numbers | Source (public) | Regenerate | Status |
|---|---|---|---|---|
| **Fig 1** pipeline/scaling | wall-clock per node-count & per dataset | `quantmsdiann-benchmarks/PXD071075_cluster_sizes/` (queue sweep) + per-dataset `nextflow_trace.txt` / `run_metadata.json` | `python -m scripts.rebuild --only performance_trace queue_sweep mdc_cluster_runtime` | trace/metadata staged from FTP |
| **Fig 2** ProteoBench | depth vs accuracy, per-version IDs | ProteoBench community submissions + our `quantmsdiann-benchmarks/proteobench/` reports | `python -m scripts.rebuild --only benchmarks proteobench_accuracy fig2_validation` | community JSONs staged |
| **Fig 3** single-cell | per-cell / completeness (both cohorts) / dynamic range / CV / totals | `quantmsdiann-benchmarks/single-cell/{PXD046357,PXD049412}/v{1_8_1,2_5_1_enterprise}/quant_tables/diann_report.{tsv,parquet}` | `python -m scripts.rebuild --only single_cell_tables single_cell_combined` | **fully reproducible** ✅ |
| **Fig 4** reanalysis / cell-line atlas | per-cohort protein groups, overlaps | `data/PXD00{3539,4701},PXD030304,...` reports (from `quantms-collections`/`quantmsdiann-benchmarks`) + per-cohort report JSONs | `python -m scripts.rebuild --only reanalysis_improvement_data reanalysis_improvement atlas` | **fully reproducible** ✅ |
| **Supp Fig S3** runtime by dataset type | per-run wall-clock by dataset type (incl. plexDIA + phospho, processed for support but not benchmarked) | `parallelism_data.tsv` (nextflow traces) | `python -m scripts.rebuild --only runtime_by_dataset` | **fully reproducible** ✅ |
| **Fig 4** spatial DVP (PXD064049) | protein groups + precursors, deposited 1.8.1 vs reanalysis, §1 report rule | `quantmsdiann-benchmarks/spatial/PXD064049/v2_5_1_enterprise/quant_tables/diann_report.parquet` + original `2025/07/PXD064049/DIANN_results.zip` (1.8.1 `MYCN_High_Low.tsv`) | `python -m scripts.rebuild --only reanalysis_improvement_data reanalysis_improvement` | **fully reproducible** ✅ |

## Original / baseline counts

The reanalysis-recovery "original" counts are now **recomputed from public
deposits** (no transcribed paper headlines), each via `scripts/rebuild.py`:

- **Guo 2019** (NCI-60, PXD003539) — recomputed from the deposited PRIDE OpenSWATH
  `feature_alignment_requant_matrix.tsv` (≥2-peptide = 4,284; all-protein = 6,556).
- **Gonçalves 2022 / ProCan** (PXD030304) — recomputed from the figshare
  peptide-count matrix (doi:10.6084/m9.figshare.19345397) (≥2-peptide = 6,698;
  all-protein = 8,498).
- **Sun 2023** (PXD004701) — recomputed from the deposited PRIDE OpenSWATH
  `1R_all_id.txt` + library (proteotypic = 6,183; all-protein = 6,696).
- **Galatidou 2024** (MSV000093870 plexDIA) — processed for workflow support only;
  not benchmarked against the deposit (DIA-NN identification counts are not comparable
  across versions), so it carries no baseline count.

Still cited as constants (no downloadable processed source): **Walzer 2022**
(PXD003539, Supplementary Table S2) and the library-size / cell-line-split
metadata. Each appears as a clearly-labelled constant with a source comment.

---
*Generators download from the FTP cache under `data/**/cache/` (git-ignored).
Run any stage from the repo root via `python -m scripts.rebuild --only <stage>`
(`--list` shows them all); `python -m scripts.rebuild --all` runs the full chain.*
