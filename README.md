# quantmsdiann manuscript: analysis code and run provenance

Analysis scripts, run inputs and figure sources for the quantmsdiann paper.
quantmsdiann is the SDRF-driven Nextflow/DIA-NN workflow at
<https://github.com/bigbio/quantmsdiann> (v2.2.0 was used for this paper).

| Folder | Content |
|---|---|
| [`sdrf/`](sdrf/) | The SDRF used by every run, verbatim; per-run Nextflow parameters (`params.json`, `launch.txt`, `diann_config.cfg`); `runs.tsv` (index of all 38 runs); the FASTA databases (`fasta/`). Supplementary Table 1 with links: [`sdrf/README.md`](sdrf/README.md). |
| [`scripts/rebuild.py`](scripts/rebuild.py) | All main-text and Supplementary Note 4-6 figures and manuscript numbers, as named stages. |
| `scripts/review_*.py` | Supplementary Figures S6-S8 and Tables S3-S5 (added in revision). |
| `scripts/build_*.py`, `stage_pride_submission.py`, `zenodo_upload.py` | Table S6 and the PRIDE / Zenodo deposition of the results. |
| [`analysis/figures/`](analysis/figures/) | Figure sources (SVG) and the small derived tables they are drawn from. |
| [`methods.md`](methods.md) | The counting rule used for every identification number. |
| [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) | Figure/table -> public source -> command. |

## Data

Every number is computed from public data. Our reanalysis results (DIA-NN reports, matrices,
logs) are on the PRIDE FTP under
`https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/`
(and `quantms-collections/` for three cell-line cohorts), one folder per dataset and DIA-NN
version; they are also deposited as a PRIDE reanalysis dataset and on Zenodo (see the paper's
Data Availability). The PRIDE dataset is available at
<https://www.ebi.ac.uk/pride/archive/projects/PXD085732>. The original deposits we compare against are on PRIDE / MassIVE; the
accessions are listed in [`sdrf/README.md`](sdrf/README.md). Large inputs are never committed:
each stage downloads what it needs into `data/` (git-ignored) on first use.

### What is in PXD085732

The PRIDE dataset holds 88 files; raw files are not re-uploaded.

| PRIDE type | Files | Content |
|---|---|---|
| EXPERIMENTAL_DESIGN | 38 | The SDRF of every run, as used (`<dataset>_<version>_<name>.sdrf.tsv`). |
| SEARCH | 35 | One `<dataset>_<version>.zip` per run with a DIA-NN search: `quant_tables/` (report, matrices, stats), `pipeline_info/`, `sdrf/` configs, `pmultiqc/`, `run_metadata.json`, `nextflow.log`. |
| OTHER | 14 | 9 `<PXD>.qpx.zip` exports (PXD046357 and PXD049412 have one per DIA-NN version), and 5 `PXD071075_v2_5_0_sweep_*cores.zip` with only `launch.txt` and `nextflow_trace.txt`. |
| RAW | 1 | `raw_files.raw.zip`: `raw_files.txt` and `raw_files.tsv`, every raw file processed with its URL in the original deposit. |

There are more SDRFs (38) than search zips (35) because the two sets cover different runs:

- 33 runs have both an SDRF and a DIA-NN search zip.
- 5 SDRFs belong to the PXD071075 node-count sweep (10, 50, 100, 200 and 300 cores, Fig. 2a).
  These runs measured wall-clock time only, so their zips keep the Nextflow traces and no DIA-NN
  report (type OTHER). The five SDRFs are the same file, one per sweep run.
- 2 search zips are the standalone DIA-NN runs of the paired comparison
  (`Paired-standalone-PXD049412`, `Paired-standalone-PXD070049`, Supplementary Note 8). Plain
  DIA-NN takes no SDRF; their quantmsdiann counterparts carry the SDRF.

`runs.tsv`, `extra.config` and Table S6 are not in the PRIDE dataset; they are in this repository
and its Zenodo record (<https://doi.org/10.5281/zenodo.23322361>).

## Setup and use

```bash
conda env create -f environment.yml          # pins numpy<2, bundles rsvg-convert
conda activate quantmsdiann
python -m scripts.rebuild --list             # stages and what they produce
python -m scripts.rebuild --all              # data -> figures -> numbers
python -m scripts.rebuild --only fig2_validation
pytest tests/
```

(pip alternative: `pip install -r analysis/requirements.txt`.)

Identification counts follow one rule ([methods.md](methods.md) §1): global counts at
`Lib.Q.Value` / `Lib.PG.Q.Value` <= 0.01, per-run counts at `Q.Value` / `PG.Q.Value` <= 0.01,
decoys dropped, no other filter. plexDIA and phosphoproteomics runs show workflow support
(runtime) only; their identification counts are not benchmarked across DIA-NN versions.
