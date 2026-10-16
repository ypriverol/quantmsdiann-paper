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
