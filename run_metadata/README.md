# Final pipeline parameters of every run

One folder per run used in the manuscript (`<dataset>/<run>/`):

| File | Content |
|---|---|
| `params.json` | Final Nextflow/nf-core parameters as resolved by the pipeline (`pipeline_info/params_*.json`) |
| `launch.txt` | The `nextflow run` command (from `pipeline_report.txt`), or the SLURM run header when that report is not available |
| `diann_config.cfg` | DIA-NN options the pipeline derived from the SDRF |

`runs.tsv` indexes all 38 runs: SDRF used (in `../sdrf/`), where the
parameters came from, and the key search parameters. Supplementary Table S6
(`analysis/figures/review/supplementary_revision/TableS6_run_parameters.tsv`)
is the readable version. `common/extra.config` is the resource-override config
passed with `-c` to every cluster run.

## Runs without a `params.json`

The phosphoproteomics runs (PXD034128, PXD034623, PXD049692 v2.5.1) and the
PXD071075 v2.5.0 node-count sweep published no `pipeline_info/` folder. For
these, `launch.txt` is the SLURM run header, which records the SDRF, FASTA,
DIA-NN version, Nextflow profiles and the extra arguments given on the command
line; every other parameter is the pipeline default.

## Rebuilding

    python -m scripts.build_run_metadata

It reads the run folders on the PRIDE FTP, plus `run_metadata/_cluster/` for
the runs above. `_cluster/` is not committed: it holds the SLURM logs and input
SDRFs copied read-only from `/hps/nobackup/juan/pride/reanalysis` on the
EMBL-EBI cluster (paths in `scripts/build_run_metadata.py`, `CLUSTER_RUNS`).
