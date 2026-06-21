"""Count precursors and protein groups from the DIA-NN *report*, not the matrices.

The `*_matrix.tsv` outputs bake in `--matrix-spec-q` (0.05 run-specific) and,
because the quantmsdiann pipeline sets `--qvalue` to 0.01 for DIA-NN 1.8.1 but
0.05 for 2.5.1/enterprise, matrix row counts are filtered at *different*
run-specific q-values per version and are NOT comparable across versions.

This module counts identifications directly from the per-precursor report
(`diann_report.parquet` for DIA-NN >= 2.x, `diann_report.tsv` for 1.8.1).

Filter rule (Vadim review, 2026-06-21 -- scientific-correctness requirement for
submission). Reported quantities fall into two classes, each with EXACTLY one
admissible filter and nothing else (no contaminant/target filter, no
positive-quantity filter, zeros counted):

  * Per-run numbers (within a single run / cell):
      - protein groups: `PG.Q.Value <= 0.01` only.
      - precursors:     `Q.Value <= precursor_q` only (run-specific `--qvalue`:
        0.01 for 1.8.1, 0.05 for >= 2.5.0, DIA-NN's per-version operating point).
  * Global numbers (dataset union / totals):
      - protein groups: `Lib.PG.Q.Value <= 0.01` only.
      - precursors:     `Lib.Q.Value <= 0.01` only.

Emitted metrics:
  * `prec_min1` / `prec_min3` -- precursors identified in >= 1 / >= 3 runs at the
    per-run `Q.Value` cut-off (replicate-reproducibility view).
  * `prec_global` -- distinct precursors with `Lib.Q.Value <= 0.01` (global total).
  * `prot_global` -- distinct protein groups with `Lib.PG.Q.Value <= 0.01`.
  * `prot_perrun_avg` / `prot_complete` -- average protein groups per run, and
    protein groups quantified in *every* run, at run-specific `PG.Q.Value <= 0.01`.
  * `peptides` -- distinct stripped peptides among `Lib.Q.Value`-passing precursors.
  * `prot_2pep` -- global protein groups (`Lib.PG.Q.Value`) with >= 2 distinct
    stripped peptides among `Lib.Q.Value`-passing precursors.

Decoys (`Decoy == 1`) are dropped everywhere: this removes the FDR null model and
is not one of the forbidden "filters".

The reports are multi-GB; this is run on the downloaded reports (from the public
PRIDE FTP) and its small output, `report_counts.tsv`, is staged into
`data/quantmsdiann_benchmarks/` and consumed by
`figure_quantmsdiann_benchmarks_vs_proteobench.py`.

Usage:
    python -m analysis.count_report_ids \
        --results-root /path/to/quantmsdiann_results \
        --out report_counts.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

Q_THRESHOLD = 0.01           # per-run global-q cut-offs and the Lib.* cut-offs
# DIA-NN's recommended run-specific precursor q-value (`--qvalue`) per release.
PRECURSOR_Q = {
    "v1_8_1": 0.01,
    "v2_5_1": 0.05,
    "v2_5_1_enterprise": 0.05,
}
DEFAULT_PRECURSOR_Q = 0.01
_NEEDED_COLS = [
    "Run", "Precursor.Id", "Protein.Group", "Stripped.Sequence",
    "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Decoy",
]

DATASET_MODULES = (
    "ProteoBench_Module_7", "PXD049412", "PXD062685", "PXD070049",
)
VERSIONS = ("v1_8_1", "v2_5_1", "v2_5_1_enterprise")


def _load_report(report_dir: Path) -> pd.DataFrame:
    """Read the DIA-NN report from `report_dir`, preferring the parquet
    (DIA-NN >= 2.x) and falling back to the classic `diann_report.tsv`
    (1.8.1). Only the columns we need are read."""
    parquet = report_dir / "diann_report.parquet"
    tsv = report_dir / "diann_report.tsv"
    if parquet.exists():
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(parquet).schema_arrow.names)
        cols = [c for c in _NEEDED_COLS if c in have]
        return pq.read_table(parquet, columns=cols).to_pandas()
    if tsv.exists():
        return pd.read_csv(
            tsv, sep="\t", usecols=lambda c: c in _NEEDED_COLS,
        )
    raise FileNotFoundError(
        f"no diann_report.parquet or .tsv in {report_dir}"
    )


def count_report(df: pd.DataFrame, precursor_q: float = DEFAULT_PRECURSOR_Q) -> dict[str, int]:
    """Compute precursor and protein-group counts from a DIA-NN report frame
    under the Vadim filter rule (see module docstring): per-run protein groups
    on `PG.Q.Value` only; per-run precursors on `Q.Value` only; global protein
    groups on `Lib.PG.Q.Value` only; global precursors on `Lib.Q.Value` only.
    No contaminant/target filter, no positive-quantity filter (zeros counted)."""
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]

    # --- per-run / replicate precursors: Q.Value only --------------------
    passing = df[df["Q.Value"] <= precursor_q]
    n_runs_per_prec = passing.groupby("Precursor.Id")["Run"].nunique()
    prec_min1 = int((n_runs_per_prec >= 1).sum())
    prec_min3 = int((n_runs_per_prec >= 3).sum())

    # --- global precursors: Lib.Q.Value only -----------------------------
    prec_global = 0
    if "Lib.Q.Value" in df.columns:
        prec_global = int(
            df.loc[df["Lib.Q.Value"] <= Q_THRESHOLD, "Precursor.Id"].nunique()
        )

    # --- global protein groups: Lib.PG.Q.Value only ----------------------
    prot_global = 0
    if "Lib.PG.Q.Value" in df.columns:
        prot_global = int(
            df.loc[df["Lib.PG.Q.Value"] <= Q_THRESHOLD, "Protein.Group"].nunique()
        )

    # --- per-run protein groups: PG.Q.Value only (zeros counted) ---------
    prot_perrun_avg = prot_complete = 0
    if "PG.Q.Value" in df.columns:
        n_runs = df["Run"].nunique()
        pg = df[df["PG.Q.Value"] <= Q_THRESHOLD][["Run", "Protein.Group"]].drop_duplicates()
        per_run = pg.groupby("Run")["Protein.Group"].nunique()
        prot_perrun_avg = int(round(per_run.mean())) if len(per_run) else 0
        in_n_runs = pg.groupby("Protein.Group")["Run"].nunique()
        prot_complete = int((in_n_runs == n_runs).sum()) if n_runs else 0

    # --- peptide-level metrics (global rule) -----------------------------
    peptides = prot_2pep = 0
    if "Stripped.Sequence" in df.columns and "Lib.Q.Value" in df.columns:
        dq = df[df["Lib.Q.Value"] <= Q_THRESHOLD]
        peptides = int(dq["Stripped.Sequence"].nunique())
        if "Lib.PG.Q.Value" in df.columns:
            dp = dq[dq["Lib.PG.Q.Value"] <= Q_THRESHOLD]
            per_prot = dp.groupby("Protein.Group")["Stripped.Sequence"].nunique()
            prot_2pep = int((per_prot >= 2).sum())

    return {
        "prec_min1": prec_min1, "prec_min3": prec_min3,
        "prec_global": prec_global, "prot_global": prot_global,
        "prot_perrun_avg": prot_perrun_avg, "prot_complete": prot_complete,
        "peptides": peptides, "prot_2pep": prot_2pep,
    }


def main(argv: list[str] | None = None) -> int:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-root", required=True, type=Path,
                    help="quantmsdiann_results dir holding <module>/<version>/quant_tables")
    ap.add_argument("--out", required=True, type=Path,
                    help="output report_counts.tsv path")
    ap.add_argument("--results-suffix", default="",
                    help="suffix for the per-version dir, e.g. '_relaxed' reads "
                         "<module>/<version-without-v><suffix>/quant_tables (the "
                         "--relaxed-prot-inf re-run dirs '1_8_1_relaxed' etc.)")
    args = ap.parse_args(argv)

    rows: list[dict] = []
    for module in DATASET_MODULES:
        for version in VERSIONS:
            if args.results_suffix:
                base = version[1:] if version.startswith("v") else version
                ver_dir = f"{base}{args.results_suffix}"
            else:
                ver_dir = version
            rdir = args.results_root / module / ver_dir / "quant_tables"
            try:
                df = _load_report(rdir)
            except FileNotFoundError as exc:
                print(f"WARN: {exc}", file=sys.stderr)
                continue
            c = count_report(df, precursor_q=PRECURSOR_Q.get(version, DEFAULT_PRECURSOR_Q))
            c.update(dataset=module, version=version)
            rows.append(c)
            print(f"{module} {version}: {c['prec_global']:,} prec (global) / "
                  f"{c['prot_global']:,} proteins (global, Lib.q<=0.01)")
    cols = ["dataset", "version", "prec_min1", "prec_min3", "prec_global",
            "prot_global", "prot_perrun_avg", "prot_complete",
            "peptides", "prot_2pep"]
    pd.DataFrame(rows)[cols].to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
