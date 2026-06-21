#!/usr/bin/env python
"""Generate the Fig. 3 single-cell input tables from the PUBLIC DIA-NN reports.

This is the reproducibility generator for ``figure_single_cell_combined.py``:
every number in Fig. 3 is derived here from the deposited DIA-NN reports, so
nothing in the figure is hand-entered.

PROVENANCE
==========
Our reanalyses (this pipeline) are published on the PRIDE FTP:

  https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/single-cell/<PXD>/<version>/quant_tables/diann_report.{parquet,tsv}

    PXD046357  HeLa Astral single-cell (Orbitrap Astral) -> "HeLa Astral SC"

  DIA-NN versions: ``v1_8_1`` (diann_report.tsv) and ``v2_5_1_enterprise``
  (diann_report.parquet).

Pipeline:        https://github.com/bigbio/quantmsdiann
SDRF tooling:    https://github.com/bigbio/sdrf-pipelines (convert-diann)
Counting logic:  analysis/count_report_ids.py (canonical, reused here).

Counting convention (Vadim filter rule; see methods.md §1)
----------------------------------------------------------
* totals are GLOBAL numbers: protein groups = distinct ``Protein.Group`` at
  ``Lib.PG.Q.Value <= 0.01`` (``prot_global``); precursors = distinct
  ``Precursor.Id`` at ``Lib.Q.Value <= 0.01`` (``prec_global``). No
  contaminant/target filter; zeros counted.
* the per-cell and completeness panels are PER-RUN numbers: protein groups per
  run at ``PG.Q.Value <= 0.01`` only (no target/global filter).

Dynamic range / CV use ``PG.MaxLFQ`` of the per-run protein groups (a
quantitative metric, not an identification count).

NOTE: this counts from the *report* (not the ``*_matrix.tsv`` files); the
matrices bake in ``--matrix-spec-q`` at a version-dependent run q-value and are
not comparable across versions (see count_report_ids.py docstring).

Outputs (data/single_cell/): mv_per_cell.tsv, mv_completeness.tsv,
mv_rank_abundance.tsv, mv_cv.tsv, sc_totals.tsv.

Run:  python -m analysis.make_single_cell_tables
      (downloads ~0.5 GB of reports once; cached under data/single_cell/cache/)
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.count_report_ids import (
    PRECURSOR_Q, DEFAULT_PRECURSOR_Q, Q_THRESHOLD, count_report,
)

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "data" / "single_cell"
CACHE_DIR = OUT_DIR / "cache"

FTP_BASE = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/single-cell"
)
# dataset display name -> PRIDE accession (for labels) and FTP sub-directory.
ACC = {"HeLa Astral SC": "PXD046357", "A549/H460 SC": "PXD049412"}
FTP_DIR = {"HeLa Astral SC": "PXD046357", "A549/H460 SC": "PXD049412"}
VERSIONS = ["1_8_1", "2_5_1_enterprise"]
FLAG = "HeLa Astral SC"  # the dataset carrying the depth/completeness panels

_COLS = [
    "Run", "Precursor.Id", "Protein.Group", "Q.Value",
    "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.MaxLFQ", "Decoy",
]


def _report_url(ftp_dir: str, version: str) -> str:
    ext = "tsv" if version == "1_8_1" else "parquet"
    return f"{FTP_BASE}/{ftp_dir}/v{version}/quant_tables/diann_report.{ext}"


def _cached_report(ftp_dir: str, version: str) -> Path:
    """Download the deposited DIA-NN report once and cache it on disk."""
    url = _report_url(ftp_dir, version)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / f"{ftp_dir}_v{version}_{Path(url).name}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"Downloading {url} (cached) ...", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=900) as resp, open(tmp, "wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)
    return dest


def _load(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(path).schema_arrow.names)
        return pq.read_table(path, columns=[c for c in _COLS if c in have]).to_pandas()
    return pd.read_csv(path, sep="\t", usecols=lambda c: c in _COLS, low_memory=False)


def _perrun_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Per-run protein-group rows under the Vadim rule: PG.Q.Value <= 1% only
    (no contaminant/target filter, no global filter). Decoys dropped."""
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    return df[df["PG.Q.Value"] <= Q_THRESHOLD]


def build() -> dict[str, pd.DataFrame]:
    per_cell, completeness, rank, cv, totals = [], [], [], [], []
    for ds, acc in ACC.items():
        for version in VERSIONS:
            df = _load(_cached_report(FTP_DIR[ds], version))
            # PXD049412 single cells only: drop the 20x/40x A549 carrier (library)
            # runs from per-cell and total counts (no-op for cohorts without them).
            df = df[~df["Run"].astype(str).str.contains("20xSC|40xSC", case=False, regex=True)]
            prot = _perrun_proteins(df)
            pgrun = prot.drop_duplicates(["Run", "Protein.Group"])
            # totals are GLOBAL numbers (Vadim rule): precursors -> prec_global
            # (Lib.Q.Value), protein groups -> prot_global (Lib.PG.Q.Value).
            # The per-cell distribution (panel A) and the completeness curve
            # below are PER-RUN numbers (PG.Q.Value only), from `pgrun`.
            c = count_report(df, precursor_q=PRECURSOR_Q.get(f"v{version}", DEFAULT_PRECURSOR_Q))
            totals.append((ds, version, c["prec_global"], c["prot_global"]))
            for run, g in pgrun.groupby("Run"):
                per_cell.append((ds, version, int(g["Protein.Group"].nunique())))
            # completeness + CV for BOTH datasets (panels B, C now show both).
            n_runs = pgrun["Run"].nunique()
            seen = pgrun.groupby("Protein.Group")["Run"].nunique()
            for mc in range(1, n_runs + 1):
                completeness.append((ds, version, mc, int((seen >= mc).sum())))
            q = pgrun.dropna(subset=["PG.MaxLFQ"]).copy()
            q["PG.MaxLFQ"] = pd.to_numeric(q["PG.MaxLFQ"], errors="coerce")
            q = q[q["PG.MaxLFQ"] > 0]
            agg = q.groupby("Protein.Group")["PG.MaxLFQ"].agg(["mean", "std", "count"])
            agg = agg[agg["count"] >= 3]
            for _, r in agg.iterrows():
                if r["mean"] > 0 and not np.isnan(r["std"]):
                    cv.append((ds, version, float(r["std"] / r["mean"])))
            # rank/dynamic-range kept for the flagship only (panel retired).
            if ds == FLAG:
                mean_int = q.groupby("Protein.Group")["PG.MaxLFQ"].mean().sort_values(ascending=False)
                for i, val in enumerate(mean_int.values, start=1):
                    if i == 1 or i % 10 == 0:
                        rank.append((version, i, float(np.log10(val))))
    return {
        "mv_per_cell.tsv": pd.DataFrame(per_cell, columns=["dataset", "version", "pg_count"]),
        "mv_completeness.tsv": pd.DataFrame(completeness, columns=["dataset", "version", "min_cells", "n_proteins"]),
        "mv_rank_abundance.tsv": pd.DataFrame(rank, columns=["version", "rank", "log10_intensity"]),
        "mv_cv.tsv": pd.DataFrame(cv, columns=["dataset", "version", "cv"]),
        "sc_totals.tsv": pd.DataFrame(totals, columns=["dataset", "version", "precursors", "proteins"]),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, frame in build().items():
        frame.to_csv(OUT_DIR / name, sep="\t", index=False)
        print(f"wrote {OUT_DIR / name} ({len(frame)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
