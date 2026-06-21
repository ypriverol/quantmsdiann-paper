"""Recompute `data/quantmsdiann_benchmarks/report_counts.tsv` from the public
PRIDE FTP reports, under the Vadim filter rule (see analysis.count_report_ids).

This replaces the old one-time cluster computation: every benchmark cohort's
full DIA-NN report is deposited on the FTP
(`quantmsdiann-benchmarks/proteobench/<dataset>/<version>/quant_tables/`), as
`diann_report.parquet` for DIA-NN >= 2.x and `diann_report.tsv` for 1.8.1. We
download whichever exists (cached under data/, gitignored), count it with the
rewritten `count_report`, and write the new-schema `report_counts.tsv`.

Run:  python -m analysis.recompute_report_counts
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from analysis.count_report_ids import (
    DATASET_MODULES, VERSIONS, PRECURSOR_Q, DEFAULT_PRECURSOR_Q,
    _load_report, count_report,
)
from analysis.figure_original_vs_quantmsdiann import download_if_missing

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "data" / "quantmsdiann_benchmarks"
FTP = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "quantmsdiann-benchmarks/proteobench"
)
OUT_COLS = [
    "dataset", "version", "prec_min1", "prec_min3", "prec_global",
    "prot_global", "prot_perrun_avg", "prot_complete", "peptides", "prot_2pep",
]


def fetch_report_dir(dataset: str, version: str) -> Path:
    """Return a local quant_tables dir holding the DIA-NN report for
    (dataset, version), downloading from the FTP if absent. 1.8.1 ships a
    `diann_report.tsv`; >= 2.x ships `diann_report.parquet`."""
    qt = ROOT / dataset / version / "quant_tables"
    parquet = qt / "diann_report.parquet"
    tsv = qt / "diann_report.tsv"
    if (parquet.exists() and parquet.stat().st_size > 0) or (
        tsv.exists() and tsv.stat().st_size > 0
    ):
        return qt
    fname = "diann_report.tsv" if version == "v1_8_1" else "diann_report.parquet"
    url = f"{FTP}/{dataset}/{version}/quant_tables/{fname}"
    print(f"  downloading {url}", flush=True)
    download_if_missing(url, qt / fname)
    return qt


def main() -> int:
    rows: list[dict] = []
    for dataset in DATASET_MODULES:
        for version in VERSIONS:
            print(f"{dataset} {version} ...", flush=True)
            qt = fetch_report_dir(dataset, version)
            df = _load_report(qt)
            c = count_report(df, precursor_q=PRECURSOR_Q.get(version, DEFAULT_PRECURSOR_Q))
            c.update(dataset=dataset, version=version)
            rows.append(c)
            print(f"  -> prec_global={c['prec_global']:,} prot_global={c['prot_global']:,} "
                  f"prot_perrun_avg={c['prot_perrun_avg']:,} prec_min1={c['prec_min1']:,}",
                  flush=True)
    out = ROOT / "report_counts.tsv"
    pd.DataFrame(rows)[OUT_COLS].to_csv(out, sep="\t", index=False)
    print(f"wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
