"""DIA-NN reports of the five cell-line cohorts, all processed the same way: quantmsdiann
de98e98 (v2.2.0dev), DIA-NN 2.5.1, Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta
(no entrapment), the same search parameters (sdrf/runs.tsv).

Each report is read from data/cell-lines-v2_5_1/<PXD>/ when present (copied from the Codon
run folders; the rebuild never purges this folder), otherwise from the PRIDE FTP path the
run is published under. Derived tables are cached under data/cell-lines-v2_5_1/integration/.
"""
from __future__ import annotations

from pathlib import Path

import fsspec

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
LOCAL = DATA / "cell-lines-v2_5_1"  # harmonised reports (<PXD>/) and derived-table caches
FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks/cell-lines"
COHORTS = ("PXD003539", "PXD004701", "PXD030304", "PXD017199", "PXD041421")


def report_url(ds: str, name: str = "diann_report.parquet") -> str:
    return f"{FTP}/{ds}/v2_5_1/quant_tables/{name}"


def open_report(ds: str, name: str = "diann_report.parquet"):
    """Binary file object for a cohort's report file (local copy first, then FTP)."""
    local = LOCAL / ds / name
    if local.exists():
        return local.open("rb")
    return fsspec.open(report_url(ds, name), block_size=2**22).open()


def source(ds: str, name: str = "diann_report.parquet") -> str:
    local = LOCAL / ds / name
    return str(local.relative_to(ROOT)) if local.exists() else report_url(ds, name)
