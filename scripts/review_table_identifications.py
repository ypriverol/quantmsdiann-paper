"""Supplementary Table 8: precursor and peptide identifications per dataset and DIA-NN version.

One rule for every row (methods.md section 1, global counts): distinct Precursor.Id and
distinct Stripped.Sequence at Lib.Q.Value <= 0.01, decoys dropped, nothing else.

Sources:
  - ProteoBench modules: data/quantmsdiann_benchmarks/report_counts.tsv (stage report_counts,
    same rule via scripts.rebuild.count_report);
  - every other row: the DIA-NN report on the PRIDE FTP (cell-line cohorts through
    scripts.cellline_reports, so a local copy is used when present). For the single-cell
    cohorts the 20xSC / 40xSC bulk-standard runs are excluded, as in Fig. 3.

Each report is streamed once (only four columns) and its counts cached as JSON under
data/identifications_cache/.

Usage: python -m scripts.review_table_identifications
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import fsspec
import pandas as pd
import pyarrow.parquet as pq

from scripts.cellline_reports import open_report

ROOT = Path(__file__).resolve().parents[1]
FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks"
CACHE = ROOT / "data" / "identifications_cache"
OUT = ROOT / "analysis/figures/review/supplementary_revision/TableS8_identifications.tsv"
Q = 0.01
COLS = ["Run", "Precursor.Id", "Stripped.Sequence", "Lib.Q.Value", "Decoy"]
SC_EXCLUDE = re.compile(r"20xSC|40xSC", re.I)
VERSION = {"v1_8_1": "1.8.1", "v2_5_1": "2.5.1", "v2_5_1_enterprise": "2.5.1 (enterprise)"}

PROTEOBENCH = [("PXD049412", "PXD049412 (Module 9)"), ("PXD062685", "PXD062685 (Module 5)"),
               ("PXD070049", "PXD070049 (Module 10)"), ("ProteoBench_Module_7", "Module 7")]
CELL_LINES = ["PXD003539", "PXD030304", "PXD004701", "PXD041421", "PXD017199"]
# (section, label, FTP folder, version, exclude bulk-standard runs)
FTP_ROWS = [
    ("Single-cell", "PXD046357", "single-cell/PXD046357", "v1_8_1", True),
    ("Single-cell", "PXD046357", "single-cell/PXD046357", "v2_5_1_enterprise", True),
    ("Single-cell", "PXD049412", "single-cell/PXD049412", "v1_8_1", True),
    ("Single-cell", "PXD049412", "single-cell/PXD049412", "v2_5_1_enterprise", True),
    ("Single-cell", "MSV000093870", "single-cell/MSV000093870", "v2_5_1", False),
    ("Phosphoproteomics", "PXD034128 (biological study)", "phospho/PXD034128-biological-study", "v2_5_1", False),
    ("Phosphoproteomics", "PXD034128 (biological study)", "phospho/PXD034128-biological-study", "v2_5_1_enterprise", False),
    ("Phosphoproteomics", "PXD034128 (high-speed DIA)", "phospho/PXD034128-highspeed-DIA", "v2_5_1", False),
    ("Phosphoproteomics", "PXD034128 (high-speed DIA)", "phospho/PXD034128-highspeed-DIA", "v2_5_1_enterprise", False),
    ("Phosphoproteomics", "PXD034623", "phospho/PXD034623", "v2_5_1", False),
    ("Phosphoproteomics", "PXD034623", "phospho/PXD034623", "v2_5_1_enterprise", False),
    ("Phosphoproteomics", "PXD049692", "phospho/PXD049692", "v2_5_1", False),
    ("Phosphoproteomics", "PXD049692", "phospho/PXD049692", "v2_5_1_enterprise", False),
]


def chunks(fh, tsv: bool):
    if tsv:  # DIA-NN 1.8.1 deposits a TSV report
        yield from pd.read_csv(fh, sep="\t", usecols=lambda c: c in COLS, chunksize=2_000_000)
        return
    pf = pq.ParquetFile(fh)
    cols = [c for c in COLS if c in pf.schema_arrow.names]
    for i in range(pf.num_row_groups):
        yield pf.read_row_group(i, columns=cols).to_pandas()


def count(fh, exclude: bool, tsv: bool = False) -> dict:
    """Distinct precursors and peptides at Lib.Q.Value <= Q, decoys dropped."""
    prec, pep = set(), set()
    for d in chunks(fh, tsv):
        if "Decoy" in d:
            d = d[d["Decoy"] == 0]
        if exclude:
            d = d[~d["Run"].astype(str).str.contains(SC_EXCLUDE)]
        d = d[d["Lib.Q.Value"] <= Q]
        prec.update(d["Precursor.Id"])
        pep.update(d["Stripped.Sequence"])
    return {"precursors": len(prec), "peptides": len(pep)}


def cached(key: str, opener, exclude: bool, tsv: bool = False) -> dict:
    path = CACHE / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    with opener() as fh:
        c = count(fh, exclude, tsv)
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(c))
    print(key, c, flush=True)
    return c


def main() -> int:
    rows = []
    rc = pd.read_csv(ROOT / "data/quantmsdiann_benchmarks/report_counts.tsv", sep="\t")
    for ds, label in PROTEOBENCH:
        for v, name in VERSION.items():
            r = rc[(rc["dataset"] == ds) & (rc["version"] == v)].iloc[0]
            rows.append(("ProteoBench", label, name, int(r["prec_global"]), int(r["peptides"])))
    for ds in CELL_LINES:
        c = cached(f"{ds}_v2_5_1", lambda ds=ds: open_report(ds), False)
        rows.append(("Bulk cell-line", ds, "2.5.1", c["precursors"], c["peptides"]))
    for section, label, folder, v, exclude in FTP_ROWS:
        tsv = v == "v1_8_1"
        url = f"{FTP}/{folder}/{v}/quant_tables/diann_report.{'tsv' if tsv else 'parquet'}"
        key = f"{folder.replace('/', '_')}_{v}"
        c = cached(key, lambda url=url: fsspec.open(url, block_size=2**22).open(), exclude, tsv)
        rows.append((section, label, VERSION[v], c["precursors"], c["peptides"]))
    out = pd.DataFrame(rows, columns=["section", "accession", "version", "precursors", "peptides"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
