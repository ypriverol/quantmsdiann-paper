"""Build the file manifest for the PRIDE reanalysis submission and the Zenodo deposit.

Crawls the published run directories on the PRIDE FTP (only the runs used in
the manuscript), classifies every file into a PRIDE file type, gets exact sizes
(HTTP HEAD), and writes:

  pride_submission/submission_files.tsv   one row per file to upload
  pride_submission/submission.px          PRIDE submission summary (metadata + file list)

Each file is renamed <dataset>_<diann-version>_<file> so files from different
runs do not collide in PRIDE's flat file list. MSstats input CSVs, QPX
exports (added as zips by scripts.build_submission_px), spectral libraries and
per-run .quant files are left out here: they are regenerable and dominate the size.

Usage:
    python -m scripts.build_submission_manifest
"""

from __future__ import annotations

import csv
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes"
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "pride_submission"

# (dataset label, version label, FTP run directory). Only runs used in the paper.
RUNS = [
    *[(f"ProteoBench-{d}", v, f"quantmsdiann-benchmarks/proteobench/{d}/{v}")
      for d in ("PXD049412", "PXD062685", "PXD070049", "ProteoBench_Module_7")
      for v in ("v1_8_1", "v2_5_1", "v2_5_1_enterprise")],
    *[(d, v, f"quantmsdiann-benchmarks/single-cell/{d}/{v}")
      for d in ("PXD046357", "PXD049412") for v in ("v1_8_1", "v2_5_1_enterprise")],
    ("MSV000093870", "v2_5_0", "quantmsdiann-benchmarks/single-cell/MSV000093870_plexDIA/v2_5_0"),
    ("PXD003539", "v2_5_1_enterprise", "quantmsdiann-benchmarks/cell-lines/PXD003539/v2_5_1_enterprise"),
    ("PXD004701", "v2_5_1", "quantmsdiann-benchmarks/cell-lines/PXD004701/v2_5_1"),
    ("PXD030304", "v2_5_1", "quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1"),
    ("PXD017199", "atlas", "quantms-collections/absolute-expression-2.0/cell-lines/PXD017199"),
    ("PXD041421", "atlas", "quantms-collections/absolute-expression-2.0/cell-lines/PXD041421"),
    *[(d.replace("-", "_"), v, f"quantmsdiann-benchmarks/phospho/{d}/{v}")
      for d in ("PXD034128-biological-study", "PXD034128-highspeed-DIA", "PXD034623", "PXD049692")
      for v in ("v2_5_1", "v2_5_1_enterprise")],
    ("PXD064049", "v2_5_1_enterprise", "quantmsdiann-benchmarks/spatial/PXD064049/v2_5_1_enterprise"),
    ("PXD071075", "v2_5_1_enterprise", "quantmsdiann-benchmarks/PXD071075-2.5.1-ent"),
]

SKIP = re.compile(r"(msstats_in\.csv|/qpx/|\.speclib$|\.quant$|/work/|/library_generation/|nf-.*-reports\.tsv$|\.predicted\.)",
                  re.I)


def classify(path: str) -> str | None:
    """PRIDE file type for a run-relative path, or None to leave it out."""
    name = path.rsplit("/", 1)[-1]
    if SKIP.search(path):
        return None
    if re.fullmatch(r"diann_report\.(parquet|tsv)", name) or name == "report.tsv":
        return "SEARCH"
    if re.search(r"(pg_matrix|pr_matrix|gg_matrix|unique_genes_matrix)\.tsv$|site_report\.parquet$|stats\.tsv$", name):
        return "RESULT"
    if name.endswith(".sdrf.tsv"):
        return "EXPERIMENTAL_DESIGN"
    if (name in ("diann_config.cfg", "diannsummary.log", "multiqc_report.html", "run_metadata.json")
            or re.match(r"params_.*\.json$", name) or "trace" in name or name == "nextflow_report.html"):
        return "OTHER"
    return None


def listing(url: str, tries: int = 5) -> list[str]:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                html = r.read().decode("utf-8", "replace")
            break
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)
    return [h for h in re.findall(r'href="([^"?#]+)"', html) if not h.startswith(("/", "http", ".."))]


def crawl(base: str, rel: str = "", depth: int = 3) -> list[str]:
    out = []
    for h in listing(f"{base}/{rel}"):
        if h.endswith("/"):
            if depth > 0 and h.rstrip("/") not in ("work", "qpx", "library_generation"):
                out += crawl(base, rel + h, depth - 1)
        else:
            out.append(rel + h)
    return out


def head_size(url: str, tries: int = 5) -> int:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.headers["Content-Length"] is not None:
                    return int(r.headers["Content-Length"])
            # Some files come back without a Content-Length: count the bytes instead.
            with urllib.request.urlopen(url, timeout=300) as r:
                return sum(len(chunk) for chunk in iter(lambda: r.read(1 << 20), b""))
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def main() -> None:
    rows = []
    for ds, ver, run in RUNS:
        base = f"{FTP}/{run}"
        try:
            files = crawl(base)
        except Exception as exc:  # missing run dir: report, keep going
            print(f"WARNING {run}: {exc}")
            continue
        for rel in files:
            ftype = classify("/" + rel)
            if ftype:
                rows.append({"dataset": ds, "diann_version": ver, "file_type": ftype,
                             "upload_name": f"{ds}_{ver}_{rel.rsplit('/', 1)[-1]}",
                             "source_url": f"{base}/{rel}"})
        print(f"{ds} {ver}: {sum(r['dataset'] == ds and r['diann_version'] == ver for r in rows)} files", flush=True)
    with ThreadPoolExecutor(4) as ex:
        for r, size in zip(rows, ex.map(lambda r: head_size(r["source_url"]), rows)):
            r["size_bytes"] = size
    # Same file name twice in one run (e.g. two traces): keep names unique.
    seen: dict[str, int] = {}
    for r in rows:
        n = seen.get(r["upload_name"], 0)
        seen[r["upload_name"]] = n + 1
        if n:
            stem, dot, ext = r["upload_name"].rpartition(".")
            r["upload_name"] = f"{stem}_{n}{dot}{ext}"
    OUT.mkdir(exist_ok=True)
    with open(OUT / "submission_files.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    tot = sum(r["size_bytes"] for r in rows)
    by = {}
    for r in rows:
        by[r["file_type"]] = by.get(r["file_type"], 0) + r["size_bytes"]
    print(f"{len(rows)} files, {tot / 1e9:.1f} GB:", {k: f"{v / 1e9:.1f} GB" for k, v in by.items()})


if __name__ == "__main__":
    main()
