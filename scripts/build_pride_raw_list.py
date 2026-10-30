"""Build the RAW-file list for the PRIDE reanalysis submission.

The submission does not re-upload any raw data. Instead it ships one text
file listing every raw file the reanalyses actually processed, each with its
public download URL in the original deposit (PRIDE FTP, MassIVE FTP, or the
ProteoBench server for Module 7).

Ground truth for "files used" is DIA-NN's own per-run stats table
(`diann_report.stats.tsv`, column File.Name) in every published run
directory, plus the run columns of the pr_matrix for the two atlas-only
cohorts whose run directories are not on the FTP (PXD017199, PXD041421).

Usage:
    python -m scripts.build_pride_raw_list --out pride_submission
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import http.client
import time
import urllib.error
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes"
PRIDE_FILES_API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}/files/all"
MASSIVE_FTP = "ftp://massive-ftp.ucsd.edu"
PROTEOBENCH_SDRF = (
    "https://raw.githubusercontent.com/bigbio/sdrf-annotated-datasets/main/"
    "datasets/ProteoBench-Astral-Module7/ProteoBench-Astral-Module7.sdrf.tsv"
)
REPO = Path(__file__).resolve().parents[1]

# Published run directories (relative to FTP) whose stats.tsv lists the runs.
RUN_DIRS = [
    "quantms-collections/absolute-expression-2.0/cell-lines/PXD003539",
    "quantmsdiann-benchmarks/PXD071075-2.5.1-ent",
    "quantmsdiann-benchmarks/cell-lines/PXD003539/v2_5_1",
    "quantmsdiann-benchmarks/cell-lines/PXD004701/v2_5_1",
    "quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1",
    *[f"quantmsdiann-benchmarks/phospho/{d}/{v}"
      for d in ("PXD034128-biological-study", "PXD034128-highspeed-DIA", "PXD034623", "PXD049692")
      for v in ("v2_5_1", "v2_5_1_enterprise")],
    *[f"quantmsdiann-benchmarks/proteobench/{d}/{v}"
      for d in ("PXD049412", "PXD062685", "PXD070049", "ProteoBench_Module_7")
      for v in ("v1_8_1", "v2_5_1", "v2_5_1_enterprise")],
    "quantmsdiann-benchmarks/single-cell/MSV000093870_plexDIA/v2_5_0",
    "quantmsdiann-benchmarks/single-cell/PXD046357/v1_8_1",
    "quantmsdiann-benchmarks/single-cell/PXD046357/v2_5_1_enterprise",
    "quantmsdiann-benchmarks/single-cell/PXD049412/v1_8_1",
    "quantmsdiann-benchmarks/single-cell/PXD049412/v2_5_1_enterprise",
    "quantmsdiann-benchmarks/spatial/PXD064049/v2_5_1_enterprise",
]
# Atlas-only cohorts: run names from the local pr_matrix header.
LOCAL_PR_MATRICES = {
    "PXD017199": REPO / "data/PXD017199/diann_report.pr_matrix.tsv",
    "PXD041421": REPO / "data/PXD041421/diann_report.pr_matrix.tsv",
}
PR_MATRIX_ANNOTATION_COLS = {
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
    "Proteotypic", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
    "Precursor.Id", "Precursor.Mz",
}
# Converted/intermediate suffixes stripped to recover the original raw stem.
SUFFIXES = (".mzml", ".raw", ".wiff2", ".wiff", ".scan", ".dia", ".zip", ".tar", ".gz", ".d")


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def stem(name: str) -> str:
    name = re.split(r"[\\/]", name.rstrip("/\\"))[-1]
    name = re.sub(r"\s+", "", name)  # ProCan zips hold names like "P02_ 8064_2 _S"
    changed = True
    while changed:
        changed = False
        for s in SUFFIXES:
            if name.lower().endswith(s):
                name, changed = name[: -len(s)], True
    return name


def accession_of(run_dir: str) -> str:
    if "ProteoBench_Module_7" in run_dir:
        return "ProteoBench_Module_7"
    return re.search(r"(PXD\d{6}|MSV\d{9})", run_dir).group(1)


def runs_from_stats(run_dir: str) -> list[str]:
    text = fetch(f"{FTP}/{run_dir}/quant_tables/diann_report.stats.tsv").decode()
    return [row["File.Name"] for row in csv.DictReader(text.splitlines(), delimiter="\t")]


def runs_from_pr_matrix(path: Path) -> list[str]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    return [c for c in header if c not in PR_MATRIX_ANNOTATION_COLS]


class HttpRangeFile(io.RawIOBase):
    """Seekable read-only view of a remote file via HTTP Range requests, so
    zipfile can read an archive's central directory without downloading it."""

    def __init__(self, url: str):
        self.url, self.pos = url, 0
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Length"])

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def read(self, n: int = -1) -> bytes:
        end = self.size if n is None or n < 0 else min(self.size, self.pos + n)
        if end <= self.pos:
            return b""
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end - 1}"})
        for attempt in range(5):  # the PRIDE FTP drops range requests now and then
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except (http.client.IncompleteRead, urllib.error.URLError):
                if attempt == 4:
                    raise
                time.sleep(5 * (attempt + 1))
        self.pos += len(data)
        return data


def zip_members(url: str) -> list[tuple[str, int]]:
    with zipfile.ZipFile(HttpRangeFile(url)) as z:
        return [(i.filename, i.file_size) for i in z.infolist()]


def pride_raw_files(acc: str, wanted: set[str]) -> dict[str, list[dict]]:
    """stem -> PRIDE RAW file records (a .wiff stem maps to .wiff + .wiff.scan).

    Some deposits pack runs into archives (ProCan per-day zips, timsTOF .d
    folders zipped by experiment). For runs not deposited as individual files,
    the archive's member list is read remotely and the run is reported as
    `<archive URL>` + `archive_member`."""
    by_stem: dict[str, list[dict]] = defaultdict(list)
    archives = []
    for f in json.loads(fetch(PRIDE_FILES_API.format(acc=acc))):
        if f.get("fileCategory", {}).get("value") != "RAW":
            continue
        ftp = next((l["value"] for l in f.get("publicFileLocations", [])
                    if l.get("name") == "FTP Protocol"), "")
        rec = {"file_name": f["fileName"], "url": ftp, "archive_member": "",
               "size_bytes": f.get("fileSizeBytes", ""), "checksum": f.get("checksum", "")}
        by_stem[stem(f["fileName"])].append(rec)
        if f["fileName"].lower().endswith(".zip"):
            archives.append(rec)
    direct = set(by_stem)
    taken: set[tuple[str, str]] = set()  # (stem, member name): first archive wins
    unreadable: list[dict] = []
    if wanted - direct:
        for a in archives:
            if stem(a["file_name"]) in wanted:
                continue  # the archive itself is a run (e.g. X.d.zip)
            seen: set[str] = set()
            try:
                members = zip_members(a["url"].replace("ftp://", "https://", 1))
            except zipfile.BadZipFile as e:
                print(f"  {acc}: cannot list {a['file_name']}: {e}", file=sys.stderr)
                unreadable.append(a)
                continue
            for member, size in members:
                if member.startswith("__MACOSX/"):
                    continue
                parts = member.rstrip("/").split("/")
                # A run is a file (X.raw, X.wiff, X.wiff.scan) or a folder (X.d/...).
                for i, part in enumerate(parts):
                    s = stem(part)
                    if s not in wanted or s in direct:
                        continue
                    top = "/".join(parts[: i + 1])
                    is_dir_run = i < len(parts) - 1 or member.endswith("/")
                    if top not in seen and (s, part) not in taken:
                        seen.add(top)
                        taken.add((s, part))
                        by_stem[s].append({"file_name": part, "url": a["url"],
                                           "archive_member": top + ("/" if is_dir_run else ""),
                                           "size_bytes": "" if is_dir_run else size,
                                           "checksum": ""})
                    break
            print(f"  {acc}: listed {a['file_name']}", file=sys.stderr)
    # Runs inside an unreadable archive (PXD030304 181231.zip is a 4 KB stub on the PRIDE
    # FTP) are not listed; they end up in raw_files_unmatched.tsv.
    for a in unreadable:
        print(f"  {acc}: skipped unreadable archive {a['file_name']}", file=sys.stderr)
    return by_stem


def massive_raw_files(acc: str) -> dict[str, list[dict]]:
    """MassIVE: list the raw/ tree over FTP (only MSV000093870 here)."""
    from ftplib import FTP_TLS, error_perm, error_temp

    by_stem: dict[str, list[dict]] = defaultdict(list)
    ftp = FTP_TLS("massive-ftp.ucsd.edu")  # MassIVE requires TLS
    ftp.login()
    ftp.prot_p()
    stack = []
    for vol in ftp.nlst():  # the dataset may be split across volumes; raw/ is on one
        try:
            if f"{vol}/{acc}/raw" in ftp.nlst(f"{vol}/{acc}"):
                stack.append(f"{vol}/{acc}/raw")
        except (error_perm, error_temp):
            continue
    while stack:
        d = stack.pop()
        for entry, facts in ftp.mlsd(d):
            path = f"{d}/{entry}"
            if facts.get("type") == "dir":
                stack.append(path)
            elif entry.lower().endswith((".raw", ".d", ".wiff", ".wiff.scan", ".mzml", ".zip")):
                by_stem[stem(entry)].append({
                    "file_name": entry, "url": f"{MASSIVE_FTP}/{path}", "archive_member": "",
                    "size_bytes": facts.get("size", ""), "checksum": "",
                })
    ftp.quit()
    return by_stem


def proteobench_raw_files() -> dict[str, list[dict]]:
    text = fetch(PROTEOBENCH_SDRF).decode()
    by_stem: dict[str, list[dict]] = defaultdict(list)
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        name, uri = row["comment[data file]"], row["comment[file uri]"]
        if not any(r["file_name"] == name for r in by_stem[stem(name)]):
            by_stem[stem(name)].append({"file_name": name, "url": uri, "archive_member": "",
                                        "size_bytes": "", "checksum": ""})
    return by_stem


def sdrf_uri_fallback(acc: str) -> dict[str, list[dict]]:
    """Runs inside archives we cannot list remotely (e.g. .7z): take the archive
    URL from the curated SDRF's comment[file uri]."""
    url = ("https://raw.githubusercontent.com/bigbio/sdrf-annotated-datasets/main/"
           f"datasets/{acc}/{acc}.sdrf.tsv")
    by_stem: dict[str, list[dict]] = defaultdict(list)
    try:
        text = fetch(url).decode()
    except OSError:
        return by_stem
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        name, uri = row.get("comment[data file]", ""), row.get("comment[file uri]", "")
        if name and uri.startswith(("http", "ftp")) and not by_stem[stem(name)]:
            member = name + ("/" if name.lower().endswith(".d") else "")
            by_stem[stem(name)].append({"file_name": name, "url": uri,
                                        "archive_member": member if stem(uri) != stem(name) else "",
                                        "size_bytes": "", "checksum": ""})
    return by_stem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=REPO / "pride_submission")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # accession -> stem -> set of run directories that processed it
    used: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for rd in RUN_DIRS:
        acc = accession_of(rd)
        for run in runs_from_stats(rd):
            used[acc][stem(run)].add(rd)
        print(f"{rd}: ok", file=sys.stderr)
    for acc, path in LOCAL_PR_MATRICES.items():
        for run in runs_from_pr_matrix(path):
            used[acc][stem(run)].add(f"atlas:{acc}")

    rows, missing = [], []
    for acc in sorted(used):
        if acc.startswith("PXD"):
            index, repo = pride_raw_files(acc, set(used[acc])), "PRIDE"
        elif acc.startswith("MSV"):
            index, repo = massive_raw_files(acc), "MassIVE"
        else:
            index, repo = proteobench_raw_files(), "ProteoBench"
        if repo == "PRIDE" and set(used[acc]) - set(index):
            for k, v in sdrf_uri_fallback(acc).items():
                index.setdefault(k, v)
        for s, dirs in sorted(used[acc].items()):
            hits = index.get(s)
            if not hits:
                missing.append((acc, s))
                continue
            for h in hits:
                rows.append({"dataset": acc, "repository": repo, **h,
                             "used_in": ";".join(sorted(dirs))})
        print(f"{acc}: {len(used[acc])} runs", file=sys.stderr)

    cols = ["dataset", "repository", "file_name", "url", "archive_member", "size_bytes", "checksum", "used_in", "note"]
    with (args.out / "raw_files.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with (args.out / "raw_files.txt").open("w") as fh:
        fh.write("# Raw files reanalysed with quantmsdiann. They are NOT re-uploaded here:\n")
        fh.write("# each is listed with its public URL in the original deposit. Where the\n")
        fh.write("# deposit packs runs into an archive, the third column gives the path\n")
        fh.write("# inside that archive. Full table (sizes, checksums, runs) in raw_files.tsv.\n")
        fh.write("# raw_file\turl\tarchive_member\n")
        for r in rows:
            fh.write(f"{r['file_name']}\t{r['url']}\t{r['archive_member']}\n")
    (args.out / "raw_files_unmatched.tsv").unlink(missing_ok=True)
    if missing:
        with (args.out / "raw_files_unmatched.tsv").open("w") as fh:
            fh.write("dataset\trun_stem\n")
            fh.writelines(f"{a}\t{s}\n" for a, s in missing)
    print(f"{len(rows)} files listed, {len(missing)} runs unmatched", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
