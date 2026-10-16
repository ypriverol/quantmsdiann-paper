"""Stage every file of the PRIDE submission into one upload folder.

Run where there is space (the upload is ~82 GB), e.g. the EMBL-EBI cluster:

    python -m scripts.stage_pride_submission --dest /hps/.../quantmsdiann_pride_upload

FTP files are downloaded (resumable, size-checked) under their upload names;
local files (SDRFs, paired standalone runs, raw-file list) are copied from this
repo; each <PXD>.qpx.zip is built from that run's qpx/ folder on the FTP, and each
<dataset>_<version>.zip run bundle from the members listed in zip_members.tsv.
--reuse <dir> takes files from an earlier staging folder instead of downloading
them again (hard links when possible).
Run it from a checkout that has pride_submission/ and
data/proteobench-review-experiments/x/ (or rsync those two folders first).
submission.px is copied last. Upload the folder with the PRIDE submission tool
or Aspera, pointing it at submission.px.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd

from scripts.build_submission_manifest import head_size, listing

REPO = Path(__file__).resolve().parents[1]
SUB = REPO / "pride_submission"


def qpx_zip(folder_url: str, out: Path, tmp: Path) -> None:
    """Download every file of an FTP qpx/ folder and store them (no recompression: parquet) in out."""
    tmp.mkdir(parents=True, exist_ok=True)
    names = [n for n in listing(folder_url) if not n.endswith("/")]
    for n in names:
        subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(tmp / n), folder_url + n], check=True)
        if (tmp / n).stat().st_size != head_size(folder_url + n):
            raise RuntimeError(f"size mismatch: {folder_url + n}")
    part = out.with_suffix(".part")
    with zipfile.ZipFile(part, "w", zipfile.ZIP_STORED, allowZip64=True) as z:
        for n in names:
            z.write(tmp / n, n)
    part.rename(out)
    shutil.rmtree(tmp)


def fetch(src: str, size: int, out: Path, reuse: Path | None, old_name: str) -> Path:
    """Put one file at out: from the reuse folder if it has it at the right size, else curl/copy."""
    prev = reuse / old_name if reuse and old_name else None
    if prev is not None and prev.exists() and prev.stat().st_size == size:
        return prev
    if src.startswith("http"):
        subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(out), src], check=True)
    else:
        shutil.copy2(src, out)
    if out.stat().st_size != size:
        raise RuntimeError(f"size mismatch: {src}")
    return out


def run_zip(name: str, members: pd.DataFrame, out: Path, tmp: Path, reuse: Path | None) -> None:
    """Build one run bundle: members stored (no recompression) under <run>/<member path>."""
    tmp.mkdir(parents=True, exist_ok=True)
    top = name.removesuffix(".zip")
    part = out.with_suffix(".part")
    with zipfile.ZipFile(part, "w", zipfile.ZIP_STORED, allowZip64=True) as z:
        for m in members.itertuples():
            path = fetch(m.source, int(m.size_bytes), tmp / "member", reuse,
                         m.old_upload_name if isinstance(m.old_upload_name, str) else "")
            z.write(path, f"{top}/{m.member}")
            if path.parent == tmp:
                path.unlink()
    part.rename(out)
    shutil.rmtree(tmp)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True, type=Path)
    ap.add_argument("--reuse", type=Path, help="earlier staging folder to take files from")
    args = ap.parse_args()
    dest, reuse = args.dest, args.reuse
    dest.mkdir(parents=True, exist_ok=True)
    files = pd.read_csv(SUB / "submission_files_final.tsv", sep="\t")
    members = pd.read_csv(SUB / "zip_members.tsv", sep="\t") if (SUB / "zip_members.tsv").exists() else None
    bad = []
    for i, r in enumerate(files.itertuples(), 1):
        out = dest / r.upload_name
        if str(r.source_url).startswith("zip:"):  # run bundle built from zip_members.tsv
            if not out.exists():
                run_zip(r.upload_name, members[members["zip"] == r.upload_name], out, dest / ".zip_tmp", reuse)
            print(f"[{i}/{len(files)}] {r.upload_name}", flush=True)
            continue
        if reuse and (reuse / r.upload_name).exists() and not out.exists() and (
                r.upload_name.endswith(".qpx.zip") or (reuse / r.upload_name).stat().st_size == r.size_bytes):
            try:
                os.link(reuse / r.upload_name, out)
            except OSError:
                shutil.copy2(reuse / r.upload_name, out)
        if r.upload_name.endswith(".qpx.zip"):  # source is a qpx/ folder: download its files, zip them
            if not out.exists():
                qpx_zip(r.source_url, out, dest / ".qpx_tmp")
            print(f"[{i}/{len(files)}] {r.upload_name}", flush=True)
            continue
        if out.exists() and out.stat().st_size == r.size_bytes:
            continue
        if r.source_url.startswith("http"):
            subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(out), r.source_url], check=True)
        else:
            shutil.copy2(r.source_url, out)
        if out.stat().st_size != r.size_bytes:
            bad.append(r.upload_name)
        print(f"[{i}/{len(files)}] {r.upload_name}", flush=True)
    shutil.copy2(SUB / "submission.px", dest / "submission.px")
    print("size mismatches:", bad or "none")


if __name__ == "__main__":
    main()
