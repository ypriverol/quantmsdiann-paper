"""Stage every file of the PRIDE submission into one upload folder.

Run where there is space (the upload is ~64 GB), e.g. the EMBL-EBI cluster:

    python -m scripts.stage_pride_submission --dest /hps/.../quantmsdiann_pride_upload

FTP files are downloaded (resumable, size-checked) under their upload names;
local files (corrected SDRFs, paired standalone runs, raw-file list) are copied
from this repo, so run it from a checkout that has pride_submission/ and
data/proteobench-review-experiments/x/ (or rsync those two folders first).
submission.px is copied last. Upload the folder with the PRIDE submission tool
or Aspera, pointing it at submission.px.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SUB = REPO / "pride_submission"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True, type=Path)
    dest = ap.parse_args().dest
    dest.mkdir(parents=True, exist_ok=True)
    files = pd.read_csv(SUB / "submission_files_final.tsv", sep="\t")
    bad = []
    for i, r in enumerate(files.itertuples(), 1):
        out = dest / r.upload_name
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
