"""Create (or resume) the Zenodo data deposit and upload pride_submission/zenodo_files.tsv.

Files with an FTP URL are streamed straight from PRIDE to Zenodo (nothing is
written to local disk); local files are read from disk. Already-uploaded files
with the same size are skipped, so the script can be re-run after a failure.
The deposit is left as a DRAFT: review it on zenodo.org and press Publish.

    export ZENODO_TOKEN=...            # zenodo.org > Applications > Personal access token (deposit:write)
    python -m scripts.zenodo_upload                 # creates a new draft
    python -m scripts.zenodo_upload --id 1234567    # resume an existing draft
    python -m scripts.zenodo_upload --sandbox       # try on sandbox.zenodo.org first
"""
from __future__ import annotations

import argparse
import json
import os
import urllib.request
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
SUB = REPO / "pride_submission"


class Sized:
    """File-like wrapper with a known length, so requests sends Content-Length."""

    def __init__(self, fh, size: int):
        self.fh, self.size = fh, size

    def __len__(self):
        return self.size

    def read(self, n=-1):
        return self.fh.read(n)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", help="existing draft deposition id")
    ap.add_argument("--sandbox", action="store_true")
    a = ap.parse_args()
    api = "https://sandbox.zenodo.org/api" if a.sandbox else "https://zenodo.org/api"
    s = requests.Session()
    s.headers["Authorization"] = f"Bearer {os.environ['ZENODO_TOKEN']}"
    meta = json.loads((REPO / "docs/zenodo_data_metadata.json").read_text())
    if a.id:
        dep = s.get(f"{api}/deposit/depositions/{a.id}").json()
    else:
        dep = s.post(f"{api}/deposit/depositions", json={}).json()
        s.put(f"{api}/deposit/depositions/{dep['id']}", json=meta).raise_for_status()
    print(f"deposition {dep['id']}: {dep['links']['html']}")
    bucket = dep["links"]["bucket"]
    done = {f["filename"]: f["filesize"] for f in s.get(f"{api}/deposit/depositions/{dep['id']}/files").json()}
    files = pd.read_csv(SUB / "zenodo_files.tsv", sep="\t")
    for i, r in enumerate(files.itertuples(), 1):
        if done.get(r.upload_name) == r.size_bytes:
            continue
        src = r.source_url
        fh = urllib.request.urlopen(src, timeout=600) if src.startswith("http") else open(src, "rb")
        with fh:
            resp = s.put(f"{bucket}/{r.upload_name}", data=Sized(fh, int(r.size_bytes)), timeout=None)
        resp.raise_for_status()
        print(f"[{i}/{len(files)}] {r.upload_name} ({r.size_bytes / 1e6:.0f} MB)", flush=True)
    print("All files uploaded. Review the draft and publish it on Zenodo.")


if __name__ == "__main__":
    main()
