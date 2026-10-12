"""Reviewer response (R1 Q1.12): peptide-sequence overlap, deposited original vs
quantmsdiann, for the two OpenSWATH cohorts (NCI-60 PXD003539, Sun PXD004701).

Original: distinct stripped sequences of quantified target rows (NCI-60 matrix)
or of transition groups at mscore <= 0.01 (Sun ids). quantmsdiann: distinct
Stripped.Sequence at Lib.Q.Value <= 0.01 in the Fig. 4 report, decoys dropped.
I/L are not merged.

Usage: python -m scripts.review_peptide_overlap
"""
from __future__ import annotations

from pathlib import Path

import fsspec
import pandas as pd
import pyarrow.parquet as pq

from scripts import rebuild as rb
from scripts.review_threshold_regimes import COHORTS, DATA, OPENSWATH_ID, strip

OUT = Path(__file__).resolve().parents[1] / "analysis/figures/review/threshold_regimes/data/peptide_overlap.tsv"


def nci60_peptides() -> set:
    path = DATA / "PXD003539/feature_alignment_requant_matrix.tsv"
    inten = [c for c in pd.read_csv(path, sep="\t", nrows=0).columns if c.startswith("Intensity_")]
    peps = set()
    for ch in pd.read_csv(path, sep="\t", usecols=["Peptide", "Protein"] + inten, dtype=str, chunksize=50000):
        ch = ch[~(ch["Peptide"].str.startswith("DECOY_", na=False) | ch["Protein"].str.contains("DECOY", case=False, na=False))]
        for pid in ch.loc[ch[inten].notna().any(axis=1), "Peptide"]:
            m = OPENSWATH_ID.match(str(pid))
            if m:
                peps.add(strip(m["modseq"]))
    return peps


def sun_peptides() -> set:
    peps = set()
    for ch in pd.read_csv(DATA / "PXD004701/1R_all_id.txt", sep="\t", usecols=["tg", "mscore"], dtype=str, chunksize=3_000_000):
        for tg in ch.loc[pd.to_numeric(ch["mscore"], errors="coerce") <= 0.01, "tg"].unique():
            m = rb._SUN_TG_RE.match(str(tg))
            if m:
                peps.add(strip(m.group(1)))
    return peps


def report_peptides(url: str) -> set:
    peps = set()
    with fsspec.filesystem("https").open(url, "rb", block_size=64 * 2**20) as fh:
        for b in pq.ParquetFile(fh).iter_batches(batch_size=2_000_000, columns=["Stripped.Sequence", "Lib.Q.Value", "Decoy"]):
            d = b.to_pandas()
            peps.update(d.loc[(d["Decoy"] == 0) & (d["Lib.Q.Value"] <= 0.01), "Stripped.Sequence"])
    return peps


def main() -> None:
    rows = []
    for ds, orig in (("PXD003539", nci60_peptides), ("PXD004701", sun_peptides)):
        o, q = orig(), report_peptides(COHORTS[ds]["report"])
        rows.append({"dataset": ds, "original": len(o), "quantmsdiann": len(q), "shared": len(o & q),
                     "original_only": len(o - q), "quantmsdiann_only": len(q - o),
                     "pct_original_recovered": round(len(o & q) / len(o) * 100, 1)})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)


if __name__ == "__main__":
    main()
