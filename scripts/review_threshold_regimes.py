"""Reviewer response (R1 Q1.10, Q1.11, Q1.16, Q1.21): counts under several
evidence thresholds for the three reanalysis-recovery cohorts of Fig. 4.

For NCI-60 (PXD003539), Sun (PXD004701) and ProCan (PXD030304) this counts,
for the deposited original and the quantmsdiann reanalysis:
  - proteins / protein groups supported by >= 1, 2, 3, 5 peptides;
  - precursors and peptide sequences, where the deposit allows it (Q1.11);
  - the peptides-per-protein distribution (Q1.21).

It also makes the Q1.16 discrepancy explicit. The NCI-60 OpenSWATH matrix keys
rows by transition group (modified sequence + charge), so the Fig. 4 baseline
(4,284 "proteins with >= 2 peptides") counts distinct precursors per protein,
while the removed Fig. S3 counted proteotypic stripped peptides at score <= 0.01
(3,644). Both definitions are reported side by side.

Original-side sources are the same deposited files Fig. 4 uses. quantmsdiann
side: the DIA-NN report parquets used by Fig. 4, counted under methods.md
section 1 (Lib.Q.Value / Lib.PG.Q.Value <= 0.01, decoys dropped, nothing
else). Reports are streamed from the PRIDE FTP and never written to disk.

Usage:
    python -m scripts.review_threshold_regimes
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from scripts import rebuild as rb

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
OUT = REPO / "analysis/figures/review/threshold_regimes"
FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantmsdiann-benchmarks"
THRESHOLDS = (1, 2, 3, 5)
PR_META = {"Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
           "Proteotypic", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id"}
COHORTS = {
    "PXD003539": {"label": "NCI-60 (Guo 2019)", "original": "OpenSWATH (deposited)",
                  "report": rb.DIANN_REPORT_PARQUET_URL},
    "PXD004701": {"label": "Breast cancer lines (Sun 2023)", "original": "OpenSWATH (deposited)",
                  "report": f"{FTP}/cell-lines/PXD004701/v2_5_1/quant_tables/diann_report.parquet"},
    "PXD030304": {"label": "ProCan-DepMapSanger", "original": "DIA-NN 1.8 (figshare)",
                  "report": f"{FTP}/cell-lines/PXD030304/v2_5_1/quant_tables/diann_report.parquet"},
}
OPENSWATH_ID = re.compile(r"^\d+_(?P<modseq>.+)_(?P<charge>\d+)_run\d+$")
UNIMOD = re.compile(r"\(UniMod:\d+\)")


def strip(modseq: str) -> str:
    return UNIMOD.sub("", modseq)


# ---------------------------------------------------------------- originals

def nci60_original() -> dict:
    """Deposited OpenSWATH matrix; quantified target rows (>= 1 Intensity_*)."""
    path = DATA / "PXD003539/feature_alignment_requant_matrix.tsv"
    header = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    inten = [c for c in header if c.startswith("Intensity_")]
    score = [c for c in header if c.startswith("score_")]
    prec_per_prot, pep_per_prot, s3_pep_per_prot = defaultdict(set), defaultdict(set), defaultdict(set)
    precursors, peptides = set(), set()
    for ch in pd.read_csv(path, sep="\t", usecols=["Peptide", "Protein"] + inten + score, dtype=str, chunksize=50000):
        decoy = ch["Peptide"].str.startswith("DECOY_", na=False) | ch["Protein"].str.contains("DECOY", case=False, na=False)
        ch = ch[~decoy]
        quant = ch[ch[inten].notna().any(axis=1)]
        for prot, pid in zip(quant["Protein"], quant["Peptide"]):
            m = OPENSWATH_ID.match(str(pid))
            if not m or not isinstance(prot, str):
                continue
            prec = f"{m['modseq']}{m['charge']}"
            prec_per_prot[prot].add(prec)
            pep_per_prot[prot].add(strip(m["modseq"]))
            precursors.add(prec)
            peptides.add(strip(m["modseq"]))
        # Removed Fig. S3 rule: proteotypic (1/...) rows detected at score <= 0.01.
        prot1 = ch[ch["Protein"].str.startswith("1/", na=False)]
        det = prot1[(prot1[score].apply(pd.to_numeric, errors="coerce") <= 0.01).any(axis=1)]
        for prot, pid in zip(det["Protein"], det["Peptide"]):
            m = OPENSWATH_ID.match(str(pid))
            if m:
                s3_pep_per_prot[prot].add(strip(m["modseq"]))
    return {"pep_per_prot": {p: len(s) for p, s in pep_per_prot.items()},
            "prec_per_prot": {p: len(s) for p, s in prec_per_prot.items()},
            "s3_pep_per_prot": {p: len(s) for p, s in s3_pep_per_prot.items()},
            "precursors": len(precursors), "peptides": len(peptides)}


def sun_original() -> dict:
    """Deposited OpenSWATH ids at mscore <= 0.01, proteotypic SwissProt peptides
    mapped through the deposited library (the Fig. 4 baseline rule)."""
    lib = pd.read_csv(DATA / "PXD004701/sun_library.tsv", sep="\t",
                      usecols=["PeptideSequence", "ProteinName", "decoy"], dtype=str)
    lib = lib[lib["decoy"].astype(str) == "0"]
    pep_acc = {}
    for pep, pn in zip(lib["PeptideSequence"], lib["ProteinName"]):
        m = re.match(r"(\d+)/(.*)", str(pn))
        if m and int(m.group(1)) == 1 and m.group(2).startswith("sp|"):
            pep_acc[pep] = m.group(2).split("|")[1]
    pep_per_prot, precursors, peptides = defaultdict(set), set(), set()
    for ch in pd.read_csv(DATA / "PXD004701/1R_all_id.txt", sep="\t", usecols=["tg", "mscore"],
                          dtype=str, chunksize=3_000_000):
        ch = ch[pd.to_numeric(ch["mscore"], errors="coerce") <= 0.01]
        for tg in ch["tg"].unique():
            m = rb._SUN_TG_RE.match(str(tg))
            if not m:
                continue
            pep, charge = m.group(1), m.group(2)
            precursors.add(f"{pep}{charge}")
            peptides.add(strip(pep))
            if pep in pep_acc:
                pep_per_prot[pep_acc[pep]].add(strip(pep))
    return {"pep_per_prot": {p: len(s) for p, s in pep_per_prot.items()},
            "precursors": len(precursors), "peptides": len(peptides)}


def procan_original() -> dict:
    """figshare peptide-count matrix: peptides per protein per sample. The union
    across samples is not recoverable, so the per-protein maximum is used (a
    lower bound on distinct peptides). No precursor/peptide identities."""
    mat = pd.read_csv(DATA / "PXD030304/peptide_counts_per_protein_per_sample.txt", sep="\t", index_col=0)
    mx = mat.apply(pd.to_numeric, errors="coerce").max(axis=0).dropna()
    return {"pep_per_prot": {p: int(v) for p, v in mx.items() if v >= 1},
            "precursors": np.nan, "peptides": np.nan}


# ---------------------------------------------------------------- quantmsdiann

def quantmsdiann_side(url: str) -> dict:
    """Stream the DIA-NN report (same file as Fig. 4) under methods.md section 1:
    precursors / peptides at Lib.Q.Value <= 0.01; protein groups at
    Lib.PG.Q.Value <= 0.01, each supported by its Lib.Q.Value <= 0.01 peptides
    (the rule behind Fig. 4's prot_2pep). Decoys dropped, nothing else."""
    import fsspec
    import pyarrow.parquet as pq
    cols = ["Protein.Group", "Proteotypic", "Stripped.Sequence", "Precursor.Id", "Lib.Q.Value", "Lib.PG.Q.Value", "Decoy"]
    pep_per_pg, pep_per_pg_pt = defaultdict(set), defaultdict(set)
    pgs, precursors, peptides = set(), set(), set()
    with fsspec.filesystem("https").open(url, "rb", block_size=64 * 2**20) as fh:
        pf = pq.ParquetFile(fh)
        for batch in pf.iter_batches(batch_size=2_000_000, columns=cols):
            d = batch.to_pandas()
            d = d[d["Decoy"] == 0]
            pgs.update(d.loc[d["Lib.PG.Q.Value"] <= 0.01, "Protein.Group"])
            p = d[d["Lib.Q.Value"] <= 0.01]
            precursors.update(p["Precursor.Id"])
            peptides.update(p["Stripped.Sequence"])
            p = p[p["Lib.PG.Q.Value"] <= 0.01].drop_duplicates(["Protein.Group", "Stripped.Sequence"])
            for pg, seq, pt in zip(p["Protein.Group"], p["Stripped.Sequence"], p["Proteotypic"]):
                pep_per_pg[pg].add(seq)
                if pt == 1:
                    pep_per_pg_pt[pg].add(seq)
    # Protein groups passing Lib.PG.Q.Value with no Lib.Q-passing peptide still count at >= 1.
    counts = {pg: len(pep_per_pg.get(pg, ())) for pg in pgs}
    return {"pep_per_prot": {pg: max(n, 1) for pg, n in counts.items()},
            "pep_per_prot_proteotypic": {p: len(s) for p, s in pep_per_pg_pt.items()},
            "precursors": len(precursors), "peptides": len(peptides)}


def at_least(counts: dict, k: int) -> int:
    return sum(1 for v in counts.values() if v >= k)


# ---------------------------------------------------------------- main

def main() -> None:
    rb.apply_house_style()
    (OUT / "data").mkdir(parents=True, exist_ok=True)
    originals = {"PXD003539": nci60_original, "PXD004701": sun_original, "PXD030304": procan_original}
    rows, dists = [], {}
    for ds, meta in COHORTS.items():
        print(f"{ds}: original", flush=True)
        o = originals[ds]()
        print(f"{ds}: quantmsdiann", flush=True)
        q = quantmsdiann_side(meta["report"])
        dists[ds] = (o["pep_per_prot"], q["pep_per_prot"])
        for side, d, rule in (("original", o["pep_per_prot"], f"{meta['original']}: distinct peptide sequences per protein"),
                              ("quantmsdiann", q["pep_per_prot"], "DIA-NN report: Lib.PG.Q <= 1% groups, distinct Lib.Q <= 1% peptides"),
                              ("quantmsdiann", q["pep_per_prot_proteotypic"], "DIA-NN report: proteotypic peptides only")):
            for k in THRESHOLDS:
                rows.append({"dataset": ds, "side": side, "rule": rule, "metric": f"proteins >= {k} peptides",
                             "value": at_least(d, k)})
        if ds == "PXD003539":
            for rule, d in (("OpenSWATH: distinct precursors per protein (Fig. 4 baseline)", o["prec_per_prot"]),
                            ("OpenSWATH: proteotypic stripped peptides, score <= 0.01 (removed Fig. S3)", o["s3_pep_per_prot"])):
                for k in THRESHOLDS:
                    rows.append({"dataset": ds, "side": "original", "rule": rule, "metric": f"proteins >= {k} peptides",
                                 "value": at_least(d, k)})
        for side, d in (("original", o), ("quantmsdiann", q)):
            for metric in ("precursors", "peptides"):
                rows.append({"dataset": ds, "side": side, "rule": "global 1% (deposit's own FDR for originals)",
                             "metric": metric, "value": d[metric]})
            pp = pd.Series(d["pep_per_prot"])
            rows.append({"dataset": ds, "side": side, "rule": "", "metric": "median peptides per protein",
                         "value": float(pp.median())})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "data/threshold_regimes.tsv", sep="\t", index=False)

    # Figure: per cohort, proteins at >= k peptides (a) and peptides-per-protein distribution (b).
    fig, axes = plt.subplots(2, 3, figsize=(22, 12), gridspec_kw={"hspace": 0.55, "wspace": 0.3})
    colors = {"original": rb.COMPARISON["original"], "quantmsdiann": rb.COMPARISON["quantmsdiann"]}
    for c, (ds, meta) in enumerate(COHORTS.items()):
        ax = axes[0][c]
        sub = df[(df["dataset"] == ds) & df["metric"].str.startswith("proteins >=")]
        width = 0.38
        for j, side in enumerate(("original", "quantmsdiann")):
            s = sub[sub["side"] == side]
            rule = (f"{meta['original']}: distinct peptide sequences per protein" if side == "original"
                    else "DIA-NN report: Lib.PG.Q <= 1% groups, distinct Lib.Q <= 1% peptides")
            s = s[s["rule"] == rule]
            vals = [s.loc[s["metric"] == f"proteins >= {k} peptides", "value"].iloc[0] for k in THRESHOLDS]
            x = np.arange(len(THRESHOLDS)) + (j - 0.5) * width
            ax.bar(x, vals, width=width * 0.92, color=colors[side],
                   label=meta["original"] if side == "original" else "quantmsdiann (DIA-NN 2.5.1)")
            for xi, v in zip(x, vals):
                ax.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=10, rotation=90)
        ax.set_xticks(range(len(THRESHOLDS)), [f"≥{k}" for k in THRESHOLDS])
        ax.set_xlabel("Peptides per protein")
        ax.set_ylabel("Proteins / protein groups")
        ax.set_title(f"{'abc'[c]}  {meta['label']}", loc="left")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
        ax.legend(loc="upper right", fontsize=12)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(rb.kfmt))

        ax = axes[1][c]
        o, q = dists[ds]
        bins = np.arange(1, 32)
        for side, d in (("original", o), ("quantmsdiann", q)):
            v = np.clip(np.array(list(d.values())), 1, 30)
            ax.hist(v, bins=bins, histtype="step", lw=2, color=colors[side], density=True,
                    label=f"{'original' if side == 'original' else 'quantmsdiann'} (median {np.median(list(d.values())):.0f})")
        ax.set_xlabel("Peptides per protein (30 = 30 or more)")
        ax.set_ylabel("Fraction of proteins")
        ax.set_title(f"{'def'[c]}  Peptides per protein", loc="left")
        ax.legend(loc="upper right", fontsize=12)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"supp_threshold_regimes.{ext}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
