"""Reviewer response (R1 Q1.24): wall-clock against number of files and raw size.

Fig. 2b shows per-dataset wall-clock annotated by sample type and instrument.
The reviewer asked to contextualise it by workload, so this plots the same
wall-clock values (analysis/figures/performance/data/parallelism_data.tsv, the
Fig. 2b audit table, with the same PXD003539 exclusion) against the number of
raw files and the total raw-data volume of the files each run processed.

Raw sizes come from pride_submission/raw_files.tsv (built by
scripts.build_pride_raw_list). Gaps there are filled as follows:
  - timsTOF .d folders packed inside a deposit zip: sum of the zip members
    under that folder (read remotely from the zip's central directory);
  - PXD062685 (all runs inside one .7z): archive size split evenly per run
    (flagged approximate);
  - ProteoBench Module 7: HTTP Content-Length of each raw file.

Usage:
    python -m scripts.review_runtime_vs_scale
"""

from __future__ import annotations

import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from scripts import rebuild as rb
from scripts.build_pride_raw_list import zip_members

REPO = Path(__file__).resolve().parents[1]
PERF = REPO / "analysis/figures/performance/data/parallelism_data.tsv"
RAW = REPO / "pride_submission/raw_files.tsv"
OUT = REPO / "analysis/figures/review/runtime_vs_scale"

# Which reanalysis (`used_in` in raw_files.tsv) each Fig. 2b trace corresponds to.
USED_IN = {
    "PXD049412": "proteobench/PXD049412/",
    "PXD062685": "proteobench/PXD062685/",
    "PXD070049": "proteobench/PXD070049/",
    "ProteoBench_Module_7": "proteobench/ProteoBench_Module_7/",
    "PXD004701": "cell-lines/PXD004701/",
    "PXD030304": "cell-lines/PXD030304/",
    "MSV000093870": "MSV000093870",
    "PXD034128": "PXD034128-biological-study/",
    "PXD034623": "phospho/PXD034623/",
    "PXD064049": "spatial/PXD064049/",
    "PXD049692": "phospho/PXD049692/",
    "PXD046357": "single-cell/PXD046357/",
    "PXD017199": "atlas:PXD017199",
    "PXD041421": "atlas:PXD041421",
    "PXD071075": "PXD071075",
}
LABEL = {"ProteoBench_Module_7": "Module 7", "MSV000093870": "MSV000093870"}
CATEGORY = {
    "PXD049412": "Benchmark (ProteoBench)", "PXD062685": "Benchmark (ProteoBench)",
    "PXD070049": "Benchmark (ProteoBench)", "ProteoBench_Module_7": "Benchmark (ProteoBench)",
    "PXD004701": "Bulk cell-line", "PXD030304": "Bulk cell-line", "PXD017199": "Bulk cell-line",
    "PXD041421": "Bulk cell-line", "MSV000093870": "plexDIA", "PXD046357": "Single-cell",
    "PXD071075": "Single-cell", "PXD034128": "Phosphoproteomics", "PXD034623": "Phosphoproteomics",
    "PXD049692": "Phosphoproteomics", "PXD064049": "Spatial DVP",
}
# Label offsets (points) where default placement collides, per panel x column.
OFFSETS = {
    "n_runs": {"PXD062685": (7, 4, "left"), "PXD070049": (7, -6, "left"), "PXD034623": (7, -12, "left"),
               "PXD049412": (7, 6, "left"), "PXD046357": (7, -14, "left"), "Module 7": (7, 0, "left"), "PXD064049": (7, 4, "left")},
    "raw_gb": {"PXD062685": (-8, 8, "right"), "MSV000093870": (7, 10, "left"), "PXD034623": (7, -12, "left"),
               "PXD004701": (7, 8, "left"), "PXD041421": (7, -12, "left"), "PXD049412": (-8, 8, "right")},
}
MARKER = {"Benchmark (ProteoBench)": "o", "Bulk cell-line": "s", "plexDIA": "P", "Single-cell": "D",
          "Phosphoproteomics": "^", "Spatial DVP": "v"}


def head_size(url: str) -> int:
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=120) as r:
        return int(r.headers["Content-Length"])


def fill_sizes(raw: pd.DataFrame) -> pd.DataFrame:
    raw = raw.copy()
    raw["size_bytes"] = pd.to_numeric(raw["size_bytes"], errors="coerce")
    raw["size_note"] = ""
    # .d folders inside zips: sum member sizes per folder (older raw_files.tsv
    # builds recorded the folder's own 0-byte directory entry instead).
    member = raw["archive_member"].fillna("")
    is_dir = member.str.endswith("/") | raw["file_name"].str.lower().str.endswith(".d")
    dirs = raw[is_dir & (raw["size_bytes"].isna() | (raw["size_bytes"] == 0)) & raw["url"].str.endswith(".zip")]
    for url, grp in dirs.groupby("url"):
        prefixes = grp["archive_member"].map(lambda m: m.rstrip("/") + "/")
        totals: dict[str, int] = defaultdict(int)
        for name, size in zip_members(url.replace("ftp://", "https://", 1)):
            for prefix in prefixes:
                if name.startswith(prefix):
                    totals[prefix] += size
        raw.loc[grp.index, "size_bytes"] = prefixes.map(totals)
    # One .7z holding all runs (PXD062685): split evenly, flag as approximate.
    for url, grp in raw[raw["size_bytes"].isna() & raw["url"].str.endswith(".7z")].groupby("url"):
        raw.loc[grp.index, "size_bytes"] = head_size(url) / len(grp)
        raw.loc[grp.index, "size_note"] = "approx: .7z archive size / runs"
    # Plain HTTP files without a recorded size (ProteoBench server).
    for i in raw.index[raw["size_bytes"].isna() & raw["archive_member"].fillna("").eq("")]:
        raw.loc[i, "size_bytes"] = head_size(raw.loc[i, "url"])
    return raw


def main() -> None:
    rb.apply_house_style()
    (OUT / "data").mkdir(parents=True, exist_ok=True)
    perf = pd.read_csv(PERF, sep="\t")
    perf = perf[perf["dataset"] != "PXD003539"].copy()  # same exclusion as Fig. 2b
    # PXD071075 appears once per executor.queueSize sweep point; keep them, tagged.
    perf["queue_size"] = perf["version"].str.extract(r"_q(\d+)$")[0].astype(float)
    raw = fill_sizes(pd.read_csv(RAW, sep="\t", dtype=str))

    rows = []
    for _, p in perf.iterrows():
        ds = p["dataset"]
        sel = raw[raw["used_in"].str.contains(USED_IN[ds], regex=False)]
        rows.append({
            "dataset": ds, "instrument": p["instrument"], "category": CATEGORY[ds],
            "n_runs": int(p["n_runs"]), "raw_gb": sel["size_bytes"].sum() / 1e9,
            "n_raw_files_sized": int(sel["size_bytes"].notna().sum()),
            "size_note": "; ".join(sorted(set(sel["size_note"]) - {""})),
            "wallclock_h": p["wallclock_seconds"] / 3600,
            "queue_size": p["queue_size"],
        })
    df = pd.DataFrame(rows)
    # PXD070049 Module 10 files are ZenoTOF 8600 (LFQ_ZenoTOF8600_*); Fig. 2b's table says 7600.
    df.loc[df["dataset"] == "PXD070049", "instrument"] = "ZenoTOF 8600"
    df["minutes_per_run"] = df["wallclock_h"] * 60 / df["n_runs"]
    df["gb_per_hour"] = df["raw_gb"] / df["wallclock_h"]
    df.sort_values("n_runs").to_csv(OUT / "data/runtime_vs_scale.tsv", sep="\t", index=False, float_format="%.4g")

    fig, axes = plt.subplots(1, 2, figsize=(17, 7.2), gridspec_kw={"wspace": 0.28})
    for ax, xcol, xlabel in ((axes[0], "n_runs", "Raw files processed"),
                             (axes[1], "raw_gb", "Raw data processed (GB)")):
        sweep = df[df["queue_size"].notna()].sort_values("queue_size")
        if len(sweep):
            x0 = sweep[xcol].iloc[0]
            ax.plot([x0, x0], [sweep["wallclock_h"].min(), sweep["wallclock_h"].max()],
                    color=rb.INSTRUMENT_COLORS["Orbitrap Eclipse"], lw=1.5, zorder=2)
            for _, r in sweep[sweep["queue_size"] != 200].iterrows():  # 200 and 300 overlap
                ax.annotate(f"{int(r['queue_size'])}", (r[xcol], r["wallclock_h"]), xytext=(-8, 0),
                            textcoords="offset points", fontsize=10, ha="right", va="center", color="#555555")
            ax.annotate("PXD071075\n(queue size sweep)", (x0, sweep["wallclock_h"].max()),
                        xytext=(7, 0), textcoords="offset points", fontsize=11, color="#333333", va="center")
        for _, r in df[df["queue_size"].isna()].iterrows():
            ax.scatter(r[xcol], r["wallclock_h"], s=150, marker=MARKER[r["category"]],
                       color=rb.INSTRUMENT_COLORS.get(r["instrument"], rb.OKABE_ITO["grey"]),
                       edgecolors="white", linewidths=1.5, zorder=3)
            lab = LABEL.get(r["dataset"], r["dataset"])
            dx, dy, ha = OFFSETS[xcol].get(lab, OFFSETS[xcol].get(r["dataset"], (7, 5, "left")))
            ax.annotate(lab, (r[xcol], r["wallclock_h"]), xytext=(dx, dy), textcoords="offset points",
                        fontsize=11, color="#333333", ha=ha)
        for _, r in sweep.iterrows():
            ax.scatter(r[xcol], r["wallclock_h"], s=90, marker=MARKER[r["category"]],
                       color=rb.INSTRUMENT_COLORS["Orbitrap Eclipse"], edgecolors="white", linewidths=1.2, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Wall-clock (h, excluding in-silico library)")
        ax.grid(True, which="major", color="#E5E5E5", lw=0.8, zorder=0)
        # Linear-scaling reference through the median throughput.
        base = df[df["queue_size"].isna() | (df["queue_size"] == df["queue_size"].max())]
        k = np.median(base["wallclock_h"] / base[xcol])
        xs = np.logspace(np.log10(df[xcol].min() / 1.5), np.log10(df[xcol].max() * 1.5), 50)
        ax.plot(xs, k * xs, color=rb.OKABE_ITO["grey"], ls=":", lw=1.4, zorder=1)
        ax.set_ylim(0.25, 60)
    axes[0].set_title("a  Wall-clock vs number of raw files", loc="left")
    axes[1].set_title("b  Wall-clock vs raw-data volume", loc="left")
    insts = df.drop_duplicates("instrument")["instrument"]
    inst_h = [Line2D([0], [0], marker="o", ls="none", markersize=10,
                     markerfacecolor=rb.INSTRUMENT_COLORS.get(i, "#9e9e9e"), markeredgecolor="white", label=i)
              for i in insts]
    cat_h = [Line2D([0], [0], marker=MARKER[c], ls="none", markersize=10, markerfacecolor="#777777",
                    markeredgecolor="white", label=c) for c in MARKER if c in set(df["category"])]
    ref_h = [Line2D([0], [0], color=rb.OKABE_ITO["grey"], ls=":", label="linear scaling at median throughput")]
    fig.legend(handles=inst_h, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.13), title="Instrument")
    fig.legend(handles=cat_h + ref_h, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25), title="Sample type")
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"supp_runtime_vs_scale.{ext}")
    print(df.sort_values("n_runs").to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
