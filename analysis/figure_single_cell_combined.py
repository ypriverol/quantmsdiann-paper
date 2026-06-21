#!/usr/bin/env python
"""Single-cell reanalysis figure - DIA-NN 1.8.1 vs 2.5.1 Enterprise, two
label-free single-cell cohorts: HeLa Astral (PXD046357) and A549/H460
(PXD049412; the 20x/40x A549 library runs are excluded from per-cell counts).
2x2 layout:
  A  (top, spanning) total precursors + total protein groups + per-cell protein
     groups (box + jitter), both cohorts x build, with per-build %change.
  B  Data-completeness curve (HeLa Astral flagship; >= N cells), y from 0.
  C  CV across cells -- quantitative precision (both cohorts; Astral solid,
     A549/H460 dashed).

Data provenance
---------------
All numbers are derived from the deposited DIA-NN reports by
``analysis/make_single_cell_tables.py`` (run it to (re)generate the inputs):
  * mv_{per_cell,completeness,rank_abundance,cv}.tsv and sc_totals.tsv
    <- PRIDE FTP quantmsdiann-benchmarks/single-cell/{PXD046357,PXD049412}/
       v{1_8_1,2_5_1_enterprise}/quant_tables/diann_report.{tsv,parquet}
       (our reanalysis; counting via analysis/count_report_ids.py under the
       methods.md filter rule: per-cell PG.Q.Value, totals Lib.*).
Pipeline: https://github.com/bigbio/quantmsdiann

Data: data/single_cell/mv_*.tsv, data/single_cell/sc_totals.tsv.

Run:  python -m analysis.figure_single_cell_combined
Out:  analysis/figures/manuscript/fig3_single_cell_combined.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()

REPO = Path(__file__).resolve().parents[1]
D = REPO / "data" / "single_cell"
OUT = REPO / "analysis" / "figures" / "manuscript" / "fig3_single_cell_combined.svg"

VERS = ["1_8_1", "2_5_1_enterprise"]
VLAB = {"1_8_1": "1.8.1", "2_5_1_enterprise": "2.5.1 Enterprise"}
VCOL = {v: fs.VERSION_COLORS[v] for v in VERS}
ACC = {"HeLa Astral SC": "PXD046357", "A549/H460 SC": "PXD049412"}
# accession by the short panel label used in the merged panel
ACC_SHORT = {"Astral": "PXD046357", "A549/H460": "PXD049412"}
DS_STYLE = {"HeLa Astral SC": "-", "A549/H460 SC": "--"}
FLAG = "HeLa Astral SC"  # cohort carrying the data-completeness panel


def _completeness(ax):
    # Flagship cohort only (HeLa Astral): the two cohorts have very different
    # cell counts (12 vs ~146), so a shared "in >= N cells" axis would distort.
    df = pd.read_csv(D / "mv_completeness.tsv", sep="\t")
    df = df[df["dataset"] == FLAG]
    for v in VERS:
        s = df[df["version"] == v].sort_values("min_cells")
        ax.plot(s["min_cells"], s["n_proteins"], linestyle="-",
                marker="o", ms=3.5, lw=1.8, color=VCOL[v])
    ax.set_xlabel("quantified in ≥ N cells"); ax.set_ylabel("protein groups")
    ax.set_ylim(bottom=0)
    ax.set_title("Data completeness")
    fs.kfmt_axis(ax.yaxis); fs.despine(ax)


def _cv(ax):
    df = pd.read_csv(D / "mv_cv.tsv", sep="\t")
    # clip at 1.5 (not 1.0): low-input single-cell data has a real high-CV tail;
    # clipping at 1.0 would pile that tail into one bin, an artefactual spike.
    bins = np.linspace(0, 1.5, 46)
    for ds in df["dataset"].unique():
        for v in VERS:
            cv = df[(df["dataset"] == ds) & (df["version"] == v)]["cv"].clip(0, 1.5)
            ax.hist(cv, bins=bins, density=True, histtype="step", linewidth=1.5,
                    color=VCOL[v], linestyle=DS_STYLE.get(ds, "-"))
    ax.set_xlabel("CV across cells"); ax.set_ylabel("density")
    ax.set_xlim(0, 1.5)
    ax.set_title("Quantitative precision")
    fs.despine(ax)


# Total counts (report-based, target-only): precursors and protein groups,
# 1.8.1 -> 2.5.1 Enterprise, read from the generated sc_totals.tsv (no hardcoding;
# see analysis/make_single_cell_tables.py).
_SHORT = {"HeLa Astral SC": "Astral", "A549/H460 SC": "A549/H460"}


def _load_totals() -> dict:
    df = pd.read_csv(D / "sc_totals.tsv", sep="\t")
    out: dict = {}
    for ds, g in df.groupby("dataset"):
        gv = g.set_index("version")
        out[_SHORT.get(ds, ds)] = {
            m: (int(gv.loc["1_8_1", m]), int(gv.loc["2_5_1_enterprise", m]))
            for m in ("precursors", "proteins")
        }
    return {k: out[k] for k in ("Astral", "A549/H460") if k in out}


_FULL = {"Astral": "HeLa Astral SC", "A549/H460": "A549/H460 SC"}


def _merged(ax):
    """Merged A+B panel: three x-sections sharing the figure, 1.8.1 vs 2.5.1
    Enterprise. (i) total precursors (left axis), (ii) total protein groups and
    (iii) per-cell protein groups (box+jitter) both on the right axis (same
    scale). Replaces the separate totals + per-cell panels."""
    bw = 0.36
    ax2 = ax.twinx()
    TOTALS = _load_totals()
    dsx = list(TOTALS)                                   # [Astral]
    percell = pd.read_csv(D / "mv_per_cell.tsv", sep="\t")
    prec_x = {d: i for i, d in enumerate(dsx)}                       # 0, 1
    prot_x = {d: i + len(dsx) + 0.6 for i, d in enumerate(dsx)}      # 2.6, 3.6
    cell_x = {d: i + 2 * len(dsx) + 1.2 for i, d in enumerate(dsx)}  # 5.2, 6.2
    rng = np.random.default_rng(0)
    for k, v in enumerate(VERS):
        idx = 0 if v == "1_8_1" else 1
        for d in dsx:
            # (i) precursors total — left axis
            xp = prec_x[d] + (k - 0.5) * bw
            hp = TOTALS[d]["precursors"][idx]
            ax.bar(xp, hp, bw, color=VCOL[v], edgecolor="white", linewidth=0.6)
            if v != "1_8_1":
                lo = TOTALS[d]["precursors"][0]
                ax.annotate(f"{round(100*(hp-lo)/lo):+d}%", (xp, hp), textcoords="offset points",
                            xytext=(0, 3), ha="center", va="bottom", fontsize=9, fontweight="bold", color=VCOL[v])
            # (ii) protein groups total — right axis
            xt = prot_x[d] + (k - 0.5) * bw
            hg = TOTALS[d]["proteins"][idx]
            ax2.bar(xt, hg, bw, color=VCOL[v], edgecolor="white", linewidth=0.6)
            if v != "1_8_1":
                lo = TOTALS[d]["proteins"][0]
                ax2.annotate(f"{round(100*(hg-lo)/lo):+d}%", (xt, hg), textcoords="offset points",
                             xytext=(0, 3), ha="center", va="bottom", fontsize=9, fontweight="bold", color=VCOL[v])
            # (iii) per-cell protein groups — right axis (box + jitter)
            xc = cell_x[d] + (k - 0.5) * bw
            vals = percell[(percell["dataset"] == _FULL[d]) & (percell["version"] == v)]["pg_count"].values
            bp = ax2.boxplot([vals], positions=[xc], widths=bw * 0.85, patch_artist=True, showfliers=False)
            fs.style_boxplot(bp, color=VCOL[v])
            ax2.scatter(xc + rng.uniform(-0.07, 0.07, len(vals)), vals, s=9, color=VCOL[v],
                        alpha=0.6, edgecolors="none", zorder=3)
    ax.axvline(len(dsx) - 0.2, color="#cccccc", linewidth=0.8)
    ax.axvline(prot_x[dsx[-1]] + 0.7, color="#cccccc", linewidth=0.8)
    ticks = [prec_x[d] for d in dsx] + [prot_x[d] for d in dsx] + [cell_x[d] for d in dsx]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{ACC_SHORT.get(d, d)}\n({d})" for d in dsx] * 3, fontsize=8.5)
    ax.set_xlim(-0.7, cell_x[dsx[-1]] + 0.7)
    ax.set_ylabel("precursors"); ax2.set_ylabel("protein groups")
    prec_max = max(TOTALS[d]["precursors"][i] for d in dsx for i in (0, 1))
    prot_max = max(TOTALS[d]["proteins"][i] for d in dsx for i in (0, 1))
    cell_max = float(percell[percell["dataset"].isin([_FULL[d] for d in dsx])]["pg_count"].max())
    ax.set_ylim(0, prec_max * 1.15); ax2.set_ylim(0, max(prot_max, cell_max) * 1.12)
    fs.kfmt_axis(ax.yaxis); fs.kfmt_axis(ax2.yaxis)
    ax.set_title("Total identifications and per-cell protein groups")
    for xs, lab in ((prec_x, "precursors\n(total)"), (prot_x, "protein groups\n(total)"),
                    (cell_x, "protein groups\n(per cell)")):
        ax.text(np.mean([xs[d] for d in dsx]), -0.16, lab, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10, fontweight="bold")
    for sp in ("top",):
        ax.spines[sp].set_visible(False); ax2.spines[sp].set_visible(False)


def render(out: Path) -> Path:
    # plexDIA (old F) moved to the reanalysis figure; dynamic range (old D)
    # dropped. Layout: merged totals+per-cell spanning the top, then the two
    # mechanistic panels (completeness, CV) below.
    fig = plt.figure(figsize=(10.5, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.42, wspace=0.28)
    ax_merged = fig.add_subplot(gs[0, :])
    ax_comp = fig.add_subplot(gs[1, 0])
    ax_cv = fig.add_subplot(gs[1, 1])
    _merged(ax_merged); _completeness(ax_comp); _cv(ax_cv)
    handles = [Line2D([0], [0], color=VCOL[v], marker="o", linewidth=2, markersize=8,
               label=f"DIA-NN {VLAB[v]}") for v in VERS]
    handles += [Line2D([0], [0], color="#555555", linestyle="-", linewidth=2, label="PXD046357 (HeLa Astral)"),
                Line2D([0], [0], color="#555555", linestyle="--", linewidth=2, label="PXD049412 (A549/H460)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.01), fontsize=11)
    for a, lab in zip([ax_merged, ax_comp, ax_cv], "ABC"):
        a.text(-0.06, 1.05, lab, transform=a.transAxes, fontsize=17, fontweight="bold", va="bottom", ha="right")
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=2.6, w_pad=2.2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out); plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
