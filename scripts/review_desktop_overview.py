"""Supplementary figure (R1 Q1.1, Q1.7): standalone DIA-NN vs quantmsdiann overview.

Publication-quality redraw of the overlap / intensity / cumulative-ID figure
made for the review experiments. One row per ProteoBench dataset:
  a  precursor overlap (area-proportional Venn)
  b  protein-group overlap
  c  log2 precursor quantity, quantmsdiann vs standalone (shared run x precursor)
  d  cumulative precursors vs global q-value threshold
  e  cumulative protein groups vs global protein-group q-value threshold

Global rule as in methods.md section 1: Lib.Q.Value / Lib.PG.Q.Value <= 0.01,
decoys dropped. Panel c uses run-level Q.Value <= 0.01 in both workflows.

Usage:
    python -m scripts.review_desktop_overview
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch

from scripts import rebuild as rb
from scripts.review_desktop_vs_quantmsdiann import DATASETS, OUT, WORKFLOWS, load, quant_pairs

Q = 0.01
COLOR = {"desktop": rb.OKABE_ITO["grey"], "quantmsdiann": rb.OKABE_ITO["blue"]}
LABEL = {"desktop": "DIA-NN 2.5.1 standalone", "quantmsdiann": "quantmsdiann (DIA-NN 2.5.1)"}


def _lens_area(r1: float, r2: float, d: float) -> float:
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return np.pi * min(r1, r2) ** 2
    a1 = r1 ** 2 * np.arccos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    a2 = r2 ** 2 * np.arccos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    return a1 + a2 - 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))


def venn(ax, a: set, b: set, title: str) -> None:
    """Two-set area-proportional Venn (a = standalone, b = quantmsdiann)."""
    n_a, n_b, n_ab = len(a), len(b), len(a & b)
    r1, r2 = np.sqrt(n_a / np.pi), np.sqrt(n_b / np.pi)
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):  # bisection on centre distance to match the shared area
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if _lens_area(r1, r2, mid) > n_ab else (lo, mid)
    d = (lo + hi) / 2
    scale = 1 / max(r1, r2)
    r1, r2, d = r1 * scale, r2 * scale, d * scale
    x1, x2 = -d / 2, d / 2
    for (x, r, wf) in ((x1, r1, "desktop"), (x2, r2, "quantmsdiann")):
        ax.add_patch(Circle((x, 0), r, facecolor=COLOR[wf], alpha=0.45, edgecolor=COLOR[wf], lw=1.8))
    ax.text((x1 + x2) / 2, 0.08, f"{n_ab:,}\nshared", ha="center", va="center", fontsize=17, fontweight="bold")
    union = n_a + n_b - n_ab
    ax.text(x1 - r1, 1.12, f"{n_a - n_ab:,} standalone only", ha="left", va="bottom", fontsize=13,
            color="#444444")
    ax.text(x2 + r2, -1.12, f"{n_b - n_ab:,} quantmsdiann only", ha="right", va="top", fontsize=13,
            color=rb.OKABE_ITO["blue"])
    ax.text(0, -0.42, f"Jaccard {n_ab / union * 100:.1f}%", ha="center", va="center", fontsize=12, color="#222222")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.45, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left")


def cumulative(ax, dfs: dict, qcol: str, idcol: str, ylabel: str, title: str) -> None:
    for wf in WORKFLOWS:
        best = dfs[wf].groupby(idcol)[qcol].min().sort_values()
        best = best[best <= 0.05]
        ax.step(best.values, np.arange(1, len(best) + 1), where="post", color=COLOR[wf], lw=2.2)
        n1 = int((best <= Q).sum())
        ax.text(0.97, 0.08 if wf == "desktop" else 0.18, f"{n1:,} at 1% FDR", transform=ax.transAxes,
                ha="right", fontsize=13, color="#444444" if wf == "desktop" else COLOR[wf])
    ax.axvline(Q, color="#777777", ls="--", lw=1)
    ax.set_xlim(0, 0.05)
    ax.set_ylim(0, None)
    ax.set_xticks([0, 0.01, 0.02, 0.03, 0.04, 0.05], ["0", "1", "2", "3", "4", "5"])
    ax.set_xlabel("Global q-value threshold (%)")
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(rb.kfmt))
    ax.set_title(title, loc="left")


def main() -> None:
    rb.apply_house_style()
    n = len(DATASETS)
    fig = plt.figure(figsize=(28, 6.4 * n))
    gs = fig.add_gridspec(n, 5, width_ratios=[1, 1, 1, 1, 1], wspace=0.34, hspace=0.42)
    for r, (ds, meta) in enumerate(DATASETS.items()):
        dfs = {wf: load(ds, wf) for wf in WORKFLOWS}
        prec = {wf: set(d.loc[d["Lib.Q.Value"] <= Q, "Precursor.Id"]) for wf, d in dfs.items()}
        pg = {wf: set(d.loc[d["Lib.PG.Q.Value"] <= Q, "Protein.Group"]) for wf, d in dfs.items()}
        L = "abcde" if r == 0 else "fghij"
        venn(fig.add_subplot(gs[r, 0]), prec["desktop"], prec["quantmsdiann"], f"{L[0]}  Precursors (1% FDR)")
        venn(fig.add_subplot(gs[r, 1]), pg["desktop"], pg["quantmsdiann"], f"{L[1]}  Protein groups (1% FDR)")

        ax = fig.add_subplot(gs[r, 2])
        m = quant_pairs(dfs["desktop"], dfs["quantmsdiann"])
        lims = (np.floor(m.min().min()), np.ceil(m.max().max()))
        ax.hexbin(m["desktop"], m["quantmsdiann"], gridsize=90, bins="log", cmap="Blues", mincnt=1,
                  extent=(*lims, *lims), linewidths=0)
        ax.plot(lims, lims, color=rb.OKABE_ITO["vermillion"], ls="--", lw=1.4)
        r_ = np.corrcoef(m["desktop"], m["quantmsdiann"])[0, 1]
        ax.text(0.04, 0.96, f"Pearson r = {r_:.3f}\nn = {len(m):,}", transform=ax.transAxes, va="top", fontsize=12)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal")
        ax.set_xlabel("Standalone log$_2$ quantity")
        ax.set_ylabel("quantmsdiann log$_2$ quantity")
        ax.set_title(f"{L[2]}  Precursor quantities", loc="left")

        cumulative(fig.add_subplot(gs[r, 3]), dfs, "Lib.Q.Value", "Precursor.Id", "Precursors",
                   f"{L[3]}  Cumulative precursors")
        cumulative(fig.add_subplot(gs[r, 4]), dfs, "Lib.PG.Q.Value", "Protein.Group", "Protein groups",
                   f"{L[4]}  Cumulative protein groups")
        fig.axes[-5].annotate(meta["label"], xy=(0, 1.22), xycoords="axes fraction", fontsize=20,
                              fontweight="bold", va="bottom")
    handles = [Patch(facecolor=COLOR[w], alpha=0.7, label=LABEL[w]) for w in WORKFLOWS] + [
        Line2D([0], [0], color=rb.OKABE_ITO["vermillion"], ls="--", label="y = x"),
        Line2D([0], [0], color="#777777", ls="--", label="1% FDR")]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.05), fontsize=16)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(OUT.parent / "supplementary_revision" / f"FigS6_standalone_vs_quantmsdiann.{ext}")


if __name__ == "__main__":
    main()
