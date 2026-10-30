"""Reviewer response (R1 Q1.9), response letter only: figure of the cell-line
cohort integration with the quantms criteria (not part of the paper or the
supplementary material).

Reads integration_fdr.tsv, written by scripts.review_integration_fdr (harmonised
DIA-NN 2.5.1 runs of the five cell-line cohorts). The quantms 2024 cell-line numbers (Dai et al.
2024, Nat Methods) are given for comparison.

Usage: python -m scripts.review_integration_figure
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from scripts.review_integration_fdr import OUT as DATA_TSV

from scripts.review_integration_fdr import TRACKED as OUT
QUANTMS_2024 = (11403, 11374)  # cell lines: merged, after 1% adjusted FDR
COHORT, MERGED, KEPT, INK, MUTED = "#9aa5b1", "#4c6a92", "#1f3a5f", "#222222", "#666666"


def main() -> int:
    df = pd.read_csv(DATA_TSV, sep="\t")
    df = df[df["rule"] == ">=2_unique_peptides"].set_index("set")
    cohorts = df.iloc[:5]
    before, after = df.loc["merged"], df.loc["decoy_1pct"]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})

    labels = list(cohorts.index) + ["Merged", "Adjusted q < 0.01"]
    values = list(cohorts.proteins) + [before.proteins, after.proteins]
    colors = [COHORT] * 5 + [MERGED, KEPT]
    y = range(len(labels))[::-1]
    a.barh(list(y), values, color=colors, height=0.7, edgecolor="white", linewidth=2)
    for yi, v in zip(y, values):
        note = ""
        if v == after.proteins:
            note = f"  ({before.proteins - after.proteins} removed)"
        a.text(v + 150, yi, f"{v:,}", va="center", fontsize=8, color=INK)
        a.text(v + 150 + 330 * len(f"{v:,}"), yi, note, va="center", fontsize=7, color=MUTED)
    a.set_yticks(list(y), labels, fontsize=8)
    a.set_xlim(0, 18500)
    a.set_xlabel("Proteins (≥2 unique peptides, 1% protein FDR)", fontsize=8)
    a.set_title("a  Cell-line cohorts integrated with the quantms criteria", loc="left", fontsize=9, color=INK)

    groups = ["quantms 2024", "quantmsdiann"]
    pairs = [QUANTMS_2024, (before.proteins, after.proteins)]
    w = 0.36
    for i, (g, (n0, n1)) in enumerate(zip(groups, pairs)):
        b.bar(i - w / 2, n0, w, color=MERGED, edgecolor="white", linewidth=2, label="Merged" if i == 0 else None)
        b.bar(i + w / 2, n1, w, color=KEPT, edgecolor="white", linewidth=2, label="Adjusted q < 0.01" if i == 0 else None)
        b.text(i, max(n0, n1) + 250, f"−{n0 - n1} ({(n0 - n1) / n0:.1%})", ha="center", fontsize=8, color=INK)
    b.set_xticks(range(2), groups, fontsize=8)
    b.set_ylim(0, 14500)
    b.set_ylabel("Proteins", fontsize=8)
    b.legend(frameon=False, fontsize=8, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.32))
    b.set_title("b  Removed by re-estimation", loc="left", fontsize=9, color=INK)

    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=8)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(OUT / f"response_R1.9_integration.{ext}", dpi=300)
    print(OUT / "response_R1.9_integration.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
