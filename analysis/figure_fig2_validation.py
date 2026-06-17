#!/usr/bin/env python
"""Fig 2 - quantmsdiann validation row (scaling + ProteoBench accuracy).

One row of three panels (the former Fig 1 b/c/d, now a standalone figure so the
workflow gets Fig 1 to itself):
  (a) wall-clock versus cluster nodes (PXD071075 single-cell sweep)
  (b) wall-clock to finish each reanalysis, one bar per dataset
  (c) ProteoBench quantification-accuracy concordance vs standalone DIA-NN

Reuses the existing per-panel renderers (composite/ax mode) so the numbers stay
identical to the standalone figures.

Out: analysis/figures/manuscript/fig2_validation.svg
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analysis import figure_style as fs
fs.apply_house_style()
from analysis.figure_queue_size_sweep import render_queue_size_sweep
from analysis.figure_performance_trace import render_parallelism_scatter
from analysis import figure_proteobench_accuracy as acc

REPO = Path(__file__).resolve().parents[1]
PERF = REPO / "analysis" / "figures" / "performance" / "data"
OUT = REPO / "analysis" / "figures" / "manuscript" / "fig2_validation.svg"


def render(out: Path) -> Path:
    dq = pd.read_csv(PERF / "queue_size_sweep.tsv", sep="\t")
    dp = pd.read_csv(PERF / "parallelism_data.tsv", sep="\t")
    # Narrower + taller than a wide strip: at \textwidth the figure is downscaled
    # less (~0.75x vs ~0.43x), so panels and fonts render markedly larger.
    fig, ax = plt.subplots(1, 3, figsize=(9.0, 5.4),
                           gridspec_kw={"width_ratios": [1.0, 1.3, 1.0]})
    render_queue_size_sweep(dq, ax=ax[0], composite=True)
    render_parallelism_scatter(dp, ax=ax[1], composite=True, show_legend=True, legend_ncol=5)
    acc.draw(ax[2], compact=True)
    for a, lab in zip(ax, "abc"):
        a.text(-0.10, 1.05, f"({lab})", transform=a.transAxes, fontsize=15,
               fontweight="bold", va="bottom", ha="right")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    print(f"wrote {render(OUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
