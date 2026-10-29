"""Final supplementary figures and tables for the revision (Figs. S7-S8, Tables S3-S5).

Reads the statistics written by the three review analysis scripts and makes
one simple figure plus one table per reviewer point:

  S-desktop  (R1 Q1.1, Q1.7, Q1.18)  review_desktop_vs_quantmsdiann
  S-runtime  (R1 Q1.24)              review_runtime_vs_scale
  S-evidence (R1 Q1.10, Q1.16, Q1.21) review_threshold_regimes

Usage (after the three analysis scripts have run):
    python -m scripts.review_supplementary
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts import rebuild as rb

REPO = Path(__file__).resolve().parents[1]
REVIEW = REPO / "analysis/figures/review"
OUT = REVIEW / "supplementary_revision"  # final supplementary items, final names
GREY, BLUE = rb.COMPARISON["original"], rb.COMPARISON["quantmsdiann"]
DS_LABEL = {"PXD049412": "Single-cell Astral\n(PXD049412)", "PXD070049": "ZenoTOF 8600\n(PXD070049)"}


def save(fig, name: str) -> None:
    for ext in ("svg", "pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)


def bars(ax, groups, left, right, labels, ylabel, fmt="{:,.0f}"):
    x = np.arange(len(groups))
    w = 0.38
    for off, vals, col, lab in ((-w / 2, left, GREY, labels[0]), (w / 2, right, BLUE, labels[1])):
        ax.bar(x + off, vals, w * 0.92, color=col, label=lab)
        for xi, v in zip(x + off, vals):
            ax.text(xi, v, fmt.format(v), ha="center", va="bottom", fontsize=11)
    ax.set_xticks(x, groups)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, max(max(left), max(right)) * 1.15)


# ---------------------------------------------------------------- Q1.1 / Q1.7

def _shared(v, ds, g) -> str:
    """Shared count, Jaccard and share of each workflow's IDs found by the other, from raw counts."""
    sh = int(v(ds, f"global {g} shared", "both"))
    d, q = int(v(ds, f"global {g} only", "desktop")), int(v(ds, f"global {g} only", "quantmsdiann"))
    return (f"{sh:,} shared; Jaccard {sh / (sh + d + q) * 100:.1f}%; "
            f"{sh / (sh + d) * 100:.1f}% of standalone / {sh / (sh + q) * 100:.1f}% of quantmsdiann found by the other")


def desktop() -> pd.DataFrame:
    s = pd.read_csv(REVIEW / "desktop_vs_quantmsdiann/data/desktop_vs_quantmsdiann_stats.tsv", sep="\t")
    v = lambda ds, metric, wf: s.loc[(s.dataset == ds) & (s.metric == metric) & (s.workflow == wf), "value"].iloc[0]
    rows = []
    for ds in DS_LABEL:
        for lvl, key in (("Precursors", "precursors"), ("Peptides", "peptides"), ("Protein groups", "protein_groups")):
            d, q = v(ds, f"per-run mean {key}", "desktop"), v(ds, f"per-run mean {key}", "quantmsdiann")
            g = {"Precursors": "Precursors", "Peptides": "Peptide sequences", "Protein groups": "Protein groups"}[lvl]
            rows.append({"Dataset": ds, "Metric": f"{lvl} per run (mean)", "DIA-NN standalone": round(d),
                         "quantmsdiann": round(q), "Difference": f"{(q - d) / d * 100:+.1f}%"})
            rows.append({"Dataset": ds, "Metric": f"{lvl}, global 1% FDR: overlap",
                         "DIA-NN standalone": int(v(ds, f"global {g} only", "desktop")),
                         "quantmsdiann": int(v(ds, f"global {g} only", "quantmsdiann")),
                         "Difference": _shared(v, ds, g)})
        rows.append({"Dataset": ds, "Metric": "Precursor quantity, log2 ratio quantmsdiann/standalone (median)",
                     "DIA-NN standalone": "", "quantmsdiann": "",
                     "Difference": f"{v(ds, 'quant median log2 diff (qmd - desktop)', 'both'):+.2f} "
                                   f"(r = {v(ds, 'quant Pearson r (log2)', 'both'):.3f})"})
        for c in "AB":
            rows.append({"Dataset": ds, "Metric": f"Median precursor CV, condition {c} (%)",
                         "DIA-NN standalone": v(ds, f"median precursor CV cond {c} (shared, complete)", "desktop"),
                         "quantmsdiann": v(ds, f"median precursor CV cond {c} (shared, complete)", "quantmsdiann"),
                         "Difference": ""})
        for sp in ("HUMAN", "YEAST", "ECOLI"):
            m = f"{sp} median |log2 A/B - expected|"
            if (s.metric == m).any() and ((s.dataset == ds) & (s.metric == m)).any():
                rows.append({"Dataset": ds, "Metric": f"{sp.capitalize()}: median |log2 A/B - expected|",
                             "DIA-NN standalone": v(ds, m, "desktop"), "quantmsdiann": v(ds, m, "quantmsdiann"),
                             "Difference": ""})
        rows.append({"Dataset": ds, "Metric": "Report rows (run x precursor)",
                     "DIA-NN standalone": int(v(ds, "report rows", "desktop")),
                     "quantmsdiann": int(v(ds, "report rows", "quantmsdiann")),
                     "Difference": f"{int(v(ds, 'report rows in one workflow only', 'both'))} rows in one workflow only"})
        rows.append({"Dataset": ds, "Metric": "Protein groups, Global.PG.Q.Value <= 1%",
                     "DIA-NN standalone": int(v(ds, "protein groups Global.PG.Q.Value <= 0.01", "desktop")),
                     "quantmsdiann": int(v(ds, "protein groups Global.PG.Q.Value <= 0.01", "quantmsdiann")), "Difference": ""})
        only_d = int(v(ds, "pg_matrix protein groups in this workflow only", "desktop"))
        only_q = int(v(ds, "pg_matrix protein groups in this workflow only", "quantmsdiann"))
        n_d = int(v(ds, "pg_matrix protein groups", "desktop"))
        rows.append({"Dataset": ds, "Metric": "pg_matrix protein groups",
                     "DIA-NN standalone": n_d, "quantmsdiann": int(v(ds, "pg_matrix protein groups", "quantmsdiann")),
                     "Difference": f"{only_d} standalone only, {only_q} quantmsdiann only "
                                   f"({(only_d + only_q) / (n_d + only_q) * 100:.1f}%, all near the 1% threshold); "
                                   f"{int(v(ds, 'pg_matrix shared quantities differing', 'both'))} differing quantities in shared groups"})
        differing = s[(s.dataset == ds) & s.metric.str.startswith("rows differing: ")]
        rows.append({"Dataset": ds, "Metric": "Report columns with any differing value",
                     "DIA-NN standalone": "", "quantmsdiann": "",
                     "Difference": "; ".join(
                         f"{r.metric.removeprefix('rows differing: ')} ({int(r.value):,} rows, max "
                         f"{v(ds, 'max abs difference: ' + r.metric.removeprefix('rows differing: '), 'both'):.4f})"
                         for r in differing.itertuples()) or "none"})
    tab = pd.DataFrame(rows)
    tab.loc[tab.Metric.str.contains("overlap"), ["DIA-NN standalone", "quantmsdiann"]] = \
        tab.loc[tab.Metric.str.contains("overlap"), ["DIA-NN standalone", "quantmsdiann"]].map(lambda x: f"{x:,} only")

    # Figure S6 for this comparison is drawn by scripts.review_desktop_overview.
    return tab


# ---------------------------------------------------------------- Q1.24

def runtime() -> pd.DataFrame:
    df = pd.read_csv(REVIEW / "runtime_vs_scale/data/runtime_vs_scale.tsv", sep="\t")
    df = df[df.queue_size.isna() | (df.queue_size == 300)].sort_values("n_runs")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df.n_runs, df.wallclock_h, s=70, color=BLUE, zorder=3)
    for _, r in df.iterrows():
        dy = -14 if r.dataset == "PXD034623" else 4
        ax.annotate(r.dataset.replace("ProteoBench_Module_7", "Module 7"), (r.n_runs, r.wallclock_h),
                    xytext=(5, dy), textcoords="offset points", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Raw files processed")
    ax.set_ylabel("Wall-clock (h)")
    ax.grid(True, color="#E5E5E5", lw=0.8, zorder=0)
    save(fig, "FigS7_runtime_vs_files")
    return pd.DataFrame({
        "Dataset": df.dataset, "Instrument": df.instrument, "Sample type": df.category,
        "Raw files": df.n_runs, "Raw data (GB)": df.raw_gb.round(1), "Wall-clock (h)": df.wallclock_h.round(2),
        "Min per file": df.minutes_per_run.round(2), "GB per hour": df.gb_per_hour.round(1),
    })


# ---------------------------------------------------------------- Q1.10 / Q1.21

def evidence() -> pd.DataFrame:
    t = pd.read_csv(REVIEW / "threshold_regimes/data/threshold_regimes.tsv", sep="\t")
    overlap = pd.read_csv(REVIEW / "threshold_regimes/data/peptide_overlap.tsv", sep="\t")
    cohorts = {"PXD003539": "NCI-60", "PXD004701": "Sun 2023", "PXD030304": "ProCan"}
    qrule = "DIA-NN report: Lib.PG.Q <= 1% groups, distinct Lib.Q <= 1% peptides"
    rows = []
    for ds, name in cohorts.items():
        sub = t[t.dataset == ds]
        orule = sub.loc[(sub.side == "original") & sub.rule.str.contains("distinct peptide sequences"), "rule"].iloc[0]
        get = lambda side, rule, metric: sub.loc[(sub.side == side) & (sub.rule == rule) & (sub.metric == metric), "value"]
        row = {"Cohort": f"{name} ({ds})"}
        for k in (1, 2, 3, 5):
            o, q = get("original", orule, f"proteins >= {k} peptides").iloc[0], get("quantmsdiann", qrule, f"proteins >= {k} peptides").iloc[0]
            row[f">={k} pep: original"], row[f">={k} pep: quantmsdiann"] = int(o), int(q)
        for m in ("precursors", "peptides"):
            o = sub.loc[(sub.side == "original") & (sub.metric == m), "value"].iloc[0]
            row[f"{m.capitalize()}: original"] = "n/a" if pd.isna(o) else int(o)
            row[f"{m.capitalize()}: quantmsdiann"] = int(sub.loc[(sub.side == "quantmsdiann") & (sub.metric == m), "value"].iloc[0])
        for side in ("original", "quantmsdiann"):
            row[f"Median peptides/protein: {side}"] = int(sub.loc[(sub.side == side) & (sub.metric == "median peptides per protein"), "value"].iloc[0])
        ov = overlap[overlap.dataset == ds]
        row["Peptides shared with original"] = "n/a" if ov.empty else f"{int(ov.shared.iloc[0]):,}"
        row["% of original peptides recovered"] = "n/a" if ov.empty else float(ov.pct_original_recovered.iloc[0])
        rows.append(row)
    tab = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw={"wspace": 0.3})
    for ax, (_, r), letter in zip(axes, tab.iterrows(), "abc"):
        ks = (1, 2, 3, 5)
        bars(ax, [f"≥{k}" for k in ks], [r[f">={k} pep: original"] for k in ks],
             [r[f">={k} pep: quantmsdiann"] for k in ks], ("Original", "quantmsdiann"), "Proteins")
        ax.set_xlabel("Peptides per protein")
        ax.set_title(f"{letter}  {r.Cohort}", loc="left")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(rb.kfmt))
    axes[0].legend(loc="upper right")
    save(fig, "FigS8_proteins_by_peptide_evidence")
    return tab


def main() -> None:
    rb.apply_house_style()
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {"TableS3_standalone_vs_quantmsdiann": desktop(), "TableS4_runtime_vs_scale": runtime(),
              "TableS5_proteins_by_peptide_evidence": evidence()}
    for name, tab in tables.items():
        tab.to_csv(OUT / f"{name}.tsv", sep="\t", index=False)
        print(f"\n== {name}\n{tab.to_string(index=False)}")


if __name__ == "__main__":
    main()
