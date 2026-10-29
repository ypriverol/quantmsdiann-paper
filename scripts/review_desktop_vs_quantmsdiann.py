"""Reviewer response (R1 Q1.1, Q1.7, Q1.18): desktop DIA-NN vs quantmsdiann.

Same raw files, same DIA-NN 2.5.0 container, same FASTA, predicted library and
search settings, on the same cluster: one run as a single standalone DIA-NN
command (`standalone_sorted.sh`, files passed with --f in name order, as
quantmsdiann orders them), one through quantmsdiann v2.3.0 (`quantmsdiann_dq.sh`).
Datasets: ProteoBench Module 9 (PXD049412, single-cell Astral) and Module 10
(PXD070049, ZenoTOF 8600).

Extends the overlap/scatter/cumulative-ID plots already made for the response
with what Reviewer 1 asked for and those plots do not show: per-run depth,
overlap at precursor / peptide / protein-group level, where workflow-specific
identifications sit in the other workflow's q-value distribution, quantitative
bias (Bland-Altman), precision (CV) and accuracy (species log2 A/B vs expected).

Counting follows methods.md section 1: per-run Q.Value / PG.Q.Value <= 0.01,
global Lib.Q.Value / Lib.PG.Q.Value <= 0.01, decoys dropped, nothing else.

Inputs (git-ignored, copied from the Codon iteration
2026-10-28_01_run_paired-standalone-vs-quantmsdiann; launch scripts in scripts/):
    data/proteobench-review-experiments/codon/<PXD>/standalone/diann_report.parquet
    data/proteobench-review-experiments/codon/<PXD>/quantmsdiann/results/quant_tables/diann_report.parquet

Usage:
    python -m scripts.review_desktop_vs_quantmsdiann
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from scripts import rebuild as rb

REPO = Path(__file__).resolve().parents[1]
X = REPO / "data/proteobench-review-experiments/codon"
OUT = REPO / "analysis/figures/review/desktop_vs_quantmsdiann"
Q = 0.01

DATASETS = {
    "PXD049412": {"label": "Single-cell Astral (PXD049412, Module 9)",
                  "short": "Single-cell Astral",
                  "cond_A": "240pg_10pg", "module": "quant_lfq_DIA_ion_singlecell"},
    "PXD070049": {"label": "ZenoTOF 8600 (PXD070049, Module 10)",
                  "short": "ZenoTOF 8600",
                  "cond_A": "Condition_A", "module": "quant_lfq_DIA_ion_ZenoTOF"},
}
WORKFLOWS = ["desktop", "quantmsdiann"]
WF_LABEL = {"desktop": "DIA-NN standalone", "quantmsdiann": "quantmsdiann"}
WF_COLOR = {"desktop": rb.OKABE_ITO["grey"], "quantmsdiann": rb.OKABE_ITO["blue"]}
SPECIES_COLOR = {"HUMAN": rb.OKABE_ITO["blue"], "YEAST": rb.OKABE_ITO["orange"],
                 "ECOLI": rb.OKABE_ITO["bluish_green"]}
FATE_CATS = ["1-2%", "2-5%", "not reported"]
FATE_LABELS = ["1-2%", "2-5%", "> 5% or not detected"]
FATE_COLORS = ["#BFD7EA", "#5B9BD5", "#D9D9D9"]
COLS = ["Run", "Precursor.Id", "Modified.Sequence", "Stripped.Sequence", "Protein.Group",
        "Protein.Names", "Decoy", "Precursor.Quantity", "Q.Value", "Lib.Q.Value",
        "PG.Q.Value", "Lib.PG.Q.Value"]


def report_path(ds: str, wf: str) -> Path:
    if wf == "desktop":
        return X / ds / "standalone/diann_report.parquet"
    return X / ds / "quantmsdiann/results/quant_tables/diann_report.parquet"


def load(ds: str, wf: str) -> pd.DataFrame:
    df = pd.read_parquet(report_path(ds, wf), columns=COLS)
    df = df[df["Decoy"] == 0].copy()
    # Run names differ by workflow (path-mangled vs plain); key on the file stem.
    df["run_key"] = df["Run"].map(lambda r: re.sub(r"^_data_raw_data_|_(raw|mzML|wiff)$", "", r))
    df["condition"] = np.where(df["run_key"].str.contains(DATASETS[ds]["cond_A"]), "A", "B")
    df["species"] = df["Protein.Names"].str.extract(r"_(HUMAN|YEAST|ECOLI)")[0]
    return df


def species_of(names: pd.Series) -> pd.Series:
    """Species of a precursor; ambiguous (multi-species) groups -> NaN."""
    sp = names.str.findall(r"_(HUMAN|YEAST|ECOLI)").map(lambda s: s[0] if len(set(s)) == 1 else None)
    return sp


# ---------------------------------------------------------------- statistics

def per_run_counts(df: pd.DataFrame) -> pd.DataFrame:
    prec = df[df["Q.Value"] <= Q]
    pg = df[df["PG.Q.Value"] <= Q]
    return pd.DataFrame({
        "precursors": prec.groupby("run_key")["Precursor.Id"].nunique(),
        "peptides": prec.groupby("run_key")["Stripped.Sequence"].nunique(),
        "protein_groups": pg.groupby("run_key")["Protein.Group"].nunique(),
    })


def global_sets(df: pd.DataFrame) -> dict[str, set]:
    prec = df[df["Lib.Q.Value"] <= Q]
    return {
        "Precursors": set(prec["Precursor.Id"]),
        "Modified peptides": set(prec["Modified.Sequence"]),
        "Peptide sequences": set(prec["Stripped.Sequence"]),
        "Protein groups": set(df.loc[df["Lib.PG.Q.Value"] <= Q, "Protein.Group"]),
    }


def specific_fate(own: pd.DataFrame, other: pd.DataFrame) -> pd.Series:
    """For precursors at Lib.Q <= 1% in `own` only: their best Lib.Q.Value in `other`."""
    own_ids = set(own.loc[own["Lib.Q.Value"] <= Q, "Precursor.Id"])
    other_best = other.groupby("Precursor.Id")["Lib.Q.Value"].min()
    only = own_ids - set(other_best[other_best <= Q].index)
    q = other_best.reindex(sorted(only))
    # Both reports stop at q = 5% (--qvalue 0.05), so "not reported" means q > 5% or undetected.
    bins = pd.cut(q, [Q, 0.02, 0.05], labels=["1-2%", "2-5%"], right=True)
    return bins.cat.add_categories("not reported").fillna("not reported").value_counts()


def quant_pairs(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """(run, precursor) quantities present at Q.Value <= 1% in both workflows."""
    key = ["run_key", "Precursor.Id"]
    qa = a[(a["Q.Value"] <= Q) & (a["Precursor.Quantity"] > 0)].groupby(key)["Precursor.Quantity"].max()
    qb = b[(b["Q.Value"] <= Q) & (b["Precursor.Quantity"] > 0)].groupby(key)["Precursor.Quantity"].max()
    m = pd.concat([np.log2(qa).rename("desktop"), np.log2(qb).rename("quantmsdiann")], axis=1, join="inner")
    return m


def precursor_matrix(df: pd.DataFrame) -> pd.DataFrame:
    d = df[(df["Q.Value"] <= Q) & (df["Precursor.Quantity"] > 0)]
    return d.pivot_table(index="Precursor.Id", columns="run_key", values="Precursor.Quantity", aggfunc="max")


def cvs(df: pd.DataFrame, keep: set | None = None) -> pd.DataFrame:
    mat = precursor_matrix(df)
    if keep is not None:
        mat = mat.loc[mat.index.intersection(keep)]
    cond = df.drop_duplicates("run_key").set_index("run_key")["condition"]
    out = {}
    for c in ("A", "B"):
        sub = mat[[r for r in mat.columns if cond.get(r) == c]].dropna()  # complete in all 3 reps
        out[c] = sub.std(axis=1, ddof=1) / sub.mean(axis=1) * 100
    return pd.DataFrame(out)


def log2fc(df: pd.DataFrame, keep: set | None = None) -> pd.DataFrame:
    mat = precursor_matrix(df)
    if keep is not None:
        mat = mat.loc[mat.index.intersection(keep)]
    cond = df.drop_duplicates("run_key").set_index("run_key")["condition"]
    a = mat[[r for r in mat.columns if cond.get(r) == "A"]]
    b = mat[[r for r in mat.columns if cond.get(r) == "B"]]
    ok = (a.notna().sum(axis=1) >= 2) & (b.notna().sum(axis=1) >= 2)
    fc = np.log2(a[ok].mean(axis=1) / b[ok].mean(axis=1))
    sp = species_of(df.drop_duplicates("Precursor.Id").set_index("Precursor.Id")["Protein.Names"]).reindex(fc.index)
    return pd.DataFrame({"log2fc": fc, "species": sp}).dropna()


# ---------------------------------------------------------------- plotting

def draw(ds: str, data: dict, axes) -> list[dict]:
    meta = DATASETS[ds]
    expected = rb.SPECIES_EXPECTED_LOG2_A_vs_B[meta["module"]]
    stats_rows: list[dict] = []
    ax_runs, ax_ovl, ax_fate, ax_ba, ax_cv, ax_fc = axes

    # (a) per-run depth
    levels = ["precursors", "peptides", "protein_groups"]
    lev_label = {"precursors": "Precursors", "peptides": "Peptides", "protein_groups": "Protein groups"}
    width = 0.38
    for j, wf in enumerate(WORKFLOWS):
        c = data["counts"][wf]
        means = [c[l].mean() for l in levels]
        x = np.arange(len(levels)) + (j - 0.5) * width
        ax_runs.bar(x, [m / 1e3 for m in means], width=width * 0.92, color=WF_COLOR[wf], label=WF_LABEL[wf])
        for xi, l in zip(x, levels):
            ax_runs.scatter(np.full(len(c), xi), c[l] / 1e3, s=10, color="#222222", zorder=3, linewidths=0)
        for l, m in zip(levels, means):
            stats_rows.append({"dataset": ds, "metric": f"per-run mean {l}", "workflow": wf, "value": round(m, 1)})
    for i, l in enumerate(levels):
        d, q = data["counts"]["desktop"][l].mean(), data["counts"]["quantmsdiann"][l].mean()
        top = max(data["counts"]["desktop"][l].max(), data["counts"]["quantmsdiann"][l].max()) / 1e3
        ax_runs.text(i, top * 1.04, f"{(q - d) / d * 100:+.1f}%", ha="center", va="bottom", fontsize=12)
    ax_runs.set_xticks(range(len(levels)), [lev_label[l] for l in levels])
    ax_runs.set_ylabel("Identified per run (thousands)")
    ax_runs.set_ylim(0, ax_runs.get_ylim()[1] * 1.08)

    # (b) global overlap per level: stacked shared / desktop-only / quantmsdiann-only
    sa, sb = data["sets"]["desktop"], data["sets"]["quantmsdiann"]
    names = list(sa)
    for i, n in enumerate(names):
        shared, only_d, only_q = len(sa[n] & sb[n]), len(sa[n] - sb[n]), len(sb[n] - sa[n])
        tot = shared + only_d + only_q
        parts = [(only_d, WF_COLOR["desktop"]), (shared, "#BFD7EA"), (only_q, WF_COLOR["quantmsdiann"])]
        left = 0.0
        for v, col in parts:
            ax_ovl.barh(i, v / tot * 100, left=left, color=col, height=0.62, edgecolor="white", linewidth=1.5)
            left += v / tot * 100
        ax_ovl.text(only_d / tot * 100 + shared / tot * 50, i, f"{shared / tot * 100:.1f}% shared",
                    ha="center", va="center", fontsize=12)
        stats_rows += [
            {"dataset": ds, "metric": f"global {n} shared", "workflow": "both", "value": shared},
            {"dataset": ds, "metric": f"global {n} only", "workflow": "desktop", "value": only_d},
            {"dataset": ds, "metric": f"global {n} only", "workflow": "quantmsdiann", "value": only_q},
            {"dataset": ds, "metric": f"global {n} Jaccard %", "workflow": "both", "value": round(shared / tot * 100, 2)},
        ]
    ax_ovl.set_yticks(range(len(names)), names)
    ax_ovl.invert_yaxis()
    ax_ovl.set_xlim(0, 100)
    ax_ovl.set_xlabel("% of union (global 1% FDR)")

    # (c) fate of workflow-specific precursors in the other workflow
    cats = FATE_CATS
    cat_col = FATE_COLORS
    for i, wf in enumerate(WORKFLOWS):
        fate = data["fate"][wf].reindex(cats).fillna(0)
        tot = fate.sum()
        left = 0.0
        for cat, col in zip(cats, cat_col):
            ax_fate.barh(i, fate[cat] / tot * 100, left=left, color=col, height=0.62, edgecolor="white", linewidth=1.5)
            left += fate[cat] / tot * 100
            stats_rows.append({"dataset": ds, "metric": f"{wf}-only precursors: other-workflow Lib.Q {cat}",
                               "workflow": wf, "value": int(fate[cat])})
        ax_fate.text(101, i, f"n={int(tot):,}", va="center", fontsize=12)
    ax_fate.set_yticks(range(2), ["Standalone\nonly", "quantmsdiann\nonly"])
    if ds == list(DATASETS)[0]:
        ax_fate.legend(handles=[Patch(color=c, label=l) for l, c in zip(FATE_LABELS, FATE_COLORS)],
                       loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, fontsize=12,
                       title="q-value in the other workflow", title_fontsize=12)
    ax_fate.invert_yaxis()
    ax_fate.set_xlim(0, 100)
    ax_fate.set_xlabel("% of workflow-specific precursors")

    # (d) Bland-Altman on shared (run, precursor) quantities
    m = data["pairs"]
    mean = (m["desktop"] + m["quantmsdiann"]) / 2
    diff = m["quantmsdiann"] - m["desktop"]
    bias, sd = diff.median(), diff.std()
    lo, hi = np.percentile(diff, [2.5, 97.5])
    ax_ba.hexbin(mean, diff, gridsize=70, bins="log", cmap="Blues", mincnt=1, linewidths=0)
    for y, ls in ((bias, "-"), (lo, "--"), (hi, "--")):
        ax_ba.axhline(y, color=rb.OKABE_ITO["vermillion"], lw=1.4, ls=ls)
    ax_ba.axhline(0, color="#444444", lw=0.8)
    ax_ba.set_ylim(-3, 3)
    ax_ba.set_xlabel("Mean log$_2$ precursor quantity")
    ax_ba.set_ylabel("log$_2$ quantmsdiann / standalone")
    within = (diff.abs() <= np.log2(1.5)).mean() * 100
    r = np.corrcoef(m["desktop"], m["quantmsdiann"])[0, 1]
    ax_ba.text(0.02, 0.97, f"median {bias:+.3f}; 95% in [{lo:.2f}, {hi:.2f}]\n"
               f"{within:.1f}% within 1.5-fold; r = {r:.3f}\nn = {len(m):,}",
               transform=ax_ba.transAxes, va="top", fontsize=11)
    stats_rows += [
        {"dataset": ds, "metric": "quant pairs n", "workflow": "both", "value": len(m)},
        {"dataset": ds, "metric": "quant median log2 diff (qmd - desktop)", "workflow": "both", "value": round(bias, 4)},
        {"dataset": ds, "metric": "quant log2 diff 2.5th pct", "workflow": "both", "value": round(lo, 3)},
        {"dataset": ds, "metric": "quant log2 diff 97.5th pct", "workflow": "both", "value": round(hi, 3)},
        {"dataset": ds, "metric": "quant % within 1.5-fold", "workflow": "both", "value": round(within, 2)},
        {"dataset": ds, "metric": "quant Pearson r (log2)", "workflow": "both", "value": round(r, 6)},
    ]

    # (e) CV distributions per condition, on precursors quantified by both
    pos = []
    for k, c in enumerate(("A", "B")):
        for j, wf in enumerate(WORKFLOWS):
            v = data["cv"][wf][c].dropna()
            p = k * 2.6 + j
            pos.append(p)
            parts = ax_cv.violinplot([v.clip(upper=100)], positions=[p], widths=0.85, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(WF_COLOR[wf]); b.set_alpha(0.75); b.set_edgecolor("none")
            med = v.median()
            ax_cv.hlines(med, p - 0.3, p + 0.3, color="#222222", lw=1.6)
            ax_cv.text(p, 103, f"{med:.1f}%", ha="center", va="bottom", fontsize=11)
            stats_rows.append({"dataset": ds, "metric": f"median precursor CV cond {c} (shared, complete)",
                               "workflow": wf, "value": round(med, 2)})
            stats_rows.append({"dataset": ds, "metric": f"precursors with CV<20% cond {c} (shared, complete)",
                               "workflow": wf, "value": int((v < 20).sum())})
    ax_cv.set_xticks([0.5, 3.1], ["Condition A", "Condition B"])
    ax_cv.set_ylabel("Precursor CV (%)")
    ax_cv.set_ylim(0, 112)

    # (f) species log2 A/B vs expected (ProteoBench), on precursors quantified by both
    sp_order = [s for s in ("HUMAN", "YEAST", "ECOLI") if s in expected]
    for k, sp in enumerate(sp_order):
        for j, wf in enumerate(WORKFLOWS):
            fc = data["fc"][wf]
            v = fc.loc[fc["species"] == sp, "log2fc"]
            p = k * 2.6 + j
            parts = ax_fc.violinplot([v.clip(-6, 6)], positions=[p], widths=0.85, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(WF_COLOR[wf]); b.set_alpha(0.75); b.set_edgecolor("none")
            ax_fc.hlines(v.median(), p - 0.3, p + 0.3, color="#222222", lw=1.6)
            eps = (v - expected[sp]).abs().median()
            stats_rows += [
                {"dataset": ds, "metric": f"{sp} median log2 A/B", "workflow": wf, "value": round(v.median(), 3)},
                {"dataset": ds, "metric": f"{sp} median |log2 A/B - expected|", "workflow": wf, "value": round(eps, 3)},
                {"dataset": ds, "metric": f"{sp} precursors in log2FC", "workflow": wf, "value": len(v)},
            ]
        ax_fc.hlines(expected[sp], k * 2.6 - 0.55, k * 2.6 + 1.55, color=rb.OKABE_ITO["vermillion"], lw=1.4, ls="--")
    ax_fc.set_xticks([k * 2.6 + 0.5 for k in range(len(sp_order))], [s.capitalize() for s in sp_order])
    ax_fc.set_ylabel("Precursor log$_2$(A/B)")
    return stats_rows


def identity_rows(ds: str) -> list[dict]:
    """Full-report comparison: rows matched on Run + Precursor.Id, every column compared."""
    a, b = (pd.read_parquet(report_path(ds, wf)) for wf in WORKFLOWS)
    m = a.merge(b, on=["Run", "Precursor.Id"], how="outer", suffixes=("_a", "_b"), indicator=True)
    rows = [{"dataset": ds, "metric": "report rows", "workflow": wf, "value": len(d)} for wf, d in zip(WORKFLOWS, (a, b))]
    rows.append({"dataset": ds, "metric": "report rows in one workflow only", "workflow": "both",
                 "value": int((m["_merge"] != "both").sum())})
    m = m[m["_merge"] == "both"]
    for c in (c for c in a.columns if c not in ("Run", "Precursor.Id")):
        x, y = m[f"{c}_a"], m[f"{c}_b"]
        diff = ~((x == y) | (x.isna() & y.isna()))
        if diff.any():
            rows.append({"dataset": ds, "metric": f"rows differing: {c}", "workflow": "both", "value": int(diff.sum())})
            rows.append({"dataset": ds, "metric": f"max abs difference: {c}", "workflow": "both",
                         "value": float((x[diff] - y[diff]).abs().max())})
    # Protein groups at Global.PG.Q.Value <= 1%, and the protein matrix (borderline groups move with PG q-values)
    for wf, d in zip(WORKFLOWS, (a, b)):
        rows.append({"dataset": ds, "metric": "protein groups Global.PG.Q.Value <= 0.01", "workflow": wf,
                     "value": d.loc[(d["Decoy"] == 0) & (d["Global.PG.Q.Value"] <= Q), "Protein.Group"].nunique()})
    pg = {wf: pd.read_csv(report_path(ds, wf).with_name("diann_report.pg_matrix.tsv"), sep="\t") for wf in WORKFLOWS}
    ids = {wf: set(m["Protein.Group"]) for wf, m in pg.items()}
    for wf, other in zip(WORKFLOWS, WORKFLOWS[::-1]):
        rows.append({"dataset": ds, "metric": "pg_matrix protein groups", "workflow": wf, "value": len(ids[wf])})
        rows.append({"dataset": ds, "metric": "pg_matrix protein groups in this workflow only", "workflow": wf,
                     "value": len(ids[wf] - ids[other])})
    run = lambda col: col.split("/")[-1].rsplit(".", 1)[0]
    qa, qb = (m.set_index("Protein.Group").rename(columns=run).iloc[:, 5:] for m in (pg["desktop"], pg["quantmsdiann"]))
    shared = sorted(ids["desktop"] & ids["quantmsdiann"])
    va, vb = qa.loc[shared, qb.columns].to_numpy(), qb.loc[shared].to_numpy()
    both = ~np.isnan(va) & ~np.isnan(vb)
    rows.append({"dataset": ds, "metric": "pg_matrix shared quantities differing", "workflow": "both",
                 "value": int((va[both] != vb[both]).sum())})
    return rows


def main() -> None:
    rb.apply_house_style()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "data").mkdir(exist_ok=True)
    fig, axes = plt.subplots(len(DATASETS), 6, figsize=(36, 6.6 * len(DATASETS)),
                             gridspec_kw={"width_ratios": [1.15, 1.25, 1.1, 1.2, 1, 1], "wspace": 0.42, "hspace": 0.55})
    rows = []
    for r, ds in enumerate(DATASETS):
        dfs = {wf: load(ds, wf) for wf in WORKFLOWS}
        pairs = quant_pairs(dfs["desktop"], dfs["quantmsdiann"])
        shared_prec = set(pairs.index.get_level_values("Precursor.Id"))
        data = {
            "counts": {wf: per_run_counts(dfs[wf]) for wf in WORKFLOWS},
            "sets": {wf: global_sets(dfs[wf]) for wf in WORKFLOWS},
            "fate": {"desktop": specific_fate(dfs["desktop"], dfs["quantmsdiann"]),
                     "quantmsdiann": specific_fate(dfs["quantmsdiann"], dfs["desktop"])},
            "pairs": pairs,
            "cv": {wf: cvs(dfs[wf], shared_prec) for wf in WORKFLOWS},
            "fc": {wf: log2fc(dfs[wf], shared_prec) for wf in WORKFLOWS},
        }
        for wf in WORKFLOWS:
            data["counts"][wf].assign(dataset=ds, workflow=wf).to_csv(
                OUT / "data" / f"{ds}_{wf}_per_run_counts.tsv", sep="\t")
        rows += draw(ds, data, axes[r])
        rows += identity_rows(ds)
        axes[r][0].annotate(DATASETS[ds]["label"], xy=(0, 1.2), xycoords="axes fraction",
                            fontsize=19, fontweight="bold", ha="left")
    titles = ["Per-run depth (1% run FDR)", "Global overlap by level", "Workflow-specific precursors",
              "Quantity agreement (shared)", "Precision (shared precursors)", "Accuracy vs expected ratio"]
    for c, t in enumerate(titles):
        axes[0][c].set_title(f"{'abcdef'[c]}  {t}", loc="left")
    handles = [Patch(color=WF_COLOR[w], label=WF_LABEL[w]) for w in WORKFLOWS] + [
        Patch(color="#BFD7EA", label="shared"),
        Line2D([0], [0], color=rb.OKABE_ITO["vermillion"], ls="--", label="expected ratio (f) / 95% range (d)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.03), fontsize=15)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"supp_desktop_vs_quantmsdiann.{ext}")
    pd.DataFrame(rows).to_csv(OUT / "data/desktop_vs_quantmsdiann_stats.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
