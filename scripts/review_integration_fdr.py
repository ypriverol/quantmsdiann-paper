"""Reviewer response (R1 Q1.9): protein-level FDR when the cell-line cohorts are
integrated, using the quantms criteria (Dai et al. 2024, Nat Methods).

All five cohorts come from the harmonised DIA-NN 2.5.1 runs (scripts.cellline_reports:
same pipeline, FASTA and parameters, no entrapment).

Per cohort (one quantmsdiann run each, FDR controlled over the whole cohort):
protein groups at Lib.PG.Q.Value <= 0.01, with and without the quantms rule of
>= 2 unique peptides (distinct proteotypic Stripped.Sequence at Lib.Q.Value
<= 0.01). Groups are counted per accession (every member of a group),
contaminants dropped.

Integration (quantms stand-in decoy procedure): in each cohort list, sorted by
q-value, stand-in decoys are injected so that the decoy/target ratio at every
target equals its q-value. Lists are merged, each accession keeps its lowest
q-value, and q-values are re-estimated as decoys / targets scoring as well or
better (made monotone). Proteins at adjusted q < 0.01 are kept.

Usage: python -m scripts.review_integration_fdr
"""
from __future__ import annotations

import shutil

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from scripts.cellline_reports import COHORTS, LOCAL, open_report, source

COLS = ["Protein.Group", "Stripped.Sequence", "Proteotypic", "Lib.Q.Value", "Lib.PG.Q.Value"]
CONTAM = ("CONTAM_", "Cont_")
OUT = LOCAL / "integration" / "integration_fdr.tsv"
# Tracked in git: summary tables, figure and README for the response.
TRACKED = LOCAL.parents[1] / "analysis/figures/review/integration"
Q = 0.01
MIN_PEPTIDES = 2


def cohort_proteins(ds: str) -> pd.DataFrame:
    """Accession, q-value and number of unique peptides for protein groups at 1%
    protein FDR (cached next to OUT as integration_<ds>.tsv)."""
    cache = OUT.with_name(f"integration_{ds}.tsv")
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    pf = pq.ParquetFile(open_report(ds))
    pg_q, peps = {}, set()
    for i in range(pf.num_row_groups):
        t = pf.read_row_group(i, columns=COLS).to_pandas()
        t = t[t["Lib.PG.Q.Value"] <= Q]
        g = t.groupby("Protein.Group")["Lib.PG.Q.Value"].min()
        for pg, q in g.items():
            pg_q[pg] = min(q, pg_q.get(pg, 1.0))
        u = t[(t["Proteotypic"] == 1) & (t["Lib.Q.Value"] <= Q)]
        peps.update(zip(u["Protein.Group"], u["Stripped.Sequence"]))
    n_pep = pd.Series([pg for pg, _ in peps]).value_counts()
    df = pd.DataFrame({"Protein.Group": list(pg_q), "q": list(pg_q.values())})
    df["n_pep"] = df["Protein.Group"].map(n_pep).fillna(0).astype(int)
    df = df.assign(accession=df["Protein.Group"].str.split(";")).explode("accession")
    df = df[~df["accession"].str.startswith(CONTAM)]
    out = df.groupby("accession", as_index=False).agg(q=("q", "min"), n_pep=("n_pep", "max")).assign(dataset=ds)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cache, sep="\t", index=False)
    return out


def stand_in_decoys(lst: pd.DataFrame) -> pd.DataFrame:
    """Decoys with scores such that decoys/targets at each target equals its q-value."""
    q = np.sort(lst["q"].to_numpy())
    n_decoys = np.floor(q * np.arange(1, len(q) + 1)).astype(int)
    new = np.diff(np.concatenate([[0], np.maximum.accumulate(n_decoys)]))
    return pd.DataFrame({"q": np.repeat(q, new)})


def integrate(lists: list[pd.DataFrame]) -> pd.DataFrame:
    targets = pd.concat(lists).groupby("accession", as_index=False)["q"].min()
    decoys = pd.concat([stand_in_decoys(x) for x in lists])
    scores = np.sort(targets["q"].to_numpy())
    n_t = np.searchsorted(scores, targets["q"], side="right")
    n_d = np.searchsorted(np.sort(decoys["q"].to_numpy()), targets["q"], side="right")
    targets["q_adj"] = n_d / n_t
    order = targets["q"].argsort()[::-1]
    targets.loc[targets.index[order], "q_adj"] = np.minimum.accumulate(targets["q_adj"].to_numpy()[order])
    return targets


def summarise(label: str, lists: list[pd.DataFrame]) -> list[dict]:
    merged = integrate(lists)
    members = pd.concat(lists).groupby("accession")["dataset"].agg(lambda d: ";".join(sorted(d)))
    merged["datasets"] = merged["accession"].map(members)
    merged.to_csv(OUT.with_name(f"integration_proteins_{label}.tsv"), sep="\t", index=False)
    rows = [{"rule": label, "set": x["dataset"].iloc[0], "proteins": len(x)} for x in lists]
    rows.append({"rule": label, "set": "merged", "proteins": len(merged)})
    kept = merged["q_adj"] < Q
    rows.append({"rule": label, "set": "decoy_1pct", "proteins": int(kept.sum())})
    by_project(label, merged, ~kept).to_csv(
        OUT.with_name(f"integration_removed_by_project_{label}.tsv"), sep="\t", index=False)
    return rows


def by_project(label: str, merged: pd.DataFrame, removed: pd.Series) -> pd.DataFrame:
    """Proteins of each cohort removed by the decoy re-estimation."""
    ex = merged.assign(ds=merged["datasets"].str.split(";")).explode("ds")
    rows = []
    for ds, g in ex.groupby("ds"):
        i = g.index
        rows.append({"rule": label, "cohort": ds, "proteins": len(g), "removed": int(removed[i].sum()),
                     "removed_only_in_this_cohort": int((removed[i] & (g["datasets"] == ds)).sum()),
                     "kept": int((~removed[i]).sum())})
    return pd.DataFrame(rows)


def main() -> int:
    lists = []
    for ds in COHORTS:
        lst = cohort_proteins(ds)
        print(ds, len(lst), source(ds), flush=True)
        lists.append(lst)
    rows = summarise("1pct_protein_fdr", lists)
    rows += summarise(f">={MIN_PEPTIDES}_unique_peptides", [x[x["n_pep"] >= MIN_PEPTIDES] for x in lists])
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    TRACKED.mkdir(parents=True, exist_ok=True)
    for f in [OUT, *sorted(OUT.parent.glob("integration_removed_by_project_*.tsv"))]:
        shutil.copy2(f, TRACKED / f.name.replace(">=", "min"))
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
