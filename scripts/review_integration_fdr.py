"""Reviewer response (R1 Q1.9): protein-level FDR when the cell-line cohorts are
integrated, using the quantms criteria (Dai et al. 2024, Nat Methods).

Per cohort (one quantmsdiann run each, FDR controlled over the whole cohort):
protein groups at Lib.PG.Q.Value <= 0.01, with and without the quantms rule of
>= 2 unique peptides (distinct proteotypic Stripped.Sequence at Lib.Q.Value
<= 0.01). Groups are counted per accession (every member of a group),
contaminants dropped.

Integration (quantms stand-in decoy procedure): in each cohort list, sorted by
q-value, stand-in decoys are injected so that the decoy/target ratio at every
target equals its q-value. Lists are merged, each accession keeps its lowest
q-value, and q-values are re-estimated as decoys / targets scoring as well or
better (made monotone). Proteins at adjusted q < 0.01 are kept. ENTRAP_
accessions (PXD017199, PXD041421 only) are kept as targets and tracked.

Entrapment control: the q-value cut-off is further lowered until the
entrapment-estimated FDP (combined estimator, scripts.review_entrapment_fdp)
of the merged list is <= 0.01, counted over proteins seen in a cohort searched
with entrapment.

Usage: python -m scripts.review_integration_fdr
"""
from __future__ import annotations

import shutil

import fsspec
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from scripts.review_entrapment_fdp import CONTAM, DATA, entrapment_ratio, fdp

FTP = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes"
REPORTS = {
    "PXD003539": f"{FTP}/quantmsdiann-benchmarks/cell-lines/PXD003539/v2_5_1_enterprise/quant_tables/diann_report.parquet",
    "PXD004701": f"{FTP}/quantmsdiann-benchmarks/cell-lines/PXD004701/v2_5_1/quant_tables/diann_report.parquet",
    "PXD030304": f"{FTP}/quantmsdiann-benchmarks/cell-lines/PXD030304/v2_5_1/quant_tables/diann_report.parquet",
    "PXD017199": f"{FTP}/quantms-collections/absolute-expression-2.0/cell-lines/PXD017199/quant_tables/diann_report.parquet",
    "PXD041421": f"{FTP}/quantms-collections/absolute-expression-2.0/cell-lines/PXD041421/quant_tables/diann_report.parquet",
}
COLS = ["Protein.Group", "Stripped.Sequence", "Proteotypic", "Lib.Q.Value", "Lib.PG.Q.Value"]
OUT = DATA.parent / "analysis/figures/review/entrapment/data/integration_fdr.tsv"
# Tracked in git (entrapment/data/ is not): summary tables, figure and README for the response.
TRACKED = DATA.parent / "analysis/figures/review/integration"
Q = 0.01
MIN_PEPTIDES = 2
ENTRAP_SEARCHED = ("PXD017199", "PXD041421")


def cohort_proteins(ds: str) -> pd.DataFrame:
    """Accession, q-value and number of unique peptides for protein groups at 1%
    protein FDR (cached in OUT's folder as integration_<ds>_all.tsv)."""
    cache = OUT.with_name(f"integration_{ds}_all.tsv")
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    local = DATA / ds / "diann_report.parquet"
    src = local.open("rb") if local.exists() else fsspec.open(REPORTS[ds], block_size=2**22).open()
    pf = pq.ParquetFile(src)
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


def entrapment_cutoff(merged: pd.DataFrame, r: float) -> tuple[float, pd.Series]:
    """Largest q-value cut-off at which the entrapment-estimated FDP of the merged
    list is <= Q, counted over proteins seen in a cohort searched with entrapment."""
    m = merged[merged["datasets"].str.contains("|".join(ENTRAP_SEARCHED))].sort_values("q")
    is_e = m["accession"].str.startswith("ENTRAP_").to_numpy()
    n_e, n_t = np.cumsum(is_e), np.cumsum(~is_e)
    ok = fdp(n_t, n_e, r) <= Q
    t = m["q"].to_numpy()[ok].max()
    return t, merged["q"] <= t


def summarise(label: str, lists: list[pd.DataFrame], r: float) -> list[dict]:
    merged = integrate(lists)
    members = pd.concat(lists).groupby("accession")["dataset"].agg(lambda d: ";".join(sorted(d)))
    merged["datasets"] = merged["accession"].map(members)
    merged.to_csv(OUT.with_name(f"integration_proteins_{label}.tsv"), sep="\t", index=False)
    entrap = merged["accession"].str.startswith("ENTRAP_")
    searched = merged["datasets"].str.contains("|".join(ENTRAP_SEARCHED))

    def row(name, n, keep):
        n_e = int((entrap & keep).sum())
        n_t = int((~entrap & keep & searched).sum())
        return {"rule": label, "set": name, "proteins": n, "entrapment": n_e,
                "entrapment_fdp": round(fdp(n_t, n_e, r), 4)}

    rows = []
    for x in lists:
        e = x["accession"].str.startswith("ENTRAP_")
        n_e = int(e.sum())
        rows.append({"rule": label, "set": x["dataset"].iloc[0], "proteins": len(x), "entrapment": n_e,
                     "entrapment_fdp": round(fdp(len(x) - n_e, n_e, r), 4) if x["dataset"].iloc[0] in ENTRAP_SEARCHED else None})
    everything = merged["q"] >= 0
    rows.append(row("merged", len(merged), everything))
    decoy = merged["q_adj"] < Q
    rows.append(row("decoy_1pct", int(decoy.sum()), decoy))
    t, entrap_ok = entrapment_cutoff(merged, r)
    both = decoy & entrap_ok
    rows.append(row(f"decoy_and_entrapment_1pct (q<={t:.5f})", int(both.sum()), both))
    by_project(label, merged, ~decoy, ~entrap_ok).to_csv(
        OUT.with_name(f"integration_removed_by_project_{label}.tsv"), sep="\t", index=False)
    return rows


def by_project(label: str, merged: pd.DataFrame, dec: pd.Series, ent: pd.Series) -> pd.DataFrame:
    """Proteins of each cohort removed by the decoy step and, after it, by the entrapment step."""
    ent = ent & ~dec
    removed = dec | ent
    is_e = merged["accession"].str.startswith("ENTRAP_")
    ex = merged.assign(ds=merged["datasets"].str.split(";")).explode("ds")
    rows = []
    for ds, g in ex.groupby("ds"):
        i = g.index
        rows.append({"rule": label, "cohort": ds, "proteins": len(g), "removed_decoy": int(dec[i].sum()),
                     "removed_entrapment": int(ent[i].sum()), "removed_total": int(removed[i].sum()),
                     "removed_only_in_this_cohort": int((removed[i] & (g["datasets"] == ds)).sum()),
                     "kept": int((~removed[i]).sum()), "entrapment_in": int(is_e[i].sum()),
                     "entrapment_kept": int((is_e[i] & ~removed[i]).sum())})
    return pd.DataFrame(rows)


def main() -> int:
    r = entrapment_ratio()
    lists = []
    for ds in REPORTS:
        lst = cohort_proteins(ds)
        print(ds, len(lst), flush=True)
        lists.append(lst)
    rows = summarise("1pct_protein_fdr", lists, r)
    rows += summarise(f">={MIN_PEPTIDES}_unique_peptides", [x[x["n_pep"] >= MIN_PEPTIDES] for x in lists], r)
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
