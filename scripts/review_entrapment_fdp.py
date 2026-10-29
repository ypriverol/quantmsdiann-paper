"""Reviewer response (R1 Q1.9): entrapment-estimated FDP for the two cohorts
searched against the human SwissProt + entrapment FASTA (PXD041421, PXD017199;
DIA-NN 2.5.0, one quantmsdiann run per cohort).

Inputs are the DIA-NN matrices (already filtered at matrix q-value 0.01).
Entrapment entries carry the ENTRAP_ prefix; contaminants (CONTAM_/Cont_) are
dropped. A row counts as entrapment when all its accessions are ENTRAP_ and as
target when none is. FDP uses the combined estimator (Wen et al. 2025):
N_e (1 + 1/r) / (N_t + N_e), with r = entrapment / target FASTA entries.

Reported per cohort and level (precursor, protein group): cohort-wide FDP,
median and range of per-run FDP, and FDP by the number of runs in which the
identification was quantified (1, 2-5, >=6).

Usage: python -m scripts.review_entrapment_fdp
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FASTA = ROOT / "sdrf/fasta/Homo-sapiens-uniprot-reviewed-entrap-contaminants-202605.fasta.gz"
OUT = ROOT / "analysis/figures/review/entrapment/data/entrapment_fdp.tsv"
COHORTS = ("PXD041421", "PXD017199")
LEVELS = {"precursor": "diann_report.pr_matrix.tsv", "protein_group": "diann_report.pg_matrix.tsv"}
CONTAM = ("CONTAM_", "Cont_")
RUN_BINS = (("1", 1, 1), ("2-5", 2, 5), (">=6", 6, None))


def entrapment_ratio() -> float:
    n_e = n_t = 0
    with gzip.open(FASTA, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            if "ENTRAP_" in line:
                n_e += 1
            elif not any(c in line for c in CONTAM):
                n_t += 1
    return n_e / n_t


def fdp(n_t, n_e, r: float):
    return n_e * (1 + 1 / r) / (n_t + n_e)


def run_columns(df: pd.DataFrame, n_runs: int) -> list[str]:
    # run columns are the trailing n_runs columns of a DIA-NN matrix
    return list(df.columns[-n_runs:])


def main() -> int:
    r = entrapment_ratio()
    rows = []
    for ds in COHORTS:
        n_runs = len(pd.read_csv(DATA / ds / f"{ds}.sdrf.tsv", sep="\t"))
        for level, fname in LEVELS.items():
            df = pd.read_csv(DATA / ds / fname, sep="\t", low_memory=False)
            acc = df["Protein.Group"].astype(str).str.split(";")
            is_contam = acc.apply(lambda a: any(x.startswith(CONTAM) for x in a))
            is_e = acc.apply(lambda a: all("ENTRAP_" in x for x in a)) & ~is_contam
            is_t = acc.apply(lambda a: not any("ENTRAP_" in x for x in a)) & ~is_contam
            seen = df[run_columns(df, n_runs)].notna()
            per_run = fdp(seen[is_t].sum(), seen[is_e].sum(), r)
            n_seen = seen.sum(axis=1)
            row = {
                "dataset": ds, "level": level, "runs": n_runs, "entrapment_ratio": round(r, 4),
                "target": int(is_t.sum()), "entrapment": int(is_e.sum()),
                "fdp_cohort": fdp(is_t.sum(), is_e.sum(), r),
                "fdp_per_run_median": per_run.median(),
                "fdp_per_run_min": per_run.min(), "fdp_per_run_max": per_run.max(),
            }
            for label, lo, hi in RUN_BINS:
                m = (n_seen >= lo) & (n_seen <= hi if hi else True)
                n_t, n_e = int((is_t & m).sum()), int((is_e & m).sum())
                row[f"target_in_{label}_runs"] = n_t
                row[f"entrapment_in_{label}_runs"] = n_e
                row[f"fdp_in_{label}_runs"] = fdp(n_t, n_e, r) if n_t + n_e else float("nan")
            rows.append(row)
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print(out.T.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
