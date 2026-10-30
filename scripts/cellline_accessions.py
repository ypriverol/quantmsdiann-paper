"""Protein accessions of the five cell-line cohorts under one rule (methods.md §1), shared
by the pan-cohort atlas (Fig. S6, the pan-cohort union) and the integration analysis
(R1.9, R1.12).

- Cohort set: every accession of the protein groups at Lib.PG.Q.Value <= 0.01 in the
  cohort's DIA-NN report, contaminant accessions dropped (scripts.review_integration_fdr
  cohort_proteins, the same lists the integration analysis merges).
- Per run: the cohort-set accessions of the protein groups at PG.Q.Value <= 0.01 in that
  run (methods.md §1 per-run rule); runs map to cell lines through the as-run SDRF
  (comment[data file]).

Reports come from the harmonised DIA-NN 2.5.1 runs (scripts.cellline_reports).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from scripts.cellline_reports import COHORTS, open_report
from scripts.review_integration_fdr import OUT, Q, cohort_proteins


def cohort_accessions(ds: str) -> set[str]:
    return set(cohort_proteins(ds)["accession"])


def _stem(name: str) -> str:
    name = Path(str(name)).name
    for ext in (".wiff.scan", ".wiff", ".raw", ".mzML", ".d", ".dia"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def accessions_per_run(ds: str) -> dict[str, set[str]]:
    """Run -> cohort-set accessions of the protein groups at PG.Q.Value <= 0.01 in that run
    (cached as per_run_<ds>.json next to the integration tables)."""
    cache = OUT.with_name(f"per_run_{ds}.json")
    if cache.exists():
        return {r: set(a) for r, a in json.loads(cache.read_text()).items()}
    keep = cohort_accessions(ds)
    pf = pq.ParquetFile(open_report(ds))
    pairs = set()
    for i in range(pf.num_row_groups):
        t = pf.read_row_group(i, columns=["Run", "Protein.Group", "PG.Q.Value"]).to_pandas()
        t = t[t["PG.Q.Value"] <= Q]
        pairs.update(zip(t["Run"], t["Protein.Group"]))
    out: dict[str, set[str]] = {}
    for run, pg in pairs:
        out.setdefault(_stem(run), set()).update(a for a in pg.split(";") if a in keep)
    cache.write_text(json.dumps({r: sorted(a) for r, a in out.items()}))
    return out


def accessions_per_cell_line(ds: str, sdrf: Path, normalise) -> dict[str, set[str]]:
    """Normalised cell line -> union of the accessions of its runs."""
    s = pd.read_csv(sdrf, sep="\t", dtype=str)
    run_cell = {_stem(f): normalise(c) for c, f in zip(s["characteristics[cell line]"], s["comment[data file]"])
                if isinstance(f, str) and normalise(c)}
    out: dict[str, set[str]] = {}
    for run, accs in accessions_per_run(ds).items():
        cell = run_cell.get(run)
        if cell:
            out.setdefault(cell, set()).update(accs)
    return out


if __name__ == "__main__":
    sets = {ds: cohort_accessions(ds) for ds in COHORTS}
    union = set().union(*sets.values())
    for ds, s in sets.items():
        print(ds, len(s))
    print("union", len(union), "in all five", len(set.intersection(*sets.values())))
