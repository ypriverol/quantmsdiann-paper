#!/usr/bin/env python
"""Recompute the staged per-cohort reanalysis protein-group count JSONs under
the methods.md §1 rule (no contaminant/target filter, no Global.* gate).

For each bulk cell-line reanalysis cohort whose figure script reads a small
``diann_report_protein_counts.json`` cache, this regenerates that cache from
the DIA-NN long-format report (parquet) using ONLY the admissible global
filters:

  * ``prot_global`` -- distinct ``Protein.Group`` with ``Lib.PG.Q.Value <= 0.01``
  * ``prot_2pep``   -- global protein groups with >= 2 distinct stripped
    peptides among ``Lib.Q.Value <= 0.01`` precursors

Decoys (``Decoy == 1``) are dropped; nothing else is filtered (zeros counted,
contaminant/entrapment rows counted). The report is read from the public PRIDE
FTP (streamed with column projection) unless a local parquet is supplied.

This is the reproducible provenance for the curated JSON caches consumed by
``figure_pxd004701_sun_vs_quantmsdiann`` and
``figure_pxd030304_procan_vs_quantmsdiann``; PXD003539 (NCI-60) computes its
counts inline from the local parquet via
``figure_original_vs_quantmsdiann.report_global_counts_diann``.

Usage:
    python -m analysis.recompute_reanalysis_pg_counts PXD004701
    python -m analysis.recompute_reanalysis_pg_counts PXD030304 --parquet /path/local.parquet
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

FTP = ("https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
       "quantmsdiann-benchmarks/cell-lines")
REPORTS = {
    "PXD004701": f"{FTP}/PXD004701/v2_5_1/quant_tables/diann_report.parquet",
    "PXD030304": f"{FTP}/PXD030304/v2_5_1/quant_tables/diann_report.parquet",
}


def recompute(parquet_source: str | Path, *, qvalue_cutoff: float = 0.01,
              batch_size: int = 2_000_000) -> dict[str, int]:
    """Return ``{"prot_global": int, "prot_2pep": int}`` under methods.md §1."""
    import pyarrow.parquet as pq

    source = str(parquet_source)
    if source.startswith(("http://", "https://")):
        import fsspec
        opener = lambda: fsspec.filesystem("https").open(source, "rb")
    else:
        opener = lambda: open(source, "rb")

    cols = ["Protein.Group", "Stripped.Sequence",
            "Lib.Q.Value", "Lib.PG.Q.Value", "Decoy"]
    lib_pg: set[str] = set()
    peps_pg: dict[str, set[str]] = collections.defaultdict(set)
    with opener() as fh:
        pf = pq.ParquetFile(fh)
        have = set(pf.schema_arrow.names)
        use = [c for c in cols if c in have]
        for batch in pf.iter_batches(batch_size=batch_size, columns=use):
            d = batch.to_pydict()
            n = len(d["Protein.Group"])
            dec = d.get("Decoy", [0] * n)
            for i in range(n):
                if dec[i] == 1:
                    continue
                lpg = d["Lib.PG.Q.Value"][i]
                pg = d["Protein.Group"][i]
                if lpg is not None and lpg <= qvalue_cutoff:
                    lib_pg.add(pg)
                    lq = d["Lib.Q.Value"][i]
                    if lq is not None and lq <= qvalue_cutoff:
                        peps_pg[pg].add(d["Stripped.Sequence"][i])
    prot_2pep = sum(1 for s in peps_pg.values() if len(s) >= 2)
    return {"prot_global": len(lib_pg), "prot_2pep": prot_2pep}


def main(argv: list[str] | None = None) -> int:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("accession", choices=sorted(REPORTS))
    ap.add_argument("--parquet", type=Path, default=None,
                    help="local parquet path (defaults to the PRIDE FTP URL)")
    args = ap.parse_args(argv)

    source = args.parquet or REPORTS[args.accession]
    print(f"Streaming {args.accession} report from {source} ...", flush=True)
    counts = recompute(source)
    out = REPO_ROOT / "data" / args.accession / "diann_report_protein_counts.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(counts))
    print(f"  prot_global={counts['prot_global']:,}  "
          f"prot_2pep={counts['prot_2pep']:,}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
