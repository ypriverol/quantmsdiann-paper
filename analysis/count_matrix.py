#!/usr/bin/env python
"""Reproducible row counting for DIA-NN ``*_pg_matrix.tsv`` / ``*_pr_matrix.tsv``.

Per methods.md §1 ("Original / deposited side of reanalysis comparisons"),
DIA-NN matrix files are ALREADY q-filtered count matrices and carry no
q-value columns, so there is no q-filter left to apply. The single
reproducible number is the count of matrix ROWS that are quantified:

  * for a ``*_pg_matrix.tsv`` a row is one protein group;
  * for a ``*_pr_matrix.tsv`` a row is one precursor.

Counting rule (no other filter is admissible):
  * NO contaminant/target filter -- ``CONTAM_`` / ``Cont_`` / ``ENTRAP_``
    prefixed rows are counted like any other;
  * NO positive-quantity filter -- a literal ``0`` quantity is "quantified"
    and is counted; only an empty / ``NA`` cell means "not measured";
  * a row counts iff it has at least one non-empty quantity across the
    sample columns (matrices can contain fully-empty rows);
  * decoys are not present in matrices, so there is nothing to drop.

The sample columns are every column that is not one of the fixed DIA-NN
metadata columns (``PG_METADATA`` / ``PR_METADATA`` below).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# Fixed DIA-NN metadata (non-sample) columns. Everything else in the header is
# a per-run / per-sample quantity column.
# Non-sample columns to EXCLUDE when present. DIA-NN versions differ slightly
# (e.g. the deposited 1.8.1 pg_matrix carries ``Protein.Ids`` but the 2.5.1
# pg_matrix does not), so this is an "exclude if present" set, not a strict
# schema. Any header column not in here is treated as a sample/quantity column.
PG_METADATA = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences",
]
PR_METADATA = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
]

# Anchor columns that MUST be present, so a wrong-kind / mis-parsed file fails
# loudly instead of silently mis-counting.
_PG_REQUIRED = ["Protein.Group", "First.Protein.Description"]
_PR_REQUIRED = ["Protein.Group", "Precursor.Id", "First.Protein.Description"]


def count_matrix_rows(matrix_path: Path, metadata_cols: list[str]) -> int:
    """Number of quantified rows in a DIA-NN pg/pr count matrix.

    A row is "quantified" iff at least one of its sample columns (every column
    not in ``metadata_cols`` that is actually present in the matrix) holds a
    non-empty, non-``NA`` value. Zeros count; contaminant/entrapment rows count.

    Raises ``ValueError`` if a required anchor column for the chosen layout is
    missing from the header (guards against mis-parsed / wrong-kind matrices).
    """
    required = _PR_REQUIRED if "Precursor.Id" in metadata_cols else _PG_REQUIRED
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"matrix {matrix_path} is missing metadata column(s): "
            f"{', '.join(missing)}"
        )
    sample_cols = [c for c in df.columns if c not in metadata_cols]
    if not sample_cols:
        return 0
    # dtype=str -> empty cells are read as NaN; the literal text "NA" must also
    # be treated as not-measured (DIA-NN writes empty, but be defensive).
    samples = df[sample_cols].replace({"NA": pd.NA, "": pd.NA})
    quantified = samples.notna().any(axis=1)
    return int(quantified.sum())
