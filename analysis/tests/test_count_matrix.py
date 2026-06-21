"""Tests for analysis.count_matrix.count_matrix_rows.

Rule under test (methods.md §1, "Original / deposited side"):
DIA-NN ``*_pg_matrix.tsv`` / ``*_pr_matrix.tsv`` are already q-filtered count
matrices with no q-value columns, so there is no q-filter to apply. The
reproducible number is the count of matrix ROWS (protein groups for a pg
matrix; precursors for a pr matrix) that are quantified.

A row counts iff it has at least one non-empty quantity across the sample
columns. There is NO contaminant/target filter and NO positive-quantity
filter: a zero quantity is "quantified" (counted); only an empty/NA cell is
"not measured". Decoys are not present in matrices.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from analysis.count_matrix import PG_METADATA, PR_METADATA, count_matrix_rows


def _write(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "matrix.tsv"
    p.write_text(textwrap.dedent(body).lstrip("\n"))
    return p


def test_pg_matrix_counts_rows_with_at_least_one_non_na(tmp_path: Path) -> None:
    matrix = _write(
        tmp_path,
        """
        Protein.Group\tProtein.Names\tGenes\tFirst.Protein.Description\tN.Sequences\tN.Proteotypic.Sequences\tRun_A\tRun_B\tRun_C
        P1\tA\tG1\td1\t3\t3\t10\t\t20
        P2\tB\tG2\td2\t1\t1\t\t\t
        P3\tC\tG3\td3\t2\t2\t\tNA\t
        P4\tD\tG4\td4\t5\t5\t\t5\t
        """,
    )
    # P1 and P4 have one real value; P2 all empty; P3 only "NA" -> 2.
    assert count_matrix_rows(matrix, PG_METADATA) == 2


def test_zero_quantity_is_counted_contaminants_not_filtered(tmp_path: Path) -> None:
    matrix = _write(
        tmp_path,
        """
        Protein.Group\tProtein.Names\tGenes\tFirst.Protein.Description\tN.Sequences\tN.Proteotypic.Sequences\tRun_A\tRun_B
        P1\tA\tG1\td1\t3\t3\t0\t
        CONTAM_P2\tB\tG2\td2\t1\t1\t5\t
        ENTRAP_P3\tC\tG3\td3\t2\t2\t\t0.0
        """,
    )
    # All three rows quantified: zero counts, contaminant + entrapment prefixes
    # are NOT filtered out.
    assert count_matrix_rows(matrix, PG_METADATA) == 3


def test_pr_matrix_counts_precursor_rows(tmp_path: Path) -> None:
    matrix = _write(
        tmp_path,
        """
        Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\tStripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\tRun_A\tRun_B
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t2\tAAAR2\t100\t
        P1\tP1\tA\tG1\td\t1\tAAAR\tAAAR\t3\tAAAR3\t\t
        P2\tP2\tB\tG2\td\t1\tBBBR\tBBBR\t2\tBBBR2\t\t200
        """,
    )
    # AAAR2 and BBBR2 quantified; AAAR3 has no value -> 2 precursor rows.
    assert count_matrix_rows(matrix, PR_METADATA) == 2


def test_raises_if_metadata_column_missing(tmp_path: Path) -> None:
    matrix = _write(
        tmp_path,
        """
        Protein.Group\tProtein.Names\tGenes\tN.Sequences\tN.Proteotypic.Sequences\tRun_A
        P1\tA\tG1\t3\t3\t10
        """,
    )
    with pytest.raises(ValueError, match="First.Protein.Description"):
        count_matrix_rows(matrix, PG_METADATA)
