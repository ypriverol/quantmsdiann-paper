"""Tests for analysis.count_report_ids.count_report (inline fixtures, no I/O).

Filter rule (Vadim review, 2026-06-21), enforced by these tests:
  * per-run protein groups: PG.Q.Value only -- NO contaminant/target filter,
    NO global filter, zeros counted.
  * per-run / replicate precursors: Q.Value only -- NO Global.Q.Value.
  * global protein groups: Lib.PG.Q.Value only.
  * global precursors: Lib.Q.Value only.
  * decoys (Decoy == 1) are dropped (FDR null model, not a "filter").
"""
from __future__ import annotations

import pandas as pd


def _row(run, pid, pg, *, q=0.001, pgq=0.001, lib_q=0.001, lib_pgq=0.001,
         decoy=0, seq=None, gq=None, pg_quantity=None):
    d = {
        "Run": run, "Precursor.Id": pid, "Protein.Group": pg,
        "Q.Value": q, "PG.Q.Value": pgq,
        "Lib.Q.Value": lib_q, "Lib.PG.Q.Value": lib_pgq, "Decoy": decoy,
    }
    if seq is not None:
        d["Stripped.Sequence"] = seq
    if gq is not None:
        d["Global.Q.Value"] = gq
    if pg_quantity is not None:
        d["PG.Quantity"] = pg_quantity
    return d


def test_per_run_proteins_pgq_only_count_contaminants():
    """Per-run protein groups filter on PG.Q.Value ONLY: contaminants are
    counted (no target filter), and a protein in every run is 'complete'."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "x1", "PGa", pgq=0.001), _row("r2", "x1", "PGa", pgq=0.001),
        _row("r1", "x4", "Cont_ALBU", pgq=0.001), _row("r2", "x4", "Cont_ALBU", pgq=0.001),
        _row("r2", "x5", "PGd", pgq=0.20),  # fails PG.Q -> excluded
    ])
    c = count_report(df)
    # r1={PGa,Cont_ALBU}=2, r2={PGa,Cont_ALBU}=2 -> avg 2 (contaminant counted)
    assert c["prot_perrun_avg"] == 2
    # complete (in both runs): PGa and Cont_ALBU
    assert c["prot_complete"] == 2


def test_per_run_proteins_count_zero_quantity():
    """A protein group with zero PG.Quantity is still counted (no
    positive-quantity filter is permitted on per-run protein numbers)."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "x1", "PGa", pgq=0.001, pg_quantity=0.0),
        _row("r2", "x1", "PGa", pgq=0.001, pg_quantity=0.0),
    ])
    c = count_report(df)
    assert c["prot_perrun_avg"] == 1
    assert c["prot_complete"] == 1


def test_per_run_precursors_qvalue_only_ignores_global_q():
    """Per-run/replicate precursor counting uses Q.Value ONLY; a bad
    Global.Q.Value must not exclude an otherwise-passing precursor."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", q=0.005, gq=0.9),
        _row("r2", "P1", "PG1", q=0.005, gq=0.9),
        _row("r3", "P1", "PG1", q=0.005, gq=0.9),
    ])
    c = count_report(df)
    assert c["prec_min1"] == 1   # counted despite Global.Q.Value=0.9
    assert c["prec_min3"] == 1


def test_per_run_precursor_q_threshold_affects_replicate_counting():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", q=0.005),
        _row("r2", "P1", "PG1", q=0.030),
        _row("r3", "P1", "PG1", q=0.030),
    ])
    assert count_report(df, precursor_q=0.01)["prec_min3"] == 0  # only r1 passes
    assert count_report(df, precursor_q=0.05)["prec_min3"] == 1  # all 3 pass


def test_global_proteins_use_lib_pg_q_only():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", lib_pgq=0.005),  # global protein counted
        _row("r1", "P2", "PG2", lib_pgq=0.05),   # fails Lib.PG.Q -> not global
    ])
    c = count_report(df)
    assert c["prot_global"] == 1   # only PG1


def test_global_precursors_use_lib_q_only():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", lib_q=0.005),
        _row("r1", "P2", "PG2", lib_q=0.05),  # fails Lib.Q -> not global
    ])
    c = count_report(df)
    assert c["prec_global"] == 1   # only P1


def test_global_counts_include_contaminants():
    """Global numbers admit only Lib.* filters -> contaminants are counted."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "P1", "PG1", lib_q=0.001, lib_pgq=0.001),
        _row("r1", "P2", "Cont_ALBU", lib_q=0.001, lib_pgq=0.001),
    ])
    c = count_report(df)
    assert c["prec_global"] == 2   # P1 + contaminant precursor
    assert c["prot_global"] == 2   # PG1 + Cont_ALBU


def test_decoys_dropped_from_every_metric():
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([_row("r1", "P1", "PG1", decoy=1)])
    c = count_report(df)
    assert c["prec_min1"] == 0
    assert c["prec_global"] == 0
    assert c["prot_global"] == 0
    assert c["prot_perrun_avg"] == 0
    assert c["prot_complete"] == 0


def test_two_peptide_proteins_global_rule():
    """prot_2pep = global protein groups (Lib.PG.Q.Value) with >= 2 distinct
    stripped peptides among Lib.Q.Value-passing precursors."""
    from analysis.count_report_ids import count_report
    df = pd.DataFrame([
        _row("r1", "p1", "PG1", seq="PEPTIDEA", lib_q=0.001, lib_pgq=0.001),
        _row("r1", "p2", "PG1", seq="PEPTIDEB", lib_q=0.001, lib_pgq=0.001),
        _row("r1", "p3", "PG2", seq="PEPTIDEC", lib_q=0.001, lib_pgq=0.001),  # 1 peptide only
    ])
    c = count_report(df)
    assert c["peptides"] == 3
    assert c["prot_2pep"] == 1   # only PG1 has >= 2 peptides
