#!/usr/bin/env python3
"""rebuild — reproduce every manuscript figure (and the PDFs) from original
sources, in one command.

All data is pulled from the public PRIDE FTP (deposited DIA-NN reports and the
original published matrices); nothing here re-runs DIA-NN from raw vendor files
(that is the upstream HPC reanalysis, out of scope). Every figure number obeys
the filter rule documented in ``methods.md`` §1.

Usage:
    python -m scripts.rebuild --all            # data prep -> figures -> PDFs
    python -m scripts.rebuild --figures-only    # skip data prep
    python -m scripts.rebuild --data-only       # only the data-prep stages
    python -m scripts.rebuild --only NAME [...]  # run specific stage(s)
    python -m scripts.rebuild --no-pdf          # skip the LaTeX/PDF build
    python -m scripts.rebuild --list            # print all stages + provenance
    python -m scripts.rebuild --keep-going      # don't stop on first failure (default)
    python -m scripts.rebuild --fail-fast       # stop at the first failing stage

Each stage is a Python module run as ``python -m <module>`` so it stays
identical to running the script by hand. The atlas (Fig S13) needs a numpy<2
environment (see methods.md / environment.yml); run rebuild inside the
``quantmsdiann`` conda env to reproduce it.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Stage = (name, module, description). Order matters: data prep before figures.
DATA_PREP: list[tuple[str, str, str]] = [
    ("report_counts", "analysis.recompute_report_counts",
     "Download FTP reports, recount under the filter rule -> report_counts.tsv"),
    ("single_cell_tables", "analysis.make_single_cell_tables",
     "Per-cell / completeness / CV tables for the single-cell figure"),
    ("phospho_tables", "analysis.make_phospho_tables",
     "Phosphopeptide / phosphosite tables (Lib.Q.Value; site Prob>=0.99)"),
]

FIGURES: list[tuple[str, str, str]] = [
    ("benchmarks", "analysis.figure_quantmsdiann_benchmarks_vs_proteobench",
     "Fig 2 — ProteoBench benchmark panels + counts.tsv"),
    ("fig2_validation", "analysis.figure_fig2_validation",
     "Fig 2 validation composite"),
    ("id_vs_epsilon", "analysis.figure_id_vs_epsilon",
     "ProteoBench id-vs-epsilon panel"),
    ("proteobench_accuracy", "analysis.figure_proteobench_accuracy",
     "ProteoBench accuracy panels"),
    ("reanalysis_improvement", "analysis.figure_reanalysis_improvement",
     "Reanalysis-recovery figure (original vs quantmsdiann)"),
    ("single_cell_combined", "analysis.figure_single_cell_combined",
     "Single-cell figure (per-cell depth, completeness, CV)"),
    ("plexdia_per_cell", "analysis.plexDIA.figure_msv000093870_oocyte_plexdia",
     "plexDIA per-cell depth (MSV000093870)"),
    ("plexdia_vs_galatidou", "analysis.plexDIA.figure_msv000093870_galatidou_vs_quantmsdiann",
     "plexDIA deposited vs quantmsdiann"),
    ("pxd003539", "analysis.figure_original_vs_quantmsdiann",
     "PXD003539 (NCI-60) panels"),
    ("pxd004701", "analysis.figure_pxd004701_sun_vs_quantmsdiann",
     "PXD004701 (Sun) panels"),
    ("pxd030304", "analysis.figure_pxd030304_procan_vs_quantmsdiann",
     "PXD030304 (ProCan) panels — streams a 2GB matrix"),
    ("pxd064049_spatial", "analysis.figure_pxd064049_spatial_vs_quantmsdiann",
     "PXD064049 spatial DVP panels"),
    ("atlas", "analysis.figure_combined_cell_lines_atlas",
     "Pan-cohort of DIA reanalyses (Fig S13) — needs numpy<2"),
    ("phospho", "analysis.figure_phospho",
     "Phosphoproteomics supplementary figure"),
    ("venn", "analysis.venn_protein_accessions",
     "Protein-accession overlap (supplementary)"),
    ("queue_sweep", "analysis.figure_queue_size_sweep",
     "Fig 1 queue-size sweep"),
    ("performance_trace", "analysis.figure_performance_trace",
     "Per-step runtime + resources (runtime_per_step.svg, resources_per_step.svg)"),
    ("mdc_cluster_runtime", "analysis.figure_mdc_cluster_runtime",
     "MDC cluster runtime"),
]

# Runs LAST: aggregates every cited number from the figure-data TSVs into one
# file (data/paper_numbers.tsv + paper/generated_numbers.tex).
COLLECT: list[tuple[str, str, str]] = [
    ("paper_numbers", "analysis.collect_paper_numbers",
     "Aggregate ALL manuscript numbers -> data/paper_numbers.tsv + paper/generated_numbers.tex"),
]

ALL_STAGES = DATA_PREP + FIGURES + COLLECT
BY_NAME = {name: (name, mod, desc) for name, mod, desc in ALL_STAGES}


def print_list() -> None:
    print("DATA PREP:")
    for name, mod, desc in DATA_PREP:
        print(f"  {name:22s} {desc}")
    print("\nFIGURES:")
    for name, mod, desc in FIGURES:
        print(f"  {name:22s} {desc}")
    print("\nNUMBERS:")
    for name, mod, desc in COLLECT:
        print(f"  {name:22s} {desc}")
    print("\nPDFs: paper/Makefile -> figures, pdf, supplementary")


def run_stage(name: str, module: str) -> tuple[str, bool, float]:
    print(f"\n=== [{name}] python -m {module} ===", flush=True)
    t0 = time.time()
    rc = subprocess.run([sys.executable, "-m", module], cwd=REPO).returncode
    dt = time.time() - t0
    ok = rc == 0
    print(f"=== [{name}] {'OK' if ok else 'FAILED (rc=%d)' % rc} in {dt:.0f}s ===",
          flush=True)
    return name, ok, dt


def run_make(target: str) -> tuple[str, bool, float]:
    print(f"\n=== [make {target}] ===", flush=True)
    t0 = time.time()
    rc = subprocess.run(["make", target], cwd=REPO / "paper").returncode
    dt = time.time() - t0
    return f"make {target}", rc == 0, dt


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="data prep + figures + PDFs")
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--data-only", action="store_true")
    ap.add_argument("--only", nargs="+", metavar="NAME", help="run specific stage(s)")
    ap.add_argument("--no-pdf", action="store_true", help="skip the PDF build")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    ap.add_argument("--fail-fast", action="store_true", help="stop at first failure")
    ap.add_argument("--keep-going", action="store_true",
                    help="continue past failures (default)")
    args = ap.parse_args(argv)

    if args.list:
        print_list()
        return 0

    if args.only:
        unknown = [n for n in args.only if n not in BY_NAME]
        if unknown:
            print(f"unknown stage(s): {', '.join(unknown)}", file=sys.stderr)
            print("use --list to see valid names", file=sys.stderr)
            return 2
        stages = [BY_NAME[n] for n in args.only]
        build_pdf = False
    elif args.data_only:
        stages, build_pdf = DATA_PREP, False
    elif args.figures_only:
        stages, build_pdf = FIGURES + COLLECT, not args.no_pdf
    elif args.all:
        stages, build_pdf = ALL_STAGES, not args.no_pdf
    else:
        ap.print_help()
        return 0

    fail_fast = args.fail_fast and not args.keep_going
    results: list[tuple[str, bool, float]] = []
    for name, mod, _ in stages:
        r = run_stage(name, mod)
        results.append(r)
        if not r[1] and fail_fast:
            break

    if build_pdf and (not fail_fast or all(ok for _, ok, _ in results)):
        for target in ("figures", "pdf", "supplementary"):
            results.append(run_make(target))

    print("\n================ SUMMARY ================")
    failed = [n for n, ok, _ in results if not ok]
    for name, ok, dt in results:
        print(f"  {'OK  ' if ok else 'FAIL'}  {name:30s} {dt:6.0f}s")
    if failed:
        print(f"\n{len(failed)} stage(s) failed: {', '.join(failed)}")
        return 1
    print("\nAll stages succeeded.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
