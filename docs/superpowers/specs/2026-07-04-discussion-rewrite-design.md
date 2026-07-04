# Discussion rewrite — design

**Date:** 2026-07-04
**File touched:** `paper/main.tex` (Discussion section, currently 4 paragraphs), `paper/references.bib`
**Goal:** Rework the Discussion so it interprets the results and projects their significance, instead of re-running the Introduction. Fold in themes the corresponding author asked for (reanalysis-as-practice, DIA-is-standard, metadata/SDRF, scale needing fast standardized solutions, method breadth, ease/provenance, knowledge-base outlook) and integrate two supplied citations.

## Key decisions (from brainstorming)
1. **Reframe forward, don't restate.** The Discussion assumes the Intro's setup and pivots to interpretation + significance. "Why not quantms" (mechanics: fixed 1.8.1, DDA/DIA conflicts, hard-coded versions, mandatory conversion) stays in the Introduction and is NOT repeated.
2. **Method breadth = inherit-from-DIA-NN framing.** Claim only demonstrated modalities as shown (single-cell, spatial DVP, plexDIA, phosphoproteomics); name immunopeptidomics as a natural same-SDRF-path extension, not a result. No immunopeptidomics overclaim.
3. Style constraints: MCP prose, **no em-dashes** (use colons/en-dashes for ranges), ~500–540 words, all numbers must match the manuscript's established figures.

## Two new references to add to `references.bib`
- **Hewapathirana et al. 2026** (openRxiv/bioRxiv), "Quantifying data reuse in proteomics using PRIDE download statistics and a semi-supervised LLM-based framework", DOI `10.64898/2026.04.16.718670`. Key: `hewapathirana2026datareuse`.
- **Perez-Riverol 2022**, "Proteomic repository data submission, dissemination, and reuse: key messages", *Expert Rev Proteomics* 19(7-12):297-310, DOI `10.1080/14789450.2022.2160324`, PMID 36529941. Key: `perezriverol2022reuse`.

## Paragraph arc (5 paragraphs)

**¶1 — Reanalysis has become how proteomics builds knowledge (hook + payoff).**
Public-data reuse is now part of the scientific process, not only to complement local data but to build more robust cross-cohort biological models `[hewapathirana2026datareuse, perezriverol2022reuse]`. DIA is now a dominant, fastest-growing share of deposits, so reusable raw data accumulates faster than it is re-mined. Land the two gains as payoff: distributed run reproduces desktop DIA-NN (Fig 2c); older/non-DIA-NN deposits recover up to +59% protein groups; gain concentrated in the ageing back-catalogue.

**¶2 — Metadata is the key that turns reanalysis into biology (SDRF).**
More identifications are necessary but not sufficient: extracting biology from pooled public data requires knowing what each sample is `[perezriverol2022reuse]`. Standards like SDRF `[dai2021sdrf]` make design machine-readable; quantmsdiann treats SDRF as the single source of truth, so one file both drives the per-run DIA-NN analysis and annotates the output. Callback: pan-cohort matrix regroupable by Cellosaurus/NCIt/Disease Ontology with no manual harmonisation.

**¶3 — Scale demands fast, standardized, provenance-complete analysis (design payoff).**
As cohorts grow, the bottleneck moves from instrument to analysis: reprocess thousands of runs quickly AND emit data + metadata in a standardized, queryable form. quantmsdiann does exactly that: single SDRF + container-pinned profiles + no manual install (ease); 2,300 runs in 2.2 h, compute-bounded (scale); all parameters captured in SDRF/`qpx` (provenance). Light nod to the DIA-NN-dedicated independent design (no quantms-limitation list).

**¶4 — Method breadth + cross-version reproducibility.**
Wrapping DIA-NN unmodified and parametrising per-run from SDRF, quantmsdiann inherits DIA-NN's method range: demonstrated on single-cell, spatial DVP, plexDIA, phosphoproteomics; other DIA-NN workflows such as immunopeptidomics `[scheid2025mhcquant2]` follow the same SDRF path unchanged. Keep the cross-version DIA-NN reproducibility point (novel, not in Intro) and the downstream handoff to scp `[vanderaa2023]`.

**¶5 — Outlook: a proteomics knowledge base.**
Turning any SDRF-annotated deposit into harmonised, queryable quantifications positions this to power DIA reanalysis at archive scale and integrate independent cohorts toward an aggregated knowledge base `[lautenbacher2022proteomicsdb, perezriverol2025pride]`.

## Deliberately excluded (already in the Introduction)
- The quantms limitation list and "moving away from quantms" mechanics.
- The "most DIA-NN users run it as a desktop app" problem.
- The heterogeneous-hardware / distributed-computing setup argument.

## Verification
- Rebuild `main.pdf`; confirm no undefined refs/citations and the two new keys resolve.
- Grep the built PDF/`.tex` for em-dashes in the new prose.
- Cross-check every number against the manuscript (59%, 2,300 runs / 2.2 h / 300 nodes).
