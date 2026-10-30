# Integrating the cell-line cohorts: protein FDR across studies (reviewer response R1.9)

Each cohort in the manuscript was analysed in its own quantmsdiann run, so no
identification is transferred between studies. Within a run, FDR is controlled over
the whole cohort: precursor and protein-group q-values are estimated once on the
cohort's empirical library (`Lib.Q.Value` and `Lib.PG.Q.Value` ≤ 0.01).

All five cell-line cohorts were processed the same way: quantmsdiann de98e98
(v2.2.0dev), DIA-NN 2.5.1, `Homo-sapiens-uniprot-reviewed-contaminants-202605.fasta`
(no entrapment) and the same search parameters (`sdrf/runs.tsv`).

To show what happens when independent cohorts are combined, we integrated the five
cohorts with the strategy of the original quantms study (Dai et al. 2024,
*Nat Methods*):

1. 1% protein FDR per cohort (`Lib.PG.Q.Value` ≤ 0.01);
2. at least two unique peptides per protein (distinct proteotypic sequences at `Lib.Q.Value` ≤ 0.01);
3. protein q-values re-estimated on the merged list with decoys, keeping proteins at adjusted q < 0.01.

**Result.** Of 11,654 proteins in the merged list, 11,546 passed 1% adjusted FDR (108
removed, 0.9%). This is comparable to the cell-line build of the quantms study
(11,403 to 11,374). Without the two-peptide rule, re-estimation removes 476 of 12,953
proteins (3.7%), so the two-peptide rule does most of the work of controlling
protein FDR across studies (Response Table 2).

## Response Table 2

Proteins counted per accession (every member of a protein group); contaminants
removed.

| Set | 1% protein FDR | + ≥2 unique peptides |
|---|---:|---:|
| PXD003539 (NCI-60) | 7,379 | 6,577 |
| PXD004701 (Sun 2023) | 8,224 | 7,138 |
| PXD030304 (ProCan) | 9,918 | 8,688 |
| PXD017199 | 10,860 | 9,993 |
| PXD041421 | 9,370 | 8,611 |
| Merged (union of the five lists) | 12,953 | 11,654 |
| Merged, decoy-adjusted q < 0.01 | 12,477 | **11,546** |

Proteins removed from each cohort are listed in
[`integration_removed_by_project_min2_unique_peptides.tsv`](integration_removed_by_project_min2_unique_peptides.tsv)
and [`integration_removed_by_project_1pct_protein_fdr.tsv`](integration_removed_by_project_1pct_protein_fdr.tsv).

![Cell-line cohorts integrated with the quantms criteria](response_R1.9_integration.svg)

**(a)** Proteins per cohort and after merging (≥2 unique peptides, 1% protein FDR).
**(b)** Proteins removed by re-estimating q-values on the merged list, for the
quantms 2024 cell-line build and for this study.

## Method

- **Per cohort.** Protein groups at `Lib.PG.Q.Value` ≤ 0.01 from each cohort's DIA-NN
  report, with the number of unique peptides (distinct proteotypic `Stripped.Sequence`
  at `Lib.Q.Value` ≤ 0.01). Each accession keeps its lowest q-value.
- **Decoy re-estimation (quantms).** In each cohort list, sorted by q-value, decoys
  are added so that the decoy/target ratio at every target equals its q-value. The
  lists are merged, each accession keeps its lowest q-value, and q-values are
  re-estimated as decoys / targets scoring as well or better (made monotone).
  Proteins at adjusted q < 0.01 are kept.

The same cohort lists define the pan-cohort union of the atlas (Fig. 4,
[`scripts/cellline_accessions.py`](../../../../scripts/cellline_accessions.py)).

## Reproduce

```bash
python -m scripts.review_integration_fdr      # tables in this folder
python -m scripts.review_integration_figure   # response_R1.9_integration.svg
```

Code: [`scripts/review_integration_fdr.py`](../../../../scripts/review_integration_fdr.py)
and [`scripts/review_integration_figure.py`](../../../../scripts/review_integration_figure.py).
Reports are read through [`scripts/cellline_reports.py`](../../../../scripts/cellline_reports.py):
a local copy in `data/cell-lines-v2_5_1/<PXD>/diann_report.parquet` when present,
otherwise `quantmsdiann-benchmarks/cell-lines/<PXD>/v2_5_1/quant_tables/diann_report.parquet`
on the PRIDE FTP.

## Files

| File | Content |
|---|---|
| `integration_fdr.tsv` | Response Table 2, both rules |
| `integration_removed_by_project_*.tsv` | Proteins of each cohort removed by the decoy step |
| `response_R1.9_integration.svg` | Figure above |

This analysis is part of the response to the reviewers only; the resource does not
claim integrated protein lists across cohorts.
