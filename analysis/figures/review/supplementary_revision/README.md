# New supplementary items for the revision

The reviewer replies and the edit list live locally in `paper/revision/response_drafts.md` (not pushed).

| File | Goes in | Answers |
|---|---|---|
| FigS6_standalone_vs_quantmsdiann (pdf, png) | Supplementary Note 8; replaces current Figs. S6 and S7 | R1.1, R1.7, R1.18 |
| TableS3_standalone_vs_quantmsdiann.tsv | Supplementary Note 8 | R1.1, R1.7, R1.18 |
| FigS7_runtime_vs_files (pdf, png) | Supplementary Note 4 | R1.24 |
| TableS4_runtime_vs_scale.tsv | Supplementary Note 4 | R1.24 |
| FigS8_proteins_by_peptide_evidence (pdf, png) | Supplementary Note 6 | R1.10, R1.21 |
| TableS5_proteins_by_peptide_evidence.tsv | Supplementary Note 6 | R1.10, R1.11, R1.12, R1.21 |
| TableS6_run_parameters.tsv | Supplementary Note 2 | R1.3, R1.6 |

Filtering (all items): global counts at Lib.Q.Value ≤ 0.01 (precursors) and
Lib.PG.Q.Value ≤ 0.01 (protein groups); per-run counts at Q.Value ≤ 0.01 and
PG.Q.Value ≤ 0.01; decoys dropped; no other filter (Supplementary Note 3).

Rebuild: `python -m scripts.review_supplementary` (after the three
`scripts/review_*` analysis scripts) regenerates and copies everything here.



Notes for whoever updates the manuscript:
- Every number in the drafts was checked against these tables. The evidence
  counts (Table S5) reproduce Fig. 4 exactly (NCI-60 6,812; ProCan 8,993 and
  10,027; Sun 8,304 and 7,381).
- Items marked [you] need facts from the authors (cluster specs, Zenodo DOI,
  quantms differences, SDRF defaults, comparison table).
- "Decide before sending" points: R1.12 (45% peptide recovery) and the
  Fig. 4 / Supplementary Note 3 wording for the NCI-60 4,284 (R1.16).
