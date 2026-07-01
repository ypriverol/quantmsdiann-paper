# Text Audit — quantmsdiann manuscript

Purpose: drive the manuscript prose to **Nature-quality** through repeated
**audit → fix → audit → fix** cycles until it converges. This file is the
standing instruction set for that work. Start a fresh session and reference it.

Target journal: **Molecular & Cellular Proteomics (MCP)**. Write the text to
**Nature quality**, not MCP-average quality. A wider, non-specialist audience
should be able to follow the narrative.

---

## Prime directive

The text is currently **VERY heavy**. The single most important job is to
**remove unreadable technical detail** and let the figures carry the numbers.
In every section ask: *is this sentence informative or useful to the reader?*
If not, cut it. Do not report every observed number; the figures speak for
themselves. The right level of detail is the one you can **justify** in each
specific case.

Do **not** invent problems that quantmsdiann supposedly solves, and do **not**
oversell the tool. Its usefulness is clear to any reader; let it stand on its
merits.

---

## The audit → fix cycle (the process to enforce)

Run this loop. Do not stop after one pass.

1. **AUDIT** — Read the full manuscript as a critical Nature editor. Produce a
   numbered list of concrete issues: heavy/unreadable passages, over-detailed
   number-dumps, undefined or misused jargon, oversold claims, invented
   problems, mis-scoped content, missing citations, and clarity failures.
   Classify each as *cut*, *rewrite*, *move*, or *add*.
2. **FIX** — Apply the fixes. Prefer cutting over rewriting; prefer moving over
   duplicating. Keep every scientific claim accurate.
3. **RE-AUDIT** — Read again from scratch. Did the fixes land? Did they create
   new roughness? List remaining issues.
4. **FIX** — Apply again.
5. **Repeat** until an audit pass yields no substantive prose issues (only
   cosmetic nits, or nothing). State explicitly when you believe it has
   converged and why.

Each cycle: report what you audited, what you changed, and what remains.

### Nature-quality bar (audit against these every pass)

- Every sentence earns its place; no filler, no hedging, no restating figures.
- Jargon is defined on first use; abbreviations are capitalised correctly and
  used consistently.
- Claims are precise and defensible; no overselling, no invented problems.
- The narrative explains **what the pipeline can do, its unique advantages and
  capabilities, and its performance**, not the implementation minutiae.
- A reader outside the immediate subfield can follow it.
- No em-dashes in prose (house style); en-dashes for ranges are fine.

---

## Specific comments to resolve

Grouped by section. Treat each as a task. Check off when fixed **and**
re-audited.

### Abstract

- [ ] Too much detail; it runs too long. Trim.
- [ ] Drop the listing of dataset types, including the special mention of DVP;
      not necessary in the abstract.
- [ ] The word **"recovered"** does not fit; replace it.

### Introduction

- [ ] Citations imbalance: **two references for timsTOF, zero for Astral**. Add
      an Astral reference (and rebalance). An Astral paper already exists in
      `references.bib` as `bubis2025astral` (Bubis et al. 2025, *Nat. Methods*);
      consider citing it or a primary Astral instrument reference.
- [ ] Replace "limits its scalability for large datasets" with
      **"limits its scalability for en masse processing of extra large
      datasets."**
- [ ] "Deployments on HPC remain difficult" is **not true**, and it oversells
      quantmsdiann. Remove/rework. Usefulness is clear without conjuring this.
- [ ] "Releases difficult to reconcile" — **not a real problem.** Anyone can
      reuse pipelines that establish exact settings; exact reanalysis *is*
      possible. Do not manufacture a problem the tool supposedly solves. Remove.

### Methods / plexDIA

- [ ] "MS1 and MS2 scan mass-window intervals" — unclear what this means.
      Clarify or remove.
- [ ] "heterogeneous plexDIA channel layouts cannot be mixed in same
      experiment" is simply **impossible**; unclear what the output would
      even be. Reframe: *any* setting heterogeneity within an experiment
      invalidates relative quantification, regardless of the software's best
      efforts. State that principle instead.
- [ ] **"psm" → "PSM"**: capitalise as an abbreviation, **define on first use**,
      and **audit the whole text for similar cases** (other abbreviations not
      capitalised/defined).
- [ ] "or a plexDIA dataset run with its channel definitions" is **not
      understandable.** Rewrite.

### Results (general)

- [ ] Too much technical detail. Instead, explain **what the pipeline can do,
      its unique advantages and capabilities, and its performance.**
- [ ] **Agent task:** examine the narrative and ensure the Results contain only
      information that is informative or useful to the reader. Nature-quality
      writing for a wider audience.

### Results — Section 3.1 ("architecture")

- [ ] **How we benchmark does not belong in 3.1.** A section declared as
      "architecture" should not carry the benchmarking methodology. Move it out
      (see the benchmarking-methodology section proposed below).
- [ ] "incomplete execution trace": the reader cannot tell what the actual
      issue was. Clarify. **Agent task:** audit the text for similar
      unexplained/obscure phrasings and fix them.

### Results — Section 3.2 (versions)

- [ ] "rises monotonically with DIA-NN version" — with only **two versions**
      tested, "monotonically" is unwarranted. Change **"rises" → "increases"**
      and drop "monotonically".

### Results — Section 3.3

- [ ] **Too detailed.** No need to report every single observed number. Figures
      speak for themselves. Apply the Nature-quality bar: the exact appropriate
      level of detail, justified in each case.

### Results — Section 3.4

- [ ] **Move FDR control and limitations out of 3.4** to *before* the numbers
      comparisons begin (see below).
- [ ] The **disclaimer itself in 3.4 is excellent**; keep it, but **specify
      which releases are meant** so the reader does not assume "everything before
      2.5". Reword along the lines of: *"Whereas older releases, such as the
      1.7-series and 1.8.1, and other search engines ..."*
- [ ] Add an extra reference to that disclaimer: Fröhlich et al. 2022,
      *Nat. Commun.* 13:2622 (DOI 10.1038/s41467-022-30094-0), which benchmarks
      many DIA analysis strategies / search engines. Already added to
      `references.bib` as `frohlich2022benchmarking` — cite this key.

### Results — level of detail on datasets

- [ ] **Do not reference every single PXD in the text** if it is already given
      in the figure legend. Instead: state the **overall conclusions**, support
      them with a reference to the figure, and give a **bulk citation** (paper
      numbers separated by commas) to the papers that generated the data.
- [ ] **DVP / 1.8.1 caveat (currently not explained):** 1.8.1 for DVP was run on
      **extra runs**, which likely increased the globally identifiable proteins
      for 1.8.1. Note in the text that **even despite this, 2.5.1 achieved higher
      numbers.**

---

## New content required (critical, for scientific correctness)

- [ ] **Add a paragraph on FDR control before any identification numbers are
      described.** It must cover:
  - loose/relaxed FDR reporting in **previous DIA-NN versions**, and
  - the **settings mismatch** between quantmsdiann and the public analyses being
    compared against.
  This is **critical** and required for scientific correctness — the reader must
  understand these caveats *before* seeing the ID counts.

### Proposed restructure: a Benchmarking Methodology section

- [ ] Create a dedicated **benchmarking methodology** section that combines:
  - the **list of benchmarks / how we benchmark** cut out of 3.1, and
  - the **FDR and comparison disclaimers** moved out of 3.4 (including the new
    FDR paragraph above).
  Place it before the identification-number comparisons.
  **Refine the document first, then execute this move.**

---

## Declarations

- [ ] **Add an AI-use disclaimer** stating that AI was used to prepare the
      manuscript. Place it in the appropriate declarations section (per MCP
      author guidelines) and keep the wording factual and concise.

---

## Discussion

- [ ] **The Discussion must not repeat Intro and Results.** Only add content
      that brings extra value. Nature-quality writing throughout.
- [ ] It is acceptable to **slim it down to a single paragraph** summarising the
      impact of the work and the outlook, if that is what serves the reader best.

---

## Figures

### Placement

- [ ] **Figures currently break the Discussion in two; fix this.** Put figures
      *before* the Discussion. Ideally place **each figure immediately after the
      Results paragraph that first mentions it.**

### Figure 4a (DVP)

- [ ] Add **"18 runs"** and **"12 runs"** directly on the figure, to the right of
      the global confident-protein numbers (trivial edit). This ties to the
      DVP / 1.8.1 extra-runs caveat noted above.

### Supplementary figures

- [ ] **Figure S1:** not readable. Add the **median numbers next to each box.**
- [ ] **Figure S1:** programmatic variable names in **CAPS are not appropriate**
      in the paper; use readable labels.
- [ ] **Figure S2:** same CAPS/variable-name issue. Also **too tall**; compute
      appropriate look and dimensions.
- [ ] **Figure S3:** label font is **tiny.** Fix.

### Whole-paper figure audit

- [ ] **Audit element proportionality across the whole paper** (fonts, label
      sizes, aspect ratios).
- [ ] **Audit every figure referenced in the main text, including all
      supplementary figures**, for readability and consistency.

---

## Working order

1. Refine/agree this plan and the restructure before large edits.
2. Add the FDR-control paragraph and stand up the benchmarking-methodology
   section (moves + new content).
3. Resolve the section-specific comments above.
4. Slim the Discussion (no repetition of Intro/Results; extra value only).
5. Fix figure placement, Figure 4a labels, and the supplementary figures; run
   the whole-paper figure/proportionality audit.
6. Run the audit → fix cycle repeatedly until convergence.
7. Report convergence with justification.

## Definition of done

- All checkboxes above resolved.
- An audit pass returns no substantive prose issues.
- Text reads at Nature quality for a wider audience; heavy technical detail
  removed; figures carry the numbers.
- Every claim accurate; no oversold statements; no invented problems.
- Once the above holds: **green light to submit the first version.**
