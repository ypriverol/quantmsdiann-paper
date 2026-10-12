"""Turn the crawled file manifest into the PRIDE submission and the Zenodo file list.

Reads pride_submission/submission_files.tsv (from scripts.build_submission_manifest)
and writes, under pride_submission/:

  sdrf/                      SDRF copies to upload: byte-identical to the as-run SDRF
                             snapshot in this repo (sdrf/)
  submission_files_final.tsv every file to upload (FTP URL or local path), with
                             its PRIDE type and whether it also goes to Zenodo
  submission.px              PRIDE submission summary (metadata + file mapping)
  zenodo_files.tsv           the Zenodo subset

Local additions: the paired standalone DIA-NN vs quantmsdiann runs (PXD049412,
PXD070049; Supplementary Note 8), from data/proteobench-review-experiments/x/.

Usage:
    python -m scripts.build_submission_px
"""

from __future__ import annotations

import io
import re
import shutil
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "pride_submission"
REVIEW = REPO / "data/proteobench-review-experiments/x"
META = REPO / "run_metadata"

# Wrong PSI-MS accessions found in as-run SDRFs (checked against OLS, 2026-10-11). Applied only to the
# instrument list in submission.px; the uploaded SDRFs stay as run.
INSTRUMENT_FIXES = {
    "AC=MS:1003356": "AC=MS:1003378",              # labelled Orbitrap Astral; 1003356 is Orbitrap Ascend
    "AC=MS:1003387": "AC=MS:1003443",              # labelled ZenoTOF 8600; 1003387 is MetaboAnnotation
    "AC=MS:1003005;NT=timsTOF SCP": "AC=MS:1003231;NT=timsTOF SCP",
    "NT=Q Exactive Plus;AC=MS:1001911": "NT=Q Exactive Plus;AC=MS:1002634",
}
# Zenodo gets everything except the single 44 GB ProCan report (in PRIDE only).
ZENODO_EXCLUDE = {"PXD030304_v2_5_1_diann_report.parquet"}

SUBMITTER = {
    "submitter_name": "Yasset Perez-Riverol", "submitter_email": "yperez@ebi.ac.uk",
    "submitter_affiliation": "European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI)",
    "submitter_pride_login": "yperez@ebi.ac.uk",
    "lab_head_name": "Yasset Perez-Riverol", "lab_head_email": "yperez@ebi.ac.uk",
    "lab_head_affiliation": "European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI)",
}
TITLE = "quantmsdiann reanalysis of public single-cell, spatial, phosphoproteomics and bulk DIA datasets"
DESCRIPTION = (
    "Reanalysis of public data-independent acquisition (DIA) datasets with quantmsdiann "
    "(https://github.com/bigbio/quantmsdiann, v2.2.0), an SDRF-driven Nextflow workflow that runs DIA-NN in "
    "parallel. Includes ProteoBench DIA benchmarks (modules 5, 7, 9, 10) with DIA-NN 1.8.1, 2.5.1 and "
    "2.5.1-enterprise; single-cell (PXD046357, PXD049412, MSV000093870), spatial (PXD064049), "
    "phosphoproteomics (PXD034128, PXD034623, PXD049692) and bulk cell-line (PXD003539, PXD004701, PXD030304, "
    "PXD017199, PXD041421) reanalyses; the PXD071075 scalability run; and paired standalone DIA-NN versus "
    "quantmsdiann runs (PXD049412, PXD070049). Raw files are not re-uploaded: raw_files.tsv lists every raw "
    "file processed with its URL in the original deposit. Analysis code: https://github.com/ypriverol/quantmsdiann-paper."
)
SAMPLE_PROTOCOL = "See the original submissions listed as reanalysed datasets; no new samples were prepared."
DATA_PROTOCOL = (
    "Each dataset was analysed in a separate quantmsdiann run (v2.2.0) from its SDRF, library-free (DIA-NN "
    "in-silico library predicted from the FASTA), with DIA-NN 1.8.1, 2.5.0, 2.5.1 or 2.5.1-enterprise as "
    "indicated in each file name (<dataset>_<version>_<file>). Per-run parameters are in the diann_config.cfg, "
    "params_*.json and diannsummary.log files. Identification counts in the manuscript use Lib.Q.Value <= 0.01 "
    "(precursors) and Lib.PG.Q.Value <= 0.01 (protein groups) for global counts, Q.Value / PG.Q.Value <= 0.01 per run."
)
KEYWORDS = "DIA, reanalysis, DIA-NN, quantmsdiann, single-cell proteomics, SDRF, Nextflow"
REANALYSIS = ["PXD003539", "PXD004701", "PXD017199", "PXD030304", "PXD034128", "PXD034623", "PXD041421",
              "PXD046357", "PXD049412", "PXD049692", "PXD062685", "PXD064049", "PXD070049", "PXD071075"]
SPECIES = [("9606", "Homo sapiens (human)"), ("4932", "Saccharomyces cerevisiae (baker's yeast)"),
           ("562", "Escherichia coli")]
MODS = [("MOD:01060", "S-carboxamidomethyl-L-cysteine"), ("MOD:00719", "L-methionine sulfoxide"),
        ("MOD:00394", "acetylated residue"), ("MOD:00696", "phosphorylated residue")]


def instruments(sdrf_dir: Path) -> dict[str, str]:
    out = {}
    for f in sdrf_dir.glob("*.tsv"):
        text = f.read_text()
        for k, v in INSTRUMENT_FIXES.items():
            text = text.replace(k, v)
        s = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
        for c in [c for c in s.columns if "instrument" in c.lower()]:
            for v in s[c].dropna().unique():
                parts = dict(p.split("=", 1) for p in v.split(";") if "=" in p)
                if "AC" in parts and "NT" in parts:
                    out[parts["AC"]] = parts["NT"]
    return out


def review_files() -> list[dict]:
    rows = []
    for ds in ("PXD049412", "PXD070049"):
        sa = REVIEW / ds / "diann"
        qm = REVIEW / ds / f"quantmsdiann/results_{ds}"
        picks = [("standalone", "SEARCH", sa / "diann_report.parquet")]
        picks += [("standalone", "RESULT", sa / f) for f in ("diann_report.pg_matrix.tsv", "diann_report.pr_matrix.tsv",
                                                            "diann_report.stats.tsv")]
        picks += [("standalone", "OTHER", sa / f) for f in ("diann_report.log.txt", "run_diann_v2_5_1.sh",
                                                           "diann_native_execution.log")]
        picks += [("quantmsdiann", "SEARCH", qm / "quant_tables/diann_report.parquet")]
        picks += [("quantmsdiann", "RESULT", qm / f"quant_tables/{f}") for f in
                  ("diann_report.pg_matrix.tsv", "diann_report.pr_matrix.tsv", "diann_report.stats.tsv")]
        picks += [("quantmsdiann", "OTHER", qm / "quant_tables/diannsummary.log"),
                  ("quantmsdiann", "OTHER", REVIEW / ds / "quantmsdiann/run_quantmsdiann.config")]
        picks += [("quantmsdiann", "OTHER", p) for p in sorted((qm / "pipeline_info").glob("params_*.json"))[-1:]]
        picks += [("quantmsdiann", "EXPERIMENTAL_DESIGN", p) for p in sorted((qm / "sdrf").glob("*.sdrf.tsv"))]
        for wf, ftype, p in picks:
            if not p.exists():
                print(f"WARNING missing local file {p}")
                continue
            rows.append({"dataset": f"Paired-{wf}-{ds}", "diann_version": "v2_5_1", "file_type": ftype,
                         "upload_name": f"Paired-{wf}-{ds}_v2_5_1_{p.name}", "source_url": str(p),
                         "size_bytes": p.stat().st_size})
    return rows


def main() -> None:
    d = pd.read_csv(OUT / "submission_files.tsv", sep="\t")
    sdrf_dir = OUT / "sdrf"
    shutil.rmtree(sdrf_dir, ignore_errors=True)
    sdrf_dir.mkdir()
    runs = pd.read_csv(META / "runs.tsv", sep="\t", dtype=str)
    # PXD071075: the paper uses the v2.5.0 node-count sweep (Fig. 2a), not the 2.5.1-enterprise FTP run.
    rows = [r for r in d.to_dict("records") if r["file_type"] != "EXPERIMENTAL_DESIGN" and r["dataset"] != "PXD071075"]
    for _, run in runs.iterrows():
        ds, ver = run["dataset"], run["run"]
        if ds.startswith("Paired-"):
            continue  # added by review_files()
        tag = f"{ds}_{ver}"
        dest = sdrf_dir / f"{tag}_{Path(run['sdrf']).name}"
        shutil.copyfile(REPO / run["sdrf"], dest)
        rows.append({"dataset": ds, "diann_version": ver, "file_type": "EXPERIMENTAL_DESIGN", "upload_name": dest.name,
                     "source_url": str(dest), "size_bytes": dest.stat().st_size,
                     "note": "as-run SDRF"})
        mdir = META / ds / ver
        on_ftp = str(run.get("params_source", "")).startswith("FTP")
        for f in ("params.json", "launch.txt") if not on_ftp else ():
            if (mdir / f).exists():
                rows.append({"dataset": ds, "diann_version": ver, "file_type": "OTHER", "upload_name": f"{tag}_{f}",
                             "source_url": str(mdir / f), "size_bytes": (mdir / f).stat().st_size,
                             "note": "run parameters recovered from the cluster run log"})
        if ds == "PXD071075":
            q = re.search(r"sweep_(\d+)cores", ver).group(1)
            trace = REPO / f"data/queue_size_sweep/q{q}/nextflow_trace.txt"
            if trace.exists():
                rows.append({"dataset": ds, "diann_version": ver, "file_type": "OTHER",
                             "upload_name": f"{tag}_nextflow_trace.txt", "source_url": str(trace),
                             "size_bytes": trace.stat().st_size, "note": "Fig. 2a node-count sweep trace"})
    for f in (META / "runs.tsv", META / "common/extra.config",
              REPO / "analysis/figures/review/supplementary_revision/TableS6_run_parameters.tsv"):
        rows.append({"dataset": "all", "diann_version": "", "file_type": "OTHER", "upload_name": f.name,
                     "source_url": str(f), "size_bytes": f.stat().st_size, "note": "run parameters for all runs"})
    rows += review_files()
    for name in ("raw_files.tsv", "raw_files.txt"):
        p = OUT / name
        rows.append({"dataset": "all", "diann_version": "", "file_type": "OTHER", "upload_name": name,
                     "source_url": str(p), "size_bytes": p.stat().st_size, "note": "raw files used, with original URLs"})
    for r in rows:
        r.setdefault("note", "")
        r["zenodo"] = "no" if r["upload_name"] in ZENODO_EXCLUDE else "yes"
    df = pd.DataFrame(rows)
    assert not df.upload_name.duplicated().any(), df[df.upload_name.duplicated()]
    df.to_csv(OUT / "submission_files_final.tsv", sep="\t", index=False)
    df[df.zenodo == "yes"].to_csv(OUT / "zenodo_files.tsv", sep="\t", index=False)

    inst = instruments(sdrf_dir)
    lines = [f"MTD\t{k}\t{v}" for k, v in SUBMITTER.items()]
    lines += [f"MTD\tproject_title\t{TITLE}", f"MTD\tproject_description\t{DESCRIPTION}",
              f"MTD\tsample_processing_protocol\t{SAMPLE_PROTOCOL}", f"MTD\tdata_processing_protocol\t{DATA_PROTOCOL}",
              f"MTD\tkeywords\t{KEYWORDS}", "MTD\tsubmission_type\tPARTIAL",
              "MTD\texperiment_type\t[PRIDE, PRIDE:0000450, Data-independent acquisition, ]"]
    lines += [f"MTD\tspecies\t[NEWT, {t}, {n}, ]" for t, n in SPECIES]
    lines += ["MTD\ttissue\t[PRIDE, PRIDE:0000442, Tissue not applicable to dataset, ]"]
    lines += [f"MTD\tinstrument\t[MS, {a}, {n}, ]" for a, n in sorted(inst.items())]
    lines += [f"MTD\tmodification\t[MOD, {a}, {n}, ]" for a, n in MODS]
    lines += [f"MTD\treanalysis\t{a}" for a in REANALYSIS]
    lines += ["FMH\tfile_id\tfile_type\tfile_path\tfile_mapping"]
    lines += [f"FME\t{i}\t{r['file_type']}\t{r['upload_name']}\t" for i, r in enumerate(rows, 1)]
    (OUT / "submission.px").write_text("\n".join(lines) + "\n")

    tot = df.size_bytes.sum() / 1e9
    z = df.loc[df.zenodo == "yes", "size_bytes"].sum() / 1e9
    print(f"PRIDE: {len(df)} files, {tot:.1f} GB | Zenodo: {(df.zenodo == 'yes').sum()} files, {z:.1f} GB")
    print("instruments:", inst)
    print(df.note.value_counts().to_string())


if __name__ == "__main__":
    main()
