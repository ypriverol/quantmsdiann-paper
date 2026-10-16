"""Turn the crawled file manifest into the PRIDE submission and the Zenodo file list.

Reads pride_submission/submission_files.tsv (from scripts.build_submission_manifest)
and writes, under pride_submission/:

  sdrf/                      SDRF copies to upload: byte-identical to the as-run SDRF
                             snapshot in this repo (sdrf/)
  submission_files_final.tsv every file to upload (FTP URL or local path), with
                             its PRIDE type and whether it also goes to Zenodo
  raw_files.raw.zip          raw_files.txt zipped, uploaded as the RAW-type file
                             (plus one <PXD>.qpx.zip row per run with a qpx/ export on the FTP;
                             the zip is built by scripts.stage_pride_submission)
  zip_members.tsv            what goes into each <dataset>_<version>.zip run bundle (one per
                             run: the run folder from the FTP, built by stage_pride_submission)
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
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from scripts.build_run_metadata import run_dir
from scripts.build_submission_manifest import SKIP, crawl, head_size, listing

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "pride_submission"
REVIEW = REPO / "data/proteobench-review-experiments/x"
SDRF_DIR = REPO / "sdrf"

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
    "quantmsdiann runs (PXD049412, PXD070049). Raw files are not re-uploaded: raw_files.raw.zip (raw_files.txt) and "
    "raw_files.tsv list every raw file processed with its URL in the original deposit. Analysis code: https://github.com/ypriverol/quantmsdiann-paper."
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


def qpx_files(rows: list[dict]) -> list[dict]:
    """One <PXD>.qpx.zip per run that has a qpx/ export on the FTP; zipped by stage_pride_submission.

    A project with qpx in more than one run gets the version in the name: <PXD>_<version>.qpx.zip.
    """
    found = []
    for (ds, ver), base in sorted({(r["dataset"], r["diann_version"]): r["source_url"].split("/quant_tables/")[0]
                                   for r in rows if "/quant_tables/" in str(r["source_url"])}.items()):
        try:
            names = [n for n in listing(f"{base}/qpx/", tries=2) if not n.endswith("/")]
        except Exception:  # no qpx/ for this run
            continue
        acc = re.search(r"(PXD\d+|MSV\d+)", ds).group(1)
        found.append((acc, ds, ver, f"{base}/qpx/", sum(head_size(f"{base}/qpx/{n}") for n in names)))
    per_acc = pd.Series([f[0] for f in found]).value_counts()
    return [{"dataset": ds, "diann_version": ver, "file_type": "RESULT",
             "upload_name": f"{acc}.qpx.zip" if per_acc[acc] == 1 else f"{acc}_{ver}.qpx.zip",
             "source_url": url, "size_bytes": size, "note": "qpx export (zipped at staging; size = sum of members)"}
            for acc, ds, ver, url, size in found]


def bundle_runs(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Replace each run's SEARCH/RESULT/OTHER files with one <dataset>_<version>.zip for PRIDE.

    For a run published on the FTP the zip holds the whole run folder (quant_tables/, pipeline_info/,
    sdrf/, pmultiqc/, run_metadata.json, nextflow.log), minus SKIP (MSstats input, spectral library,
    .quant, work/) and qpx/ (its own zip). Local extras (recovered params, sweep traces, paired runs)
    go in as well. The global files (runs.tsv, extra.config, Table S6, raw_files.tsv) go into
    run_metadata.zip. SDRFs, qpx zips and raw_files.raw.zip stay separate. Returns (rows, zip members).
    """
    def run_key(r):
        return (r["dataset"], r["diann_version"])
    # The global OTHER files (runs.tsv, extra.config, Table S6, raw_files.tsv) form one more bundle.
    per_run = [r for r in rows if r["file_type"] in ("SEARCH", "RESULT", "OTHER")
               and not r["upload_name"].endswith(".qpx.zip")]
    keep = [r for r in rows if r not in per_run]
    bases = {run_key(r): r["source_url"].split("/quant_tables/")[0]
             for r in per_run if str(r["source_url"]).startswith("http") and "/quant_tables/" in r["source_url"]}
    old_name = {r["source_url"]: r["upload_name"] for r in per_run}
    members = []
    for key in sorted({run_key(r) for r in per_run}):
        zname = "run_metadata.zip" if key[0] == "all" else f"{key[0]}_{key[1]}.zip"
        mine = [r for r in per_run if run_key(r) == key]
        prefix = "" if key[0] == "all" else f"{key[0]}_{key[1]}_"
        if key in bases:
            base = bases[key]
            rels = [rel for rel in crawl(base) if not SKIP.search("/" + rel)]
            with ThreadPoolExecutor(8) as ex:
                sizes = list(ex.map(lambda rel: head_size(f"{base}/{rel}"), rels))
            for rel, size in zip(rels, sizes):
                url = f"{base}/{rel}"
                members.append({"zip": zname, "member": rel, "source": url, "size_bytes": size,
                                "old_upload_name": old_name.get(url, "")})
            local = [r for r in mine if not str(r["source_url"]).startswith("http")]
            for r in local:
                members.append({"zip": zname, "member": "recovered/" + r["upload_name"].removeprefix(prefix),
                                "source": r["source_url"], "size_bytes": r["size_bytes"],
                                "old_upload_name": r["upload_name"]})
        else:
            for r in mine:
                members.append({"zip": zname, "member": r["upload_name"].removeprefix(prefix),
                                "source": r["source_url"], "size_bytes": r["size_bytes"],
                                "old_upload_name": r["upload_name"]})
        ftype = "SEARCH" if any(r["file_type"] == "SEARCH" for r in mine) else "OTHER"
        keep.append({"dataset": key[0], "diann_version": key[1], "file_type": ftype, "upload_name": zname,
                     "source_url": f"zip:{zname}", "size_bytes": sum(m["size_bytes"] for m in members if m["zip"] == zname),
                     "note": "run bundle (zipped at staging; size = sum of members)"})
    return keep, members


def main() -> None:
    d = pd.read_csv(OUT / "submission_files.tsv", sep="\t")
    sdrf_dir = OUT / "sdrf"
    shutil.rmtree(sdrf_dir, ignore_errors=True)
    sdrf_dir.mkdir()
    runs = pd.read_csv(SDRF_DIR / "runs.tsv", sep="\t", dtype=str)
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
        mdir = run_dir(ds, ver)
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
    for f in (SDRF_DIR / "runs.tsv", SDRF_DIR / "extra.config",
              REPO / "analysis/figures/review/supplementary_revision/TableS6_run_parameters.tsv"):
        rows.append({"dataset": "all", "diann_version": "", "file_type": "OTHER", "upload_name": f.name,
                     "source_url": str(f), "size_bytes": f.stat().st_size, "note": "run parameters for all runs"})
    rows += review_files()
    rows += qpx_files(d.to_dict("records"))
    # PRIDE needs at least one RAW-type file: the raw file list goes in as raw_files.raw.zip.
    raw_zip = OUT / "raw_files.raw.zip"
    with zipfile.ZipFile(raw_zip, "w", zipfile.ZIP_DEFLATED) as z:
        info = zipfile.ZipInfo("raw_files.txt", date_time=(2026, 1, 1, 0, 0, 0))  # fixed date: same bytes on rebuild
        info.compress_type = zipfile.ZIP_DEFLATED
        z.writestr(info, (OUT / "raw_files.txt").read_bytes())
    for ftype, p in (("RAW", raw_zip), ("OTHER", OUT / "raw_files.tsv")):
        rows.append({"dataset": "all", "diann_version": "", "file_type": ftype, "upload_name": p.name,
                     "source_url": str(p), "size_bytes": p.stat().st_size, "note": "raw files used, with original URLs"})
    for r in rows:
        r.setdefault("note", "")
        # qpx zips are built at PRIDE staging; zenodo_upload streams single FTP files, so PRIDE only.
        r["zenodo"] = "no" if r["upload_name"] in ZENODO_EXCLUDE or r["upload_name"].endswith(".qpx.zip") else "yes"
    flat = pd.DataFrame(rows)
    assert not flat.upload_name.duplicated().any(), flat[flat.upload_name.duplicated()]
    # Zenodo keeps the flat files (zenodo_upload streams single files); PRIDE gets one zip per run.
    flat[flat.zenodo == "yes"].to_csv(OUT / "zenodo_files.tsv", sep="\t", index=False)
    rows, members = bundle_runs(rows)
    for r in rows:
        r.setdefault("zenodo", "no")
    df = pd.DataFrame(rows)
    assert not df.upload_name.duplicated().any(), df[df.upload_name.duplicated()]
    df.to_csv(OUT / "submission_files_final.tsv", sep="\t", index=False)
    pd.DataFrame(members).to_csv(OUT / "zip_members.tsv", sep="\t", index=False)

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
    zen = flat[flat.zenodo == "yes"]
    print(f"PRIDE: {len(df)} files, {tot:.1f} GB ({len(members)} files inside {len(df[df.source_url.str.startswith('zip:')])} "
          f"run zips) | Zenodo: {len(zen)} files, {zen.size_bytes.sum() / 1e9:.1f} GB")
    print("instruments:", inst)
    print(df.note.value_counts().to_string())


if __name__ == "__main__":
    main()
