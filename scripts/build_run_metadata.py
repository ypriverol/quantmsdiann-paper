"""Snapshot the SDRFs and the final pipeline parameters of every run in the manuscript.

Writes (committed to the repo, so later edits to sdrf-annotated-datasets or the
FTP do not change what the paper used):

  sdrf/<dataset>/<file>.sdrf.tsv           as-run SDRFs, verbatim (one copy per distinct content)
  sdrf/<sdrf folder>/<run>/
      params.json                          final Nextflow parameters (pipeline_info/params_*.json)
      launch.txt                           launch command (pipeline_report.txt) or SLURM run header
      diann_config.cfg                     DIA-NN options derived from the SDRF
  sdrf/extra.config                        resource overrides passed with -c to every cluster run
  sdrf/runs.tsv                            index: run -> SDRF, sources, key parameters (Table S6)

Sources, in order: the run's own folder on the PRIDE FTP; for runs whose FTP
folder has no sdrf/ or pipeline_info/ (phosphoproteomics, PXD071075 sweep), the
files copied from the cluster run folders (data/_cluster_run_files/, fetched once
with ssh; see sdrf/README.md); for the paired standalone/quantmsdiann
runs, data/proteobench-review-experiments/x/.

Usage:
    python -m scripts.build_run_metadata
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import urllib.request
from pathlib import Path

import pandas as pd

from scripts.build_submission_manifest import FTP, RUNS, crawl

REPO = Path(__file__).resolve().parents[1]
SDRF_DIR = REPO / "sdrf"
CL = REPO / "data/_cluster_run_files"
REVIEW = REPO / "data/proteobench-review-experiments/codon"

# Runs with no sdrf/ or pipeline_info/ on the FTP: (sdrf, launch log, params) in _cluster/.
P = "quantms-datasets/projects/phospho-proteomics"
LOGS = "logs/phospho-proteomics"
CLUSTER_RUNS = {
    ("PXD034128_biological_study", "v2_5_1"): (f"{P}/PXD034128/PXD034128-phospho-biological-study.sdrf.tsv",
                                               f"{LOGS}/PXD034128/v2_5_1/phospho-biological-study/slurm_49181394.out"),
    ("PXD034128_biological_study", "v2_5_1_enterprise"): (f"{P}/PXD034128/PXD034128-phospho-biological-study.sdrf.tsv",
                                                          f"{LOGS}/PXD034128/v2_5_1_enterprise/phospho-biological-study/slurm_49181395.out"),
    ("PXD034128_highspeed_DIA", "v2_5_1"): (f"{P}/PXD034128/PXD034128-phospho-highspeed-DIA.sdrf.tsv",
                                            f"{LOGS}/PXD034128/v2_5_1/phospho-highspeed-DIA/slurm_49181396.out"),
    ("PXD034128_highspeed_DIA", "v2_5_1_enterprise"): (f"{P}/PXD034128/PXD034128-phospho-highspeed-DIA.sdrf.tsv",
                                                       f"{LOGS}/PXD034128/v2_5_1_enterprise/phospho-highspeed-DIA/slurm_49181398.out"),
    ("PXD034623", "v2_5_1"): (f"{P}/PXD034623/PXD034623-M2-Galectin1-phospho-DIA.sdrf.tsv",
                              f"{LOGS}/PXD034623/v2_5_1/M2-Galectin1-phospho-DIA/slurm_49181399.out"),
    ("PXD034623", "v2_5_1_enterprise"): (f"{P}/PXD034623/PXD034623-M2-Galectin1-phospho-DIA.sdrf.tsv",
                                         f"{LOGS}/PXD034623/v2_5_1_enterprise/M2-Galectin1-phospho-DIA/slurm_49181400.out"),
    ("PXD049692", "v2_5_1"): (f"{P}/PXD049692/PXD049692-NK-fibrin-IL15-phospho-diaPASEF.sdrf.tsv",
                              f"{LOGS}/PXD049692/v2_5_1/NK-fibrin-IL15-phospho-diaPASEF/slurm_49181401.out"),
    **{("PXD071075", f"v2_5_0_sweep_{q}cores"): (
        "quantms-test-datasets/benchmarks/dia/OrbitrapEclipse/PXD071075/PXD071075.sdrf.tsv",
        f"logs/PXD071075/v2_5_0_sweep_{q}cores/slurm_{j}.out")
       for q, j in (("010", 27462633), ("050", 28039038), ("100", 27477233), ("200", 27903443), ("300", 27489842))},
}
KEYS = ["diann_version", "database", "allowed_missed_cleavages", "min_peptide_length", "max_peptide_length",
        "min_precursor_charge", "max_precursor_charge", "min_pr_mz", "max_pr_mz", "min_fr_mz", "max_fr_mz",
        "max_mods", "variable_mods", "mod_localization", "precursor_mass_tolerance", "fragment_mass_tolerance",
        "mass_acc_automatic", "pg_level", "normalize", "matrix_qvalue", "matrix_spec_q", "met_excision"]
EXTRA_ARG = re.compile(r"--(\w+)\s+(\S+)")


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def run_dir(dataset: str, run: str) -> Path:
    """Folder for one run's parameter files, next to its SDRF in sdrf/."""
    for prefix, folder in (("ProteoBench-", "ProteoBench"), ("Paired-", "Paired_standalone_vs_quantmsdiann"),
                           ("PXD034128_", "PXD034128")):
        if dataset.startswith(prefix):
            base = SDRF_DIR / folder / dataset[len(prefix):]
            break
    else:
        base = SDRF_DIR / dataset
    return base if dataset.startswith("Paired-") else base / run


def store_sdrf(dataset: str, name: str, content: bytes, seen: dict) -> str:
    """Keep one verbatim copy per distinct SDRF content; return its repo path."""
    h = hashlib.sha256(content).hexdigest()
    if h not in seen:
        folder = ("ProteoBench" if dataset.startswith("ProteoBench-") else
                  "Paired_standalone_vs_quantmsdiann" if dataset.startswith("Paired-") else dataset.split("_")[0])
        base = SDRF_DIR / folder / name
        if base.exists() and hashlib.sha256(base.read_bytes()).hexdigest() != h:
            base = base.with_name(base.name.replace(".sdrf.tsv", f"_{h[:8]}.sdrf.tsv"))
        base.parent.mkdir(parents=True, exist_ok=True)
        base.write_bytes(content)
        seen[h] = str(base.relative_to(REPO))
    return seen[h]


def launch_from_report(text: str) -> str:
    m = re.search(r"The command used to launch the workflow was as follows:\s*\n(.+?)\n\s*\n", text, re.S)
    return m.group(1).strip() if m else ""


def header_from_log(text: str) -> str:
    keep = [l for l in text.splitlines()[:40] if re.match(r"(=+|SLURM job|Pipeline dir|SDRF|FASTA|RAW|Work dir|Results dir|"
                                                        r"DIA-NN version|Profiles|Input type|Extra NF|N E X T F L O W|Launching)", l)]
    return "\n".join(keep)


def key_params(params: dict) -> dict:
    out = {k: params.get(k, "") for k in KEYS}
    out["database"] = Path(str(out["database"])).name if out["database"] else ""
    return out


def from_log(header: str) -> dict:
    """Run-specific parameters recoverable from a SLURM header (others are pipeline defaults)."""
    out = {k: "" for k in KEYS}
    for line in header.splitlines():
        if line.startswith("FASTA"):
            out["database"] = line.rsplit("/", 1)[-1]
        if line.startswith("DIA-NN version"):
            out["diann_version"] = line.split(":", 1)[1].strip().replace("_", ".", 2).replace("_enterprise", "-enterprise")
        if line.startswith("Extra NF args"):
            for k, v in EXTRA_ARG.findall(line):
                if k in out:
                    out[k] = v
    return out


def sdrf_mods(path: str) -> str:
    """Modifications declared in the SDRF, which is what quantmsdiann passes to DIA-NN."""
    s = pd.read_csv(REPO / path, sep="\t", dtype=str)
    mods = set()
    for c in [c for c in s.columns if "modification parameters" in c.lower()]:
        for v in s[c].dropna().unique():
            f = dict(x.split("=", 1) for x in v.split(";") if "=" in x)
            if "NT" in f:
                site = f.get("TA") or f.get("PP", "")
                mods.add(f"{f['NT']} ({site}, {f.get('MT', '').lower()})")
    return "; ".join(sorted(mods))


def table_s6(df: pd.DataFrame) -> None:
    """Supplementary Table S6: one row per run, the parameters a reader needs to reproduce it."""
    rng = lambda a, b: "" if not (str(a) and str(b)) or "nan" in (str(a), str(b)) else f"{a}-{b}"
    t = pd.DataFrame({
        "Dataset": df["dataset"], "Run": df["run"], "DIA-NN": df["diann_version"],
        "FASTA": df["database"],
        "Missed cleavages": df["allowed_missed_cleavages"],
        "Peptide length": [rng(a, b) for a, b in zip(df["min_peptide_length"], df["max_peptide_length"])],
        "Precursor charge": [rng(a, b) for a, b in zip(df["min_precursor_charge"], df["max_precursor_charge"])],
        "Precursor m/z": [rng(a, b) for a, b in zip(df["min_pr_mz"], df["max_pr_mz"])],
        "Fragment m/z": [rng(a, b) for a, b in zip(df["min_fr_mz"], df["max_fr_mz"])],
        "Max variable mods": df["max_mods"], "Modifications (from SDRF)": df["sdrf"].map(sdrf_mods),
        "Normalisation": df["normalize"],
        "SDRF (this repo)": df["sdrf"], "Parameters from": df["params_source"],
    }).fillna("")
    t = t.replace({"nan": ""})
    t.loc[t["Parameters from"].str.startswith("SLURM"), "Parameters from"] = (
        "SLURM run log (run-specific arguments shown; blank = pipeline default)")
    out = REPO / "analysis/figures/review/supplementary_revision/TableS6_run_parameters.tsv"
    t.to_csv(out, sep="\t", index=False)


def main() -> None:
    if not CL.exists():
        raise SystemExit("data/_cluster_run_files/ missing: fetch it first (see sdrf/README.md)")
    seen: dict[str, str] = {}
    rows = []
    for ds, ver, run in RUNS:
        if ds == "PXD071075":
            continue  # the paper's PXD071075 runs are the v2.5.0 sweep (below)
        out = run_dir(ds, ver)
        out.mkdir(parents=True, exist_ok=True)
        base = f"{FTP}/{run}"
        files = crawl(base)
        sdrfs = [f for f in files if f.startswith("sdrf/") and f.endswith(".sdrf.tsv")]
        params = sorted(f for f in files if re.match(r"pipeline_info/params_.*\.json$", f))
        row = {"dataset": ds, "run": ver, "ftp_folder": run}
        if sdrfs:
            row["sdrf"] = store_sdrf(ds, sdrfs[0].split("/")[-1], fetch(f"{base}/{sdrfs[0]}"), seen)
            row["sdrf_source"] = "as-run (FTP sdrf/)"
        if "sdrf/diann_config.cfg" in files:
            (out / "diann_config.cfg").write_bytes(fetch(f"{base}/sdrf/diann_config.cfg"))
        if params:
            p = json.loads(fetch(f"{base}/{params[-1]}"))
            (out / "params.json").write_text(json.dumps(p, indent=2))
            row.update(key_params(p))
            row["params_source"] = f"FTP {params[-1]}"
        if "pipeline_info/pipeline_report.txt" in files:
            cmd = launch_from_report(fetch(f"{base}/pipeline_info/pipeline_report.txt").decode("utf-8", "replace"))
            (out / "launch.txt").write_text(cmd + "\n")
        if (ds, ver) in CLUSTER_RUNS:
            sdrf, log = CLUSTER_RUNS[(ds, ver)]
            if "sdrf" not in row:
                row["sdrf"] = store_sdrf(ds, Path(sdrf).name, (CL / sdrf).read_bytes(), seen)
                row["sdrf_source"] = "as-run input (cluster quantms-datasets/projects)"
            if "params_source" not in row:
                hdr = header_from_log((CL / log).read_text(errors="replace"))
                (out / "launch.txt").write_text(hdr + "\n")
                row.update(from_log(hdr))
                row["params_source"] = f"SLURM run log {Path(log).name} (run-specific args; rest = pipeline defaults)"
        rows.append(row)
        print(ds, ver, row.get("sdrf_source", "NO SDRF"), "|", row.get("params_source", "NO PARAMS"), flush=True)

    for (ds, ver), (sdrf, log) in CLUSTER_RUNS.items():
        if ds != "PXD071075":
            continue
        out = run_dir(ds, ver)
        out.mkdir(parents=True, exist_ok=True)
        hdr = header_from_log((CL / log).read_text(errors="replace"))
        (out / "launch.txt").write_text(hdr + "\n")
        row = {"dataset": ds, "run": ver, "ftp_folder": "quantmsdiann-benchmarks/PXD071075_cluster_sizes",
               "sdrf": store_sdrf(ds, Path(sdrf).name, (CL / sdrf).read_bytes(), seen),
               "sdrf_source": "as-run input (cluster quantms-test-datasets)",
               **from_log(hdr), "params_source": f"SLURM run log {Path(log).name} (run-specific args; rest = pipeline defaults)"}
        rows.append(row)

    for ds in ("PXD049412", "PXD070049"):  # paired standalone vs quantmsdiann (Supplementary Note 8)
        # Codon iteration 2026-10-28_01_run_paired-standalone-vs-quantmsdiann (DIA-NN 2.5.0, quantmsdiann v2.3.0)
        qm = REVIEW / ds / "quantmsdiann"
        out = run_dir(f"Paired-{ds}", "")
        out.mkdir(parents=True, exist_ok=True)
        for f in ("standalone_sorted.sh", "quantmsdiann_dq.sh", "env.sh"):
            shutil.copy2(REVIEW / "scripts" / f, out / f)
        shutil.copy2(qm / "run.config", out / "quantmsdiann_run.config")
        p = sorted((qm / "results/pipeline_info").glob("params_*.json"))[-1]
        params = json.loads(p.read_text())
        (out / "params.json").write_text(json.dumps(params, indent=2))
        shutil.copy2(qm / "results/sdrf/diann_config.cfg", out / "diann_config.cfg")
        shutil.copy2(qm / f"results/sdrf/{ds}.sdrf.tsv", out / "as_run.sdrf.tsv")
        # The deposited design is the annotated SDRF; as_run.sdrf.tsv is the pipeline input (cluster paths).
        sdrf = SDRF_DIR / f"Paired_standalone_vs_quantmsdiann/{ds}_local.sdrf.tsv"
        rows.append({"dataset": f"Paired-{ds}", "run": "quantmsdiann v2_3_0 + standalone DIA-NN 2.5.0",
                     "ftp_folder": "", "sdrf": store_sdrf(f"Paired-{ds}", sdrf.name, sdrf.read_bytes(), seen),
                     "sdrf_source": "annotated; pipeline input in as_run.sdrf.tsv", **key_params(params),
                     "params_source": f"{p.name}; standalone command in standalone_sorted.sh"})

    shutil.copy2(CL / "quantms-datasets/scripts/extra.config", SDRF_DIR / "extra.config")
    df = pd.DataFrame(rows)
    # FASTAs are stored once, gzipped, in sdrf/fasta/ (checksums in SHA256SUMS.txt).
    df.insert(df.columns.get_loc("database") + 1, "fasta", "sdrf/fasta/" + df["database"] + ".gz")
    df.to_csv(SDRF_DIR / "runs.tsv", sep="\t", index=False)
    table_s6(df)
    print(f"\n{len(df)} runs; {len(seen)} distinct SDRFs; missing SDRF: {df.sdrf.isna().sum()}; "
          f"missing params: {df.params_source.isna().sum()}")


if __name__ == "__main__":
    main()
