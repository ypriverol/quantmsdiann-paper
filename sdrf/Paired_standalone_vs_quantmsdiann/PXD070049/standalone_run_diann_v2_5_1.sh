#!/bin/bash

# ==============================================================================
# DIA-NN Standalone Benchmark Script (Strict Match with quantmsdiann Logs)
# ==============================================================================

BASE_DIR="/home/yueqx/Data_Disk/proteogenomics/quantmsdiann/tests/data_test/benchmarks/ProteoBench/PXD070049"
OUTPUT_DIR="${BASE_DIR}/test/diann/20260928_reanalyse/output"

mkdir -p "$OUTPUT_DIR"

RUN_LOG="diann_native_execution.log"
RESOURCE_LOG="diann_native_resources.log"

echo "=========================================================="
echo "Starting DIA-NN native benchmark test..."
echo "Logs are saved in: $RUN_LOG and $RESOURCE_LOG"
echo "=========================================================="

/usr/bin/time -v -o "$RESOURCE_LOG" docker run --rm \
    -v ${BASE_DIR}/raws/DIA-ZenoSWATH_mzML:/data/raw_data \
    -v ${BASE_DIR}/fasta:/data/fasta \
    -v ${OUTPUT_DIR}:/data/output \
    ghcr.io/bigbio/diann:2.5.1 \
    diann \
    --dir /data/raw_data \
    --fasta /data/fasta/ProteoBenchFASTA_MixedSpecies_HYE.fasta \
    --temp /data/output \
    --out /data/output/diann_report.tsv \
    --threads 18 \
    --fasta-search \
    --cut K*,R*,!*P \
    --missed-cleavages 2 \
    --fixed-mod Carbamidomethyl,57.021464,C \
    --var-mod Oxidation,15.994915,M \
    --var-mod Acetyl,42.010565,*n \
    --var-mods 2 \
    --mass-acc 20.0 \
    --mass-acc-ms1 10.0 \
    --min-pr-mz 400 \
    --max-pr-mz 2400 \
    --min-fr-mz 200 \
    --max-fr-mz 1800 \
    --min-pep-len 7 \
    --max-pep-len 30 \
    --min-pr-charge 2 \
    --max-pr-charge 4 \
    --predictor \
    --verbose 3 \
    --met-excision \
    --quick-mass-acc \
    --individual-mass-acc \
    --individual-windows \
    --rt-profiling \
    --pg-level 2 \
    --matrices \
    --qvalue 0.05 \
    --matrix-qvalue 0.01 \
    --matrix-spec-q 0.05 \
    --gen-spec-lib \
    --use-quant \
    --direct-quant \
    --reanalyse > "$RUN_LOG" 2>&1

echo "DIA-NN test run finished!"

# --qvalue: [precursor_qvalue] auto by diann_version (<2.5 -> 0.01, >=2.5 -> 0.05).