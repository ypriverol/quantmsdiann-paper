#!/bin/bash
# Standalone DIA-NN 2.5.0: one command over all runs (in-silico library + two passes with --reanalyse).
# Files are passed with --f in name order (as quantmsdiann sorts them), not --dir: the order changes the empirical library.
# Usage: standalone.sh <dataset> <replicate>
set -euo pipefail
source $(dirname $0)/env.sh
DS=$1; REP=$2; dataset $DS
OUT=$W/results/$DS/standalone_sorted_r$REP
mkdir -p $OUT
exec > $W/logs/standalone_sorted_${DS}_r${REP}_${SLURM_JOB_ID}.log 2>&1
FILES=$(ls $IN | grep -E '\.(raw|mzML)$' | LC_ALL=C sort | sed "s|^|--f $IN/|" | tr '\n' ' ')
start=$(date +%s)
singularity exec -B /hps $DIANN_IMG diann \
  $FILES --fasta $FASTA --temp $OUT --out $OUT/diann_report.tsv --threads $THREADS \
  --fasta-search --cut 'K*,R*,!*P' --missed-cleavages 2 \
  --fixed-mod Carbamidomethyl,57.021464,C --var-mod Oxidation,15.994915,M --var-mod 'Acetyl,42.010565,*n' --var-mods 2 \
  --mass-acc 20 --mass-acc-ms1 10 --window $WINDOW \
  --min-pr-mz 400 --max-pr-mz $MAXPR --min-fr-mz 200 --max-fr-mz 1800 --min-pep-len 7 --max-pep-len 30 \
  --min-pr-charge 2 --max-pr-charge 4 --predictor --verbose 3 --met-excision --rt-profiling --pg-level 2 $NORM \
  --matrices --qvalue 0.05 --matrix-qvalue 0.01 --matrix-spec-q 0.05 --gen-spec-lib --use-quant --direct-quant --reanalyse
echo "wall_clock_min $(( ($(date +%s) - start) / 60 ))"
