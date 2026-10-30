#!/bin/bash
# quantmsdiann v2.3.0 with DIA-NN 2.5.0, local executor on one node, calibration turned off so that
# every DIA-NN step uses the same fixed tolerances and scan window as the standalone run.
# Usage: quantmsdiann.sh <dataset>
set -euo pipefail
source $(dirname $0)/env.sh
DS=$1; dataset $DS
OUT=$W/results/$DS/quantmsdiann_dq
rm -rf $OUT/.nextflow* $OUT/results
mkdir -p $OUT $W/tmp/work_${DS}_dq $W/tmp/sing_cache
exec > $W/logs/quantmsdiann_${DS}_${SLURM_JOB_ID}.log 2>&1
export NXF_ASSETS=$W/tmp/nxf_home/assets NXF_SINGULARITY_CACHEDIR=$W/tmp/sing_cache NXF_SINGULARITY_LIBRARYDIR=$SING NXF_OPTS="-Xms1g -Xmx4g"
# Compute nodes have no internet: keep the default NXF_HOME (shared framework + plugin cache);
# only the pinned pipeline code lives in the iteration (NXF_ASSETS). nf-schema 2.6.1 needs Nextflow >= 25.10.
NXF=/hps/software/users/juan/pride/anaconda3/envs/nextflow/bin/nextflow
# Same predicted library as the standalone runs (the only step that is not bit-identical across machines)
SPECLIB=$W/results/$DS/standalone_r1/diann_report-lib.predicted.speclib
cat > $OUT/run.config <<CFG
process {
    executor = 'local'
    resourceLimits = [ cpus: 32, memory: '120.GB', time: '48.h' ]
    withLabel: 'diann' { container = '$DIANN_IMG' }
    withName: '.*(INSILICO_LIBRARY_GENERATION|PRELIMINARY_ANALYSIS|ASSEMBLE_EMPIRICAL_LIBRARY|INDIVIDUAL_ANALYSIS|FINAL_QUANTIFICATION)' { cpus = $THREADS; memory = '56.GB' }
}
singularity { enabled = true; autoMounts = true; runOptions = '-B /hps' }
docker.enabled = false
CFG
cd $OUT
start=$(date +%s)
# --direct-quant in every DIA-NN step (the pipeline adds it only to the final step when QuantUMS is off).
# Typed parameters: Nextflow >= 25.10 passes CLI values such as "--normalize false" as the string
# "false", which the pipeline's Groovy checks treat as true.
cat > $OUT/params.json <<JSON
{
  "input": "$W/scripts/sdrf/$DS.sdrf.tsv", "database": "$FASTA", "outdir": "$OUT/results", "speclib": "$SPECLIB",
  "diann_version": "2.5.0", "precursor_qvalue": 0.05, "matrix_qvalue": 0.01, "matrix_spec_q": 0.05,
  "mass_acc_automatic": false, "mass_acc_ms2": 20, "mass_acc_ms1": 10,
  "scan_window_automatic": false, "scan_window": $WINDOW, "quick_mass_acc": false,
  "min_pr_mz": 400, "max_pr_mz": $MAXPR, "min_fr_mz": 200, "max_fr_mz": 1800,
  "min_peptide_length": 7, "max_peptide_length": 30, "min_precursor_charge": 2, "max_precursor_charge": 4,
  "allowed_missed_cleavages": 2, "max_mods": 2, "met_excision": true, "pg_level": 2, "normalize": $NORM_NF,
  "use_quant": true, "quantums": false, "enable_pmultiqc": false,
  "extra_args": "--direct-quant"
}
JSON
$NXF run bigbio/quantmsdiann -r v2.3.0 -c $OUT/run.config -w $W/tmp/work_${DS}_dq -params-file $OUT/params.json -ansi-log false
echo "wall_clock_min $(( ($(date +%s) - start) / 60 ))"
