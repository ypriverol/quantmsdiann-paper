# Shared settings for the paired standalone vs quantmsdiann comparison.
W=/hps/nobackup/juan/pride/reanalysis/quantms-work/quantmsdiann/2026-10-28_01_run_paired-standalone-vs-quantmsdiann
SING=/hps/nobackup/juan/pride/reanalysis/singularity
DIANN_IMG=$SING/ghcr.io-bigbio-diann-2.5.0.img
WIFF_IMG=$SING/ghcr.io-bigbio-wiffconverter-sif-0.10.img
THREADS=16
# Per dataset: input folder, FASTA, max precursor m/z, scan window, normalisation flag
dataset() {
  case $1 in
    PXD049412) IN=$W/tmp/raw/PXD049412; FASTA=$W/tmp/fasta/ProteoBenchFASTA_DDAQuantification_noecoli.fasta
               MAXPR=800; WINDOW=6; NORM=--no-norm; NORM_NF=false ;;
    PXD070049) IN=$W/tmp/mzml/PXD070049; FASTA=$W/tmp/fasta/ProteoBenchFASTA_MixedSpecies_HYE.fasta
               MAXPR=2400; WINDOW=5; NORM=; NORM_NF=true ;;
  esac
}
