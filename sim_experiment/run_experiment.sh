#!/usr/bin/env bash
#
# End-to-end CAMUS simulated experiment:
#   1. simulate 18 samples of paired fastq.gz (2 cell types, 3 groups, 2 timepoints, 3 reps)
#   2. demultiplex every sample with CAMUS (group 3 is separated by cell-type barcode)
#   3. run MAGeCK RRA (initial vs D21) per cell type/group and compare group3 vs group1/2
#
# Run from the repo:  bash sim_experiment/run_experiment.sh
#
# Tunables (environment variables):
#   READS      reads per sample           (default 5000000)
#   WORKERS    parallel workers           (default 4)
#   FDR        significance threshold     (default 0.05)
#   MAGECK_BIN explicit path to mageck    (optional; else taken from PATH)
#   CONDA_ENV  conda env holding mageck   (optional; e.g. CONDA_ENV=mageck)
#
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$HERE/data"
DEMUX="$HERE/demux"
RESULTS="$HERE/results"

READS="${READS:-5000000}"
WORKERS="${WORKERS:-4}"
FDR="${FDR:-0.05}"

echo "=============================================================="
echo "CAMUS simulated experiment"
echo "  reads/sample : $READS"
echo "  workers      : $WORKERS"
echo "  FDR          : $FDR"
echo "=============================================================="

echo; echo "[1/3] Simulating fastq.gz ..."
python3 "$HERE/simulate_experiment.py" --out "$DATA" --reads-per-sample "$READS" --workers "$WORKERS"

echo; echo "[2/3] Demultiplexing with CAMUS ..."
python3 "$HERE/run_demux_all.py" --data "$DATA" --out "$DEMUX" --workers "$WORKERS"

echo; echo "[3/3] MAGeCK RRA + group comparison ..."
EXTRA=()
[ -n "${MAGECK_BIN:-}" ] && EXTRA+=(--mageck-binary "$MAGECK_BIN")
[ -n "${CONDA_ENV:-}" ] && EXTRA+=(--conda-env "$CONDA_ENV")
python3 "$HERE/analyze_experiment.py" --data "$DATA" --demux "$DEMUX" --out "$RESULTS" --fdr "$FDR" "${EXTRA[@]}"

echo
echo "=============================================================="
echo "Done."
echo "  fastq.gz        : $DATA/fastq_files/"
echo "  demux counts    : $DEMUX/"
echo "  results+report  : $RESULTS/REPORT.md"
echo "=============================================================="
