#!/bin/bash
#SBATCH -J sp6_slowseq_array
#SBATCH -o sp6_slowseq_array-%A_%a.out
#SBATCH -e sp6_slowseq_array-%A_%a.err
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH -p small
#SBATCH --array=0-54   # set to (total number of FASTA files - 1)

set -euo pipefail
module purge 2>/dev/null || true

# ======== SET PATHS ========
MODEL_DIR="${MODEL_DIR:-/path/to/sequential_models_signalp6_parent}"
INDIR=/path/to/longest_isoform_proteins
OUTBASE=/path/to/signalp6_output
CONDA=/path/to/miniconda3/bin/conda
ENV=signalp6_slow
# =================================================

# Validate MODEL_DIR
test -d "$MODEL_DIR/sequential_models_signalp6" || {
  echo "[FATAL] Missing: $MODEL_DIR/sequential_models_signalp6"
  exit 2
}

# Input directory detections
mapfile -t FILES < <(find "$INDIR" -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.faa" \) | sort)
TOTAL=${#FILES[@]}
echo "[INFO] Found $TOTAL FASTA files"
echo "[INFO] Running array index: $SLURM_ARRAY_TASK_ID / $((TOTAL - 1))"

# Guard against out-of-range index
if [ "$SLURM_ARRAY_TASK_ID" -ge "$TOTAL" ]; then
  echo "[FATAL] Array index $SLURM_ARRAY_TASK_ID exceeds file count $TOTAL"
  exit 1
fi

# Set input/output for array task
f="${FILES[$SLURM_ARRAY_TASK_ID]}"
sample="$(basename "$f")"
sample="${sample%.*}"
outdir="${OUTBASE}/${sample}"
mkdir -p "$outdir"

echo "[RUN] Sample : $sample"
echo "      FASTA  : $f"
echo "      OUTDIR : $outdir"

# SignalP6 Run slow-sequential
"$CONDA" run -n "$ENV" signalp6 \
  --fastafile "$f" \
  --organism euk \
  --mode slow-sequential \
  --format txt \
  --model_dir "$MODEL_DIR" \
  --output_dir "$outdir"

echo "[DONE] $sample"
