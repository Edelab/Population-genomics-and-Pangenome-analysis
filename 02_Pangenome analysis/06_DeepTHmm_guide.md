# DeepTMHMM GPU array pipeline for transmembrane topology prediction on SignalP-processed proteins

> This pipeline was run on the Manitou HPC cluster (Université Laval) 
> using GPU nodes via SLURM. The DeepTMHMM Singularity image 
> (`deeptmhmm-1.0.24.sif`) was obtained from:
> https://biolib.com/DTU/DeepTMHMM/
> Input protein FASTA files were pre-processed using SignalP 6.0 

```bash
#!/bin/bash
#SBATCH -J pb_deepthmm
#SBATCH -D /path/to/directory/effector_prediction
#SBATCH -o Pb-%A_%a.out
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH -p small
#SBATCH --time=1-00:00
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=email@ulaval.ca

# === USER SETTINGS ===
INPUT_DIR="/path/to/directory/slow_signalp_processed_proteins"
OUT_BASE="/path/to/directory/outputs_thmm"
PATTERN="*_processed_entries.fasta"
SIF="/prg/singularity/images/deeptmhmm-1.0.24.sif"

# optional: tiny CPU threads
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

# === DISCOVER FILES ===
mkdir -p "$OUT_BASE"
mapfile -t FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "$PATTERN" | sort)
N=${#FILES[@]}
if (( N == 0 )); then
  echo "[ERROR] No files matched '$PATTERN' in $INPUT_DIR" >&2
  exit 2
fi

# require array submission
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=1-$N (don’t run with bash)}"
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if (( IDX < 0 || IDX >= N )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID out of range (have $N files)" >&2
  exit 3
fi

FASTA="${FILES[$IDX]}"
SAMPLE="$(basename "$FASTA")"
# Assumes filenames end in _processed_entries.fasta
SAMPLE="${SAMPLE%_processed_entries.fasta}"

mkdir -p "$OUT_BASE/$SAMPLE"
cd "$OUT_BASE/$SAMPLE"

echo "[INFO] Task $SLURM_ARRAY_TASK_ID/$N  Sample=$SAMPLE"
singularity run --nv -B /project:/project "$SIF" --fasta "$FASTA"


# command used to run this - sbatch --array=1-55 deepthmm_array.sh

```

## Attribution
Hallgren et al. (2022) DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks.
*bioRxiv*. https://doi.org/10.1101/2022.04.08.487609
