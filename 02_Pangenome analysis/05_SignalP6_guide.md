# SignalP 6.0 slow-sequential installation and SLURM array pipeline for signal peptide prediction

### This documentation captures everything that worked for us to install **SignalP6** and reliably run it on many protein FASTA files (CPU-mode).

> This pipeline ran on the Manitou HPC cluster (Université Laval) on CPU nodes via SLURM.

## 1. Prerequisites

- Linux HPC cluster with SLURM
- Conda (miniconda or mambaforge) available on the cluster
- SignalP 6.0 slow-sequential tarball: `signalp-6.0h.slow_sequential.tar.gz`
- Python 3.10 (do not use 3.12+)

> **Download:** SignalP 6.0 requires academic registration and is available upon request at:
> https://services.healthtech.dtu.dk/services/SignalP-6.0/

## 2. Create a clean environment (Python 3.10)

> SignalP 6 requires **torch 1.x**, which in turn expects **NumPy 1.x**. Don’t use Python 3.12+ for this.

```bash
conda create -n signalp6-310 python=3.10 pip -y
conda activate signalp6-310

# Make sure NumPy is 1.x (Torch 1.x expects this)
pip install 'numpy<2,>=1.23' --upgrade

# (Recommended) Install PyTorch 1.13 (CPU)
pip install 'torch==1.13.1'

```


## 3. Install SignalP 6 package

Unpack the official archive (`signalp-6.0g.fast.tar.gz`), then:

```bash
# Path where you unpacked the tarball
cd /path/to/signalp6_slow/signalp-6-package

# Install the Python package
pip install .

```


## 4. Place the model weights (slow-sequential mode)

```bash
SITE_DIR="/home/<USER>/miniconda3/envs/signalp6-310/lib/python3.10/site-packages/signalp"
mkdir -p "$SITE_DIR/model_weights"
cp -v /path/to/signalp-6-package/models/sequential_models_signalp6/* "$SITE_DIR/model_weights/"
```

## 5. Quick single-file test

```bash
conda activate signalp6-310

# Set MODEL_DIR
export MODEL_DIR="/path/to/sequential_models_signalp6_parent"

signalp6 --fastafile /path/to/input.fasta --organism euk --mode slow-sequential --format txt --model_dir "$MODEL_DIR" --output_dir /path/to/output_dir

```

## 6) Final production run — SLURM array (slow-sequential mode)

### Slow-sequential model weights must be extracted and MODEL_DIR set before submitting.

```bash
# Extract slow-sequential weights
mkdir -p ~/signalp6_models/_tmp_extract
tar -xzf ~/signalp-6.0h.slow_sequential.tar.gz -C ~/signalp6_models/_tmp_extract

# Detect MODEL_DIR (parent of 'sequential_models_signalp6')
FOUND=$(find ~/signalp6_models/_tmp_extract -type d -name 'sequential_models_signalp6' -print -quit)
export MODEL_DIR="$(dirname "$FOUND")"
ls -l "$MODEL_DIR/sequential_models_signalp6"
```

### For an array batch jobs for multiple isolates

```bash
#!/bin/bash
#SBATCH -J sp6_slowseq_arr7
#SBATCH -o sp6_slowseq_arr7-%A_%a.out
#SBATCH -e sp6_slowseq_arr7-%A_%a.err
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH -p small
#SBATCH --array=0-6   # matches the 7 FILES listed below

set -euo pipefail
module purge 2>/dev/null || true

# Set paths before running
MODEL_DIR="${MODEL_DIR:-/path/to/sequential_models_signalp6_parent}"
INDIR=/path/to/longest_isoform_proteins
OUTBASE=/path/to/signalp6_output

test -d "$MODEL_DIR/sequential_models_signalp6" || {
  echo "[FATAL] Missing: $MODEL_DIR/sequential_models_signalp6"; exit 2;
}

CONDA=/prg/miniconda/3-py3.10/bin/conda
ENV=signalp6_slow
RUN="$CONDA run -n $ENV signalp6"

mkdir -p "$OUTBASE"

# Edit this list to match your samples
FILES=(
  sample_1_proteins.fa
  sample_2_proteins.fa
  sample_3_proteins.fa
  sample_4_proteins.fa
  sample_5_proteins.fa
  sample_6_proteins.fa
  sample_7_proteins.fa
)

f="${INDIR}/${FILES[$SLURM_ARRAY_TASK_ID]}"
sample="${FILES[$SLURM_ARRAY_TASK_ID]%.*}"
outdir="${OUTBASE}/${sample}"
mkdir -p "$outdir"

$RUN \
  --fastafile "$f" \
  --organism euk \
  --mode slow-sequential \
  --format txt \
  --model_dir "$MODEL_DIR" \
  --output_dir "$outdir"
```

## Attribution
Teufel et al. (2022) SignalP 6.0 predicts all five types of signal peptides using protein language models.
*Nature Biotechnology* 40, 1023–1025. https://doi.org/10.1038/s41587-021-01156-3

