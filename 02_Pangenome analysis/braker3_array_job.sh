#!/bin/bash

#SBATCH -D /path/to/array_job
#SBATCH -J braker_array
#SBATCH -o braker_%A_%a.out
#SBATCH -e braker_%A_%a.err
#SBATCH -c 10
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --time=7-00:00:00
#SBATCH --mem=200G
#SBATCH --array=0-13  # For 14 samples

# module load singularity

# Export Augustus config path
export AUGUSTUS_CONFIG_PATH=/path/to/config

# Get sample name
SAMPLES_FILE="polishing_samples.txt"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLES_FILE)

# Define paths
BASE_DIR="/path/to/array_job"
WORKDIR="${BASE_DIR}/${sample}/braker3"
GENOME="${BASE_DIR}/${sample}/repeatmasker_output/${sample}_softmasked.fasta"
PROTEINS="/path/to/proteins.fa"
RNASEQ_DIR="/path/to/rnaseq_data/rnaseq"

# Create BRAKER3 working directory
mkdir -p "$WORKDIR"

# Run BRAKER3 with Singularity
singularity exec -B ${PWD}:${PWD} -B /path/to/config:/opt/Augustus/config /path/to/braker3.sif braker.pl \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --useexisting \
    --genome="$GENOME" \
    --species="$sample" \
    --prot_seq="$PROTEINS" \
    --rnaseq_sets_ids=S001A93_ID1,S001E7A_ID2,21dpi_ID3,7dpi_ID4 \
    --rnaseq_sets_dirs="$RNASEQ_DIR" \
    --workingdir="$WORKDIR" \
    --threads=10 \
    --busco_lineage=eukaryota_odb10
