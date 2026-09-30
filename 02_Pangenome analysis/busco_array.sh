#!/bin/bash

#SBATCH -D /path/to/directory/array_job
#SBATCH -J busco
#SBATCH -o logs/busco_%A_%a.out
#SBATCH -e logs/busco_%A_%a.err
#SBATCH -c 6
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G
#SBATCH --array=0-13  # Adjust based on your sample count

# Load modules
module load miniconda/3-py3.10 busco/5.8.2

# Initialize conda environment
eval "$(conda shell.bash hook)"
conda activate busco-5.8.2

# Make sure log directory exists
mkdir -p logs

# Read the sample name from file
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" polishing_samples.txt)

# Set paths
basedir="${sample}"
input_fasta="${basedir}/hapog_3/hapog_results/${sample}_hapog_3.fasta"
output_base="${basedir}/busco_output"
output_name="${sample}_busco"
lineage_path="/path/to/directory/busco_downloads/lineages/eukaryota_odb10"

# Make sure output directory exists
mkdir -p "$output_base"

# Check if input exists
if [[ ! -f "$input_fasta" ]]; then
  echo "ERROR: Input genome file $input_fasta not found for $sample"
  exit 1
fi

# Run BUSCO
echo "Running BUSCO for $sample ..."
busco -i "$input_fasta" \
      -o "$output_name" \
      --out_path "$output_base" \
      -l "$lineage_path" \
      --mode genome \
      -c 6
