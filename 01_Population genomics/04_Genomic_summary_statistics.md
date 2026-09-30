# Pixy complete pipeline for nucleotide diversity (π) and absolute divergence (Dxy) analysis


## 1. Installation

```bash
# Create and activate conda environment
conda create -n pixy_env python=3.8 -y
conda activate pixy_env

# Install pixy via conda-forge
conda install -c conda-forge pixy -y

# Verify installation
pixy --version

# Install htslib tools (required for VCF indexing)
conda install -c bioconda htslib samtools bcftools -y
```


## 2. Input File Preparation

Pixy requires an **all-sites VCF** (invariant + variant sites).


## 3. Running Pixy

```bash

# Computes pi, and dxy only because fst didn't work for pixy if data is in haploid mode 
# Window size: 10,000 bp (10 kb)

pixy --stats pi dxy \
     --vcf all_sites.vcf.gz \
     --populations Pb_populations.txt \
     --window_size 10000 \
     --n_cores 8 \
     --output_folder ./pixy_out \
     --output_prefix Pb_genome_masked

# Output files produced:
#   ./pixy_out/Pb_genome_masked_pi.txt
#   ./pixy_out/Pb_genome_masked_dxy.txt

```



## 4. Analysis and Visualization Scripts


```bash

# Mean pi per cluster

echo "=== Mean pi per cluster (pixy recommended) ==="
awk 'NR>1 && $7!="NA" && $8!="NA" {
    diffs[$1] += $7
    comps[$1] += $8
}
END {
    for (c in diffs)
        printf "%s\tmean_pi=%.8f\n", c, diffs[c]/comps[c]
}' ./path/to/Pb_genome_masked_pi.txt

# Mean dxy between cluster pairs

echo "=== Mean dXY between clusters (pixy recommended) ==="
awk 'NR>1 && $7!="NA" && $8!="NA" {
    diffs[$1"_vs_"$2] += $7
    comps[$1"_vs_"$2] += $8
}
END {
    for (p in diffs)
        printf "%s\tmean_dxy=%.8f\n", p, diffs[p]/comps[p]
}' ./path/to/Pb_genome_masked_dxy.txt

```


