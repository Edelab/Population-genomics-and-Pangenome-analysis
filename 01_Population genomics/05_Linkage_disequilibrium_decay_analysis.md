# Linkage disequilibrium (LD) decay analysis workflow

This repository contains the complete workflow for linkage disequilibrium (LD) decay analysis used to assess reproductive modes and population structure in *P. brassicae*. The pipeline uses PLINK for LD computation with native haploid support.


### Key Features

- Native haploid support via PLINK
- Nonlinear least squares (NLS) exponential decay modeling
- Median-based binning for robust statistics


### 1. Population Clustering

Isolates were grouped into clusters based on PCA-based population structure analysis. Create separate VCF files for each cluster:

```bash
# Extract samples for each cluster
vcftools --vcf all_samples.vcf \
  --keep cluster1_samples.txt \
  --recode --out cluster1

vcftools --vcf all_samples.vcf \
  --keep cluster2_samples.txt \
  --recode --out cluster2

vcftools --vcf all_samples.vcf \
  --keep cluster3_samples.txt \
  --recode --out cluster3
```

### 2. Cluster 3 Subpopulation Analysis

For fine-scale analysis within clusters:

```bash
# Split Cluster 3 into subpopulations A and B
vcftools --vcf cluster3.recode.vcf \
  --keep cluster3_popA_samples.txt \
  --recode --out cluster3_popA

vcftools --vcf cluster3.recode.vcf \
  --keep cluster3_popB_samples.txt \
  --recode --out cluster3_popB
```

---

## 3. LD Computation with PLINK

### Standard Workflow for Each Cluster 

```bash

# Step 1: Convert VCF to binary format with haploid support
plink --vcf cluster1_ld_ready.vcf.gz \
  --allow-extra-chr \
  --chr-set -20 \
  --double-id \
  --make-bed \
  --out ./cluster1/c1_haploid

# Step 2: Compute pairwise LD (PLINK v1.90b5.3 64-bit)
 plink --bfile ./cluster1/c1_haploid 
 --allow-extra-chr
  --chr-set -20
  --ld-window 999999
  --ld-window-kb 100
  --ld-window-r2 0
  --out ./cluster1/c1_100kb_ld
  --r2        

```

### Key PLINK Parameters Explained

| Parameter | Purpose |
|-----------|---------|
| `--chr-set -20` | Enable haploid mode for 20 chromosomes |
| `--ld-window-kb 100` | Restrict analysis to 100 kb windows |
| `--ld-window 999999` | Remove SNP-count limitations |
| `--ld-window-r2 0` | Retain all r² values including low-LD pairs |
| `--r2 dprime` | Compute both r² and D' statistics (also ran the dprime statistics as well)|

