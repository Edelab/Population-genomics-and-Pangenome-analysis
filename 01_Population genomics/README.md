# Whole-genome sequence (WGS) variant analysis of _Plasmodiophora brassicae_

This workflow outlines the steps for performing quality control, read processing, mapping, variant calling, post-variant analysis, and selection analysis for whole-genome sequencing (WGS) data of _Plasmodiophora brassicae_. Each step is described with the necessary commands and their functions to ensure reproducibility.

---

## Prerequisites

- Infrastructure Environments and programming language
  - Conda
  - HPC
  - Bash
  - Python
  - R
  - C++  

- Software and tools required:
  - FastQC
  - Trimmomatic
  - BWA-mem2
  - SAMtools
  - Sambamba
  - BCFtools
  - VCFtools
  - PLINK
  - ADMIXTURE
  - MAFFT
  - IQ-TREE
  - SplitsTree
  - SnpEff
  - Pixy
    
- Input files:
  - Paired-end Illumina reads (`Read_1.fq` and `Read_2.fq`)
  - Reference genome (`genome#.fa`)

---
