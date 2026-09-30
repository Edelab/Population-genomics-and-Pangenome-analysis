## ClubGenoEvo Project

# Population genomics and Pangenome analysis

## Global population genomics and pangenome of the clubroot pathogen link clonal expansion and effector diversification to resistance breakdown in canola

A reproducible collection of workflows for analysing short- and long-read sequencing
datasets from the clubroot pathogen _Plasmodiophora brassicae_. The repository couples a
[Population genomics](./01_Population%20genomics) pipeline (read pre-processing → variant calling → population
structure, diversity, and selection) with a [Pangenome analysis](./02_Pangenome%20analysis) pipeline (long-read assembly →
pangenome graph construction → effector characterisation), together with downstream post-processing
of SNP and structural-variant (SV) call sets.

---

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Citation](#citation)
- [License](#license)

---

## Overview

This repository documents the full analytical workflow behind the study above. It analyses
whole-genome sequencing (WGS) data for the obligate biotrophic protist _P. brassicae_,
combining genome-wide variant analysis across a global panel of isolates with a pangenome
built from long-read assemblies. The aims are to (i) resolve population structure and
clonal expansion, (ii) quantify genomic diversity and signatures of selection, and
(iii) build a pangenome to characterise effector repertoires and structural variation linked to the breakdown
of clubroot resistance.

Each workflow folder contains numbered Markdown documents that walk through the steps in
order, alongside the scripts used to run them on local Linux and high-performance computing (HPC) environments.

## Repository structure

```
Population-genomics-and-Pangenome-analysis/
│
├── Population genomics/            # SNP-level population analyses (short reads)
│   ├── 01_WGS_pre-processing_and_variant_calling_analysis.md 
│   ├── 02_Population_structure_analysis.md
│   ├── 03_Phylogenetic_analysis.md
│   ├── 04_Genomic_summary_statistics.md
│   ├── 05_Linkage_disequilibrium_decay_analysis.md
│   ├── 06_SNP_Annotations.md
│   ├── Other_file_handling_array_jobs_and_selection_related_scripts 
│   └── README.md
│
├── Pangenome analysis/             # Long-read assembly + graph pangenome
│   ├── 01_Long_read_WGS_preprocessing_assembly_and_annotation.md
│   ├── 02_Pangenome_building_with_Minigraph-cactus.md
│   ├── 03_Gene_based_pangenome_workflow.md
│   ├── 04_Phylogenomics_single_copy_orthologue.md
│   ├── 05_SignalP6_guide.md
│   ├── 06_DeepTMHMM_guide.md
│   ├── Other_file_handling_and_array_job_scripts
│   └── README.md
│
├── Citation
├── LICENSE
└── README.md
```


## Citation

If you use these workflows, please cite:

> Javed et al., (2026). Global population genomics and pangenome of the clubroot pathogen link clonal expansion and effector diversification to resistance breakdown in canola.
> preprint (Coming soon)


## License

This project is released under the terms of the [LICENSE](./LICENSE) file in this repository.
The code and workflows were developed in [EdeLab](https://github.com/EdeLab) and are released under the MIT License to support reproducibility and reuse.

 
