# Pangenome analysis of — _Plasmodiophora brassicae_

This directory contains the full bioinformatics pipeline used to build and 
analyse the _P. brassicae_ pangenome and downstream functional annotation across 55 isolates, from raw long-read 
assembly to effector and transmembrane protein prediction.


### Pipeline guides


| `01_Long_read_WGS_preprocessing_assembly_and_annotation.md` | ONT read QC, de novo assembly with Hifiasm and Canu, contamination removal, repeat masking, and gene prediction with BRAKER3 |

| `02_Pangenome_building_with_Minigraph-cactus.md` | Pangenome graph construction with Minigraph-Cactus, SV detection, vcfbub and vcfwave processing, and SVTYPE annotation |

| `03_Gene_based_pangenome_workflow.md` | Gene-based pangenome analysis using OrthoFinder

| `04_Phylogenomics_single_copy_orthologue.md` | Phylogenetic tree construction from single-copy orthologue alignments and IQ-TREE2 with _Spongospora subterranea_ outgroup |

| `05_SignalP6_guide.md` | Signal peptide prediction using SignalP 6.0 slow-sequential mode via SLURM array on 55 isolates |

| `06_DeepTHmm_guide.md` | Transmembrane topology prediction using DeepTMHMM on SignalP-processed proteins via GPU SLURM array |


## Prerequisites

The main tools used across this pipeline are listed below. Full installation 
notes are in the relevant guide for each tool.

- Hifiasm, Canu
- NCBI FCS-GX 
- RepeatModeler2, RepeatMasker, EarlGrey
- BRAKER3 
- Minigraph-Cactus
- OrthoFinder
- MAFFT, IQ-TREE2
- SignalP 6.0
- DeepTMHMM
