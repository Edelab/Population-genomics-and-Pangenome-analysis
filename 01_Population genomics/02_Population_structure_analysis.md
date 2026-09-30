## Population structure and genomics analysis of collected P. brassicae isolates from different regions/countries including provinces

1. Population structure analysis and LD pruning pipeline

```bash
#!/bin/bash


VCFIN="/path/to/directory/samples.vcf.gz"
DIROUT="09_LD"

plink --vcf $VCFIN --out $DIROUT/LDprunned_$(basename ${VCFIN/.vcf.gz}) --indep-pairwise 50 20 0.2 --double-id  --allow-extra-chr --set-missing-var-ids @:#

sed -i "s/:/\t/" $DIROUT/LDprunned_$(basename ${VCFIN/.vcf.gz}).prune.in
bcftools view -T $DIROUT/LDprunned_$(basename ${VCFIN/.vcf.gz}).prune.in -Oz -o $DIROUT/LDprunned_$(basename $VCFIN) $VCFIN

bcftools view -h $DIROUT/LDprunned_$(basename $VCFIN) >$DIROUT/NONCODING_LDprunned_$(basename ${VCFIN/.gz})

bedtools intersect -b /path/to/directory/pb3A_sorted_filtered_annotations_NCBI.gff3  -a $DIROUT/LDprunned_$(basename $VCFIN) -v >> $DIROUT/NONCODING_LDprunned_$(basename ${VCFIN/.gz})
sed "s/chromosome//" $DIROUT/NONCODING_LDprunned_$(basename ${VCFIN/.gz}) > $DIROUT/num_NONCODING_LDprunned_$(basename ${VCFIN/.gz})

```

2.  Identify and remove clonal samples from LD-pruned VCF

```bash
#!/bin/bash


vcf="/path/to/directory/141.vcf.gz"


plink --vcf $vcf --allow-extra-chr --distance ibs square --double-id --out 09_LD/individual_IBS

Rscript 01_script/00_utils/ID_clonal_samples.R 
bcftools view -S Clonal_samples -Oz -o 09_LD/indep_samples.vcf.gz $vcf
plink --vcf $vcf --out /path/to/directory/ --double-id --allow-extra-chr --chr-set -20 --recode A
plink --vcf /path/to/directory/indep_samples.vcf.gz --out /path/to/directory/indep_samples --double-id --allow-extra-chr --chr-set -20 --recode A

```

3. PCA of LD-pruned genotypes for 141 isolates and clonal-removed dataset

```R
library(magrittr)
library(dplyr)
library(ggplot2)

setwd("/path/to/directory/vcf/")

# ---- Function: impute NA to the most common genotype per SNP ----
Input_012 <- function(Gentable){
  for(i in 1:ncol(Gentable)){
    Gentable[,i][is.na(Gentable[,i])] <- as.numeric(names(which.max(table(Gentable[,i]))))
  }
  return(Gentable)
}

tidy_names <- function(x) toupper(as.character(x))

# Read data 
all_samples <- read.table("./path/to/directory/all_samples.raw", h = T) %>%
  select(-c(FID, PAT, MAT, SEX, PHENOTYPE))

indep_samples <- read.table("./indep_samples.raw", h = T) %>%
  select(-c(FID, PAT, MAT, SEX, PHENOTYPE))

# Order by isolate name so scores rows stay aligned with IID
all_samples   <- all_samples[order(all_samples$IID), ]
indep_samples <- indep_samples[order(indep_samples$IID), ]

# ---- Impute + PCA (covariance PCA, scale = F, matching your main figure) 
imputed_all   <- Input_012(all_samples[, -1])
imputed_indep <- Input_012(indep_samples[, -1])

pca_all   <- prcomp(imputed_all,   scale = F)
pca_indep <- prcomp(imputed_indep, scale = F)

# Variance explained = eigenvalue / sum(eigenvalues) = sdev^2 / sum(sdev^2) 

ve_all   <- pca_all$sdev^2   / sum(pca_all$sdev^2)   * 100
ve_indep <- pca_indep$sdev^2 / sum(pca_indep$sdev^2) * 100

# ---- Region assignment by name prefix (case-insensitive; matches your Plotly scheme) 
assign_region <- function(x){
  x <- as.character(x)
  dplyr::case_when(
    grepl("^AUS",                x, ignore.case = TRUE) ~ "Asia & Oceania",
    grepl("^(DAR|CN|JP)",        x, ignore.case = TRUE) ~ "Asia & Oceania",
    grepl("^(US|BZ|MX|CL)",      x, ignore.case = TRUE) ~ "America & South America",
    grepl("^SA01",               x, ignore.case = TRUE) ~ "Africa",
    grepl("^sa_",                x, ignore.case = TRUE) ~ "Africa",
    grepl("^AU[0-9]",            x, ignore.case = TRUE) ~ "Europe",
    grepl("^(DK|GM|EE|FI|F1|NO|FR|IT|PL)", x, ignore.case = TRUE) ~ "Europe",
    grepl("ABC",                 x, ignore.case = TRUE) ~ "Canada",
    grepl("^(17126D|QC|SK|ON|AB_|TH|4A|5X|6|8)", x, ignore.case = TRUE) ~ "Canada",
    TRUE ~ "Unassigned"
  )
}

region_levels <- c("Canada", "Europe", "America & South America",
                   "Africa", "Asia & Oceania", "Unassigned")

# Fill (interior) + darker border colours, from your Plotly script
region_fill <- c(
  "Canada"                  = "#70b4c4",
  "Europe"                  = "#f9dc86",
  "America & South America" = "#eaa0b2",
  "Africa"                  = "#a1b995",
  "Asia & Oceania"          = "#ca938c",
  "Unassigned"              = "#cccccc"
)
region_border <- c(
  "Canada"                  = "#4da3b8",
  "Europe"                  = "#fabb38",
  "America & South America" = "#d34765",
  "Africa"                  = "#497337",
  "Asia & Oceania"          = "#933b33",
  "Unassigned"              = "#666666"
)

# Build score data frames (region from original IID; Isolate tidied for display)

scores_all <- as.data.frame(pca_all$x) %>%
  mutate(Region  = factor(assign_region(all_samples$IID), levels = region_levels),
         Isolate = tidy_names(all_samples$IID))

scores_indep <- as.data.frame(pca_indep$x) %>%
  mutate(Region  = factor(assign_region(indep_samples$IID), levels = region_levels),
         Isolate = tidy_names(indep_samples$IID))

# Sanity check 
cat("\n== 141-sample region counts ==\n"); print(table(scores_all$Region))
cat("\n== 15-sample region counts ==\n");  print(table(scores_indep$Region))
unplaced <- scores_all$Isolate[scores_all$Region == "Unassigned"]
if (length(unplaced) > 0) { cat("\nUNASSIGNED (verify):\n"); print(unplaced) }

# Plot builder
#      (x = -PC1, y = PC2, no axis reversal)
make_pca_plot <- function(df, ve, pt_size, x_by, y_by){
  ggplot(df, aes(x = -PC1, y = PC2, fill = Region, colour = Region)) +
    geom_point(shape = 21, size = pt_size, stroke = 1) +
    scale_fill_manual(values = region_fill,   drop = TRUE, name = "Region") +
    scale_colour_manual(values = region_border, drop = TRUE, guide = "none") +
    scale_x_continuous(breaks = scales::breaks_width(x_by)) +
    scale_y_continuous(breaks = scales::breaks_width(y_by)) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#e5e5e5", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.ticks       = element_blank()
    )
}

p_all   <- make_pca_plot(scores_all,   ve_all,   pt_size = 7, x_by = 15, y_by =7)
p_indep <- make_pca_plot(scores_indep, ve_indep, pt_size = 7, x_by = 15, y_by = 7)


print(p_all)
print(p_indep)


```

4. ADMIXTURE analysis of 141 isolates and clonal-removed dataset across K1-K15 with 10 replicates

```bash

#!/bin/env bash

# Constants

vcf="/path/to/directory/all_samples.vcf.gz"
indep_vcf="/path/to/directory/indep_samples.vcf.gz"
OUTDIR="10_admixture"

#code
plink --vcf $vcf --out $OUTDIR/path/to/directory/141_noncoding_LD_pruned --make-bed  --double-id --allow-extra-chr
plink --vcf $indep_vcf --out $OUTDIR/path/to/directory/noncoding_indep_LD_pruned --make-bed  --double-id --allow-extra-chr
cd $OUTDIR

for i in 1 2 3 4 5 6 7 8 9 10
	do

	for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 
	do	
		echo "STARTING ADMIXTURE"
		echo $(pwd)
		admixture -s time --haploid="*" --cv ./path/to/directory/141_noncoding_LD_pruned.bed $k >>${i}_141_noncoding_LD_pruned.log
                admixture -s time --haploid="*" --cv ./path/to/directory/noncoding_indep_LD_pruned.bed $k >> ${i}_noncoding_indep_LD_pruned.log
		
		mv 141_noncoding_LD_pruned.${k}.Q all_samples/admixture/141_noncoding_LD_pruned.${k}.${i}.Q
                mv 141_noncoding_LD_pruned.${k}.P 141_samples/admixture/141_noncoding_LD_pruned.${k}.${i}.P
		
                mv noncoding_indep_LD_pruned.${k}.Q indep_samples/admixture/noncoding_indep_LD_pruned.${k}.${i}.Q
                mv noncoding_indep_LD_pruned.${k}.P indep_samples/admixture/noncoding_indep_LD_pruned.${k}.${i}.P

	done
done

mv *${i}_141_noncoding_LD_pruned.log normal/admixture/
mv *${i}_noncoding_indep_LD_pruned.log strict/admixture/

cd ..

```

5. PERMANOVA of population structure by cluster origin and host across 141 isolates

```bash

## PERMANOVA WITH ALL SAMPLES
# Permanova
Permanova.global <- adonis2(imputed_SNPs_all_samples~( Cluster* Detailed_Origin * Host  ), data = all_samples_meta, method = "euc", by = "terms", permutations = 10000)
Permanova.global

# Cluster 3
Cluster3_Metadata <- all_samples_meta[which(all_samples_meta$Cluster == 3),]
Cluster3_data <- imputed_SNPs_all_samples[which(all_samples_meta$Cluster == 3),]
permanova.Cluster3 <-adonis2(Cluster3_data~( Detailed_Origin* Host), data = Cluster3_Metadata, method = "euc", by = "terms", permutations = 10000)

# Cluster 2
Cluster2_Metadata <- all_samples_meta[which(all_samples_meta$Cluster == 2),]
Cluster2_data <- imputed_SNPs_all_samples[which(all_samples_meta$Cluster == 2),]
permanova.Cluster2 <- adonis2(Cluster2_data~( Detailed_Origin* Host), data = Cluster2_Metadata, method = "euc", by = "terms", permutations = 10000)

# Cluster 1
Cluster1_Metadata <- all_samples_meta[which(all_samples_meta$Cluster == 1),]
Cluster1_data <- imputed_SNPs_all_samples[which(all_samples_meta$Cluster == 1),]
permanova.Cluster1 <- adonis2(Cluster1_data~( Detailed_Origin* Host), data = Cluster1_Metadata, method = "euc", by = "terms", permutations = 10000)

```



