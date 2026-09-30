# Whole-Genome Seqeuencing Analysis(WGS) of _Plasmodiophora brassicae_ Genome 

1. Quality checking of Illumina data (Paired end) using FastQC

```bash
fastqc path_to_directory/Read_1.fq path_to_directory/Read_2.fq  -o path_to_output_directory
```
2. Adapter removal and unwanted/low qaulity reads trimming using Trimmomatic [Default Parameters]

```bash
java -jar /prg/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 20 path_to_directory/read_1.fq path_to_directory/read_2.fq sample#_forward_paired.fq sample#_forward_unpaired.fq sample#_reverse_paired.fq sample#_reverse_unpaired.fq ILLUMINACLIP:/path_to_adapter_directory/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 
```
3. Indexing and Mapping of trimmed reads with reference genome with BWA_mem2 with adding read groups

```bash
bwa-mem2 index /path_to_reference_genome/genom#.fa
bwa-mem2 mem -R "@RG\tID:sample_ID\tSM:sample_name\tPL:illumina" -t 15 /path_to_indexed_reference_genome/genome#.fa /path_to_trimmed_reads/output_sample#_forward_paired.fq /path_to_trimmed_reads/output_reverse_paired.fq > /path_to_output_mapping/sample#_mapping.sam

# [-R = read group contains information such as the sample name, ID and platform].
```

3. Convert sam to bam with samtools, sorting and marking duplicates using Sambamba 

```bash
samtools view -Sbh sample#.sam > sample#.bam

# Sorting and duplicates marking by Sambamba

sambamba sort -o sample#_mapped_sorted.bam sample#_mapped.bam

# MarkDuplicates by by Sambamba

sambamba markdup sample#_mapped_sorted.bam sample#_mapped_sorted.dd.bam

## sambamba generates all index files automatically from sorted and duplicates files

```
 
4. Checking Alignment Quality and average depth of coverage using samtools

```bash
samtools stats sample#_sorted_markdup.bam > sample#.stats.txt

samtools flagstat sample#_sorted_markdup.bam > sample#_qc.txt

samtools depth sample#_sorted_markdup.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'     

```

5. Pre-processing before using bcftools  mpileup

```bash
# [ It requires all sorted_markdup.bam files into index bam.bai, create .dict and indexed fa.fai file from reference fasta file with same name]

samtools faidx referenc_genome_pb3A.fa

bcftools mpileup --threads 30 -a AD,DP,SP,ADF,ADR -f /home/edelab/asim_work/reference_genome_pb3A/Pb3A_genomic.fna /data1/Asim_backup_files/Asim_Data/pb_international_batch/dk21975/dk21975_mapped_sorted.dd.bam -q 5 -d 20000 |bcftools call --threads 30 -m -G- --ploidy 1 -a GQ,GP -Oz -o dk21975_vcf.gz


# For multiple samples

bcftools mpileup --threads 30 -a AD,DP,SP,ADF,ADR -f /home/edelab/asim_work/reference_genome_pb3A/Pb3A_genomic.fna -q 5 -d 20000 -b final_bam_files.txt | bcftools call --threads 30 -m -G- --ploidy 1 -a GQ,GP -Oz -o pb_2025_06_vcf.gz


```

6. Post-variant analysis for parameters selection before variant filtration using bcftools

```bash
# counting SNPs
bcftools view -H variant_sample#.raw.snps.vcf.gz | wc -l


# Variant quality scan for filtration

bcftools view variant_sample#.raw.snps.vcf.gz | vcfrandomsample -r 0.012 > sample_subset.vcf         

```
7. Variant SNPs/INDELs filtering using multiple parameters

7.1. Filter minor allele count (MAC)

```bash

vcf=$1
bcftools filter -i 'MAC >1' -Oz -o 04_filtered_vcfs/MAC_1_biallelic_NOmissing_10_all_samples_filterd_cov0.15_3.0_0.1_0.9_40.0.vcf.gz

```
7.2. Filter missingness

```bash
#!/bin/bash

# modules
module load bcftools/1.15

# variables
vcf=$1
bad_samples=$2
max_missing=$3

# constant
out="04_filtered_vcfs/NOmissing_$(basename $vcf)"

# run

bcftools view --force-samples -S ^$bad_samples $vcf -q 0.00000001:minor | bcftools filter -i "F_MISSING < $max_missing" -Oz -o $out
bcftools index $out

echo "done"
```
7.3. Filter vcf file based on coverage, genotype quality and allelic ratio of coverage

```bash

#!/usr/bin/env python
#Python script to filter vcf file based on coverage, genotype quality and allelic ratio of coverage
#Fitlering is done on an individual bases
#Coverage threshold is applied in % of the median for lower and upper bound
#Genotype not passing the filter are assigned missing value
#usage: filtering_vcf.py --vcf <file> --median <file> --lcov <float> --highcov <float> --lratio <float> -- hratio <float> --gq <float>

import gzip as gz
import argparse as arg
from collections import defaultdict

## Define class
class vcf_parser_line():
    sample_list = []
    #filtering_status = defaultdict(lambda: defaultdict(int))

    def __init__(self, vcf_line: str):
        self.line = vcf_line.strip().split("\t")

        if vcf_line.startswith("#CHROM"):
            vcf_parser_line.sample_list = self.line[9:]

        if vcf_line.startswith("#"):
            self.header = True
            return

        self.header = False
        self.linestart = self.line[0:8]
        self.info = self.line[8].split(":")
        self.info_dict = defaultdict(list)
       
        for ind in self.line[9:]:
             for pos, field in enumerate(ind.split(":")):
                 self.info_dict[self.info[pos]].append(field)


    def rebuild_line(self):
        geno = "\t".join([":".join(x) for x in zip(*self.info_dict.values())])
        return "\t".join(self.linestart) + "\t" + ":".join(self.info) + "\t" + geno

    def filter_lcov(self, median_dict, thresh):
        if "DP" not in self.info:
            print("cannot filter on Coverage")
            return

        for ind in range(len(vcf_parser_line.sample_list)):
            if float(self.info_dict["DP"][ind]) < thresh * median[vcf_parser_line.sample_list[ind]]:
                self.info_dict["GT"][ind] = "."
                #vcf_parser_line.filtering_status[vcf_parser_line.sample_list[ind]]["lowCov"]+=1
    
    def filter_hcov(self, median_dict, thresh):
        if "DP" not in self.info:
            print("cannot filter on Coverage")
            return
        
        for ind in range(len(vcf_parser_line.sample_list)):
            if float(self.info_dict["DP"][ind]) > (thresh * median[vcf_parser_line.sample_list[ind]]):
                self.info_dict["GT"][ind] = "."
                #vcf_parser_line.filtering_status[vcf_parser_line.sample_list[ind]]["highCov"]+=1

    def filter_ratio(self, thresh):
        if "AD" not in self.info:
            print("cannot filter on ratio")
            return
    
        for ind in range(len(vcf_parser_line.sample_list)):
            AD = [int(x) for x in self.info_dict["AD"][ind].split(",")]
            AD = AD[0] / sum(AD) if sum(AD) >0 else 0
           
            if (AD <= thresh[1] and AD >= thresh[0]):
                self.info_dict["GT"][ind] = "."
                #vcf_parser_line.filtering_status[vcf_parser_line.sample_list[ind]]["allelicRatio"]+=1

    def filter_GQ(self, thresh):
        if "GQ" not in self.info:
            print("cannot filter on ratio")
            return

        for ind in range(len(vcf_parser_line.sample_list)):
            GQ = float(self.info_dict["GQ"][ind])
            if GQ < thresh:
                self.info_dict["GT"][ind] = "."
                #vcf_parser_line.filtering_status[vcf_parser_line.sample_list[ind]]["GQ"]+=1

#    @staticmethod
#    def parse_dict():
#        line = ""
#        for ind in vcf_parser_line.filtering_status:
#            subdict = vcf_parser_line.filtering_status[ind]
#            line += "\n".join([f"{ind}\t{metric}\t{subdict[metric]}" for metric in subdict])
#            line += "\n"
#        return line

def parse_median_table(file):
    with open(file, "r") as f:
       a = [z.strip().split("\t") for z in f.readlines()]
       return{x: int(y) for x, y, *_, in a}

def myopen(file, mode = "rt"):
    if file.endswith(".gz"):
        return gz.open(file, mode)
    return open(file, mode)


## Define argument parser
parser = arg.ArgumentParser()
parser.add_argument("--vcf", type = str, required = True)
parser.add_argument("--median", type = str, required = True)
parser.add_argument("--lcov", default = 0.1, type = float)
parser.add_argument("--hcov", default = 3, type = float)
parser.add_argument("--lratio", default = 0, type = float)
parser.add_argument("--hratio", default = 1, type = float)
parser.add_argument("--gq", default = 50, type = float)
args = parser.parse_args()
median = parse_median_table(args.median)
out = "04_filtered_vcfs/all_samples_filterd_cov" + str(args.lcov) + "_" + str(args.hcov) + "_" + str(args.lratio)+ "_" + str(args.hratio) + "_" + str(args.gq) + ".vcf.gz"
indat = myopen(args.vcf)
output = myopen(out, "wt")

for line in indat:
    l = vcf_parser_line(line)
    if l.header:
        output.write(line)
        continue

    if len(l.linestart[3]) >1 or len(l.linestart[4]) >1:
        continue

    l.filter_lcov(median, args.lcov)
    l.filter_hcov(median, args.hcov)
    if l.linestart[4] != "." : #check if polymorpic, then filter genotype
        l.filter_ratio([args.lratio, args.hratio])
        l.filter_GQ(args.gq)
    output.write(l.rebuild_line() + "\n")

```

8. Extract filtering results

 ```bash

#!/bin/bash

#variables
vcf=$1

#script
filecode=$(echo $(basename $vcf) | sed "s/filterd_cov//" | sed "s/\.vcf\.gz//")
vcftools --gzvcf $vcf --out 98_sumstats/$filecode --missing-indv
vcftools --gzvcf $vcf --out 98_sumstats/$filecode --missing-site

```

9. biallelic, multiallelic quanitification and filteration

```bash

## SNPs profiling using bcftools and counting the biallelic and multiallelic SNPs. 

# Total SNPs count in good genomes (18) = 205643

bcftools view -m2 -M2 -v snps global_pb_SNPs.vcf.gz | bcftools stats | grep 'number of SNPs:'   # Here used the good genome dataset (18 genomes) and counted the number of	201699 biallelic SNPs.

bcftools view -m3 -v snps pb_filtered_snps.vcf | bcftools stats | grep 'number of SNPs:'  # Here used the good genome dataset (18 genomes) and counted the number of 992 multiallelic SNPs.

bcftools view -m2 -v snps pb_filtered_snps.vcf | bcftools stats | grep 'number of SNPs:'  # command to count the total number of biallelic & multiallelic SNPs.

# To filter out the only specific SNPs then use alongwith the output command

bcftools view -m2 -M2 -v snps pb_filtered_snps.vcf -Oz -o pb_filtered_bialleic_snps.vcf.gz

```

10. To sort the filtered vcf.gz file

```bash
bcftools sort your_file.vcf -Oz -o sorted_file.vcf.gz

# To index the final vcf.gz file
tabix -p vcf your_file.vcf.gz

```

11. Script that extracts core quality variables per sample

```bash

#!/bin/bash

##Script that extract core quality variables per samples:
# Extract coverage, allelic ratio and genotype quality
#Usage: extract_sumstats.sh <path_to_vcf>

vcf=$1


echo -e "SAMPLE\tCOV\tALL1\tALL2\tGT\tQUAL" >98_sumstats/sumstat_quality.txt

bcftools query -f "[%SAMPLE\t%DP\t%AD{0}\t%AD{1}\t%GT\t%GQ\n]" $vcf >>98_sumstats/sumstat_quality.txt
echo done

```


