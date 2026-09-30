# SnpEff database construction and VCF annotation pipeline


## Convert the gtf file to gff3 file format

```bash
agat_convert_sp_gxf2gxf.pl --gff pb3A_gene.gtf -o pb3A.gff3

```

## Run SnpEff pipeline

```bash

# 1. Add a custom snpEff.config
echo "Pb_genome.genome : Plasmodiophora brassicae" >> /path/to/directory/snpEff.config

# 2. Set up database directory
mkdir -p /path/to/directory/Pb_genome
cp Pb3A_genomic.fna /path/to/directory/Pb_genome/sequences.fa
cp pb3A.gff3 /path/to/directory/Pb_genome/genes.gff

mkdir -p logs

# 3. Build SnpEff database
java -jar /prg/snpEff/5.2e/snpEff.jar build \
  -config /path/to/directory/snpEff.config \
  -gff3 -v Pb_genome \
  2> logs/snpeff_build.log

# 4. Annotate VCF
java -Xmx8g -jar /prg/snpEff/5.2e/snpEff.jar ann \
  -config /path/to/directory/snpEff.config \
  -v Pb_genome \
  -stats snpeff_summary.html \
  samples.vcf.gz \
  > annotated.vcf \
  2> logs/snpeff_ann.log

bgzip annotated.vcf
bcftools index --tbi annotated.vcf.gz

```
