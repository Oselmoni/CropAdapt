#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

### Zip VCF file
bgzip -c DATA/RAW/teparybean/TASSEL_all.vcf > DATA/RAW/teparybean/TASSEL_all.vcf.gz

### indefx VCF file
bcftools index DATA/RAW/teparybean/TASSEL_all.vcf.gz

### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>5%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/teparybean/GENES.tsv DATA/RAW/teparybean/TASSEL_all.vcf.gz -Oz -o DATA/PROCESSED_VCF/teparybean/processed.vcf.gz

## remove intermediary files
rm DATA/RAW/teparybean/TASSEL_all.vcf.gz*