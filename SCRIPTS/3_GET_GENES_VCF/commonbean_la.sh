#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

### gunzip vcf files
for vcf in DATA/RAW/commonbean_la/628sample.Chr*.recode.vcf; do
    bgzip -c $vcf > ${vcf}.gz
done

### index vcf files
for vcf in DATA/RAW/commonbean_la/628sample.Chr*.recode.vcf.gz; do
    bcftools index $vcf
done

### concatenate vcf files 
bcftools concat -Oz -o DATA/RAW/commonbean_la/628sample.recode.vcf.gz $(ls DATA/RAW/commonbean_la/628sample.Chr*.recode.vcf.gz | sort -V)

### Index concatenated vcf files 
bcftools index DATA/RAW/commonbean_la/628sample.recode.vcf.gz

### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>1%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.01&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/commonbean_la/GENES.tsv DATA/RAW/commonbean_la/628sample.recode.vcf.gz -Oz -o DATA/PROCESSED_VCF/commonbean_la/processed.vcf.gz

### remove intermediary files 
rm DATA/RAW/commonbean_la/628sample.Chr*.recode.vcf.gz*
rm DATA/RAW/commonbean_la/628sample.recode.vcf.gz*