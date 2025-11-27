#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait


### concatenate vcf files 
bcftools concat -Oz -o DATA/RAW/soybean/AnLab_1.5K.vcf_.gz $(ls DATA/RAW/soybean/Chr*.AnLab_1.5K.vcf_.gz | sort -V)

### Index concatenated vcf files 
bcftools index DATA/RAW/soybean/AnLab_1.5K.vcf_.gz

### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>5%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/soybean/GENES.tsv DATA/RAW/soybean/AnLab_1.5K.vcf_.gz -Oz -o DATA/PROCESSED_VCF/soybean/processed.vcf.gz

## remove temporary files
rm DATA/RAW/soybean/AnLab_1.5K.vcf_.gz