#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait


### Index vcf file 
bcftools index --threads 4 DATA/RAW/sunflower_wr/Annuus.ann_env.tranche90_snps_bi_AN50_AF99.vcf.gz

### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>5%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' --threads 4 -R DATA/PROCESSED_VCF/sunflower_wr/GENES.tsv DATA/RAW/sunflower_wr/Annuus.ann_env.tranche90_snps_bi_AN50_AF99.vcf.gz -Oz -o DATA/PROCESSED_VCF/sunflower_wr/processed.vcf.gz