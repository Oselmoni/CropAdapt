#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait



# ### convert plink format to vcf
#plink --file DATA/RAW/rice_la/NB-core_v4 --recode vcf --out DATA/RAW/rice_la/NB-core_v4.vcf

### Compress vcf
#bgzip -c DATA/RAW/rice_la/NB-core_v4.vcf.vcf > DATA/RAW/rice_la/NB-core_v4.vcf.gz

### Index vcf
#bcftools index --threads 6 DATA/RAW/rice_la/NB-core_v4.vcf.gz

# ### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>1%, save as vcf.gz
bcftools view --threads 6 -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/rice_la/GENES.tsv DATA/RAW/rice_la/NB-core_v4.vcf.gz -Oz -o DATA/PROCESSED_VCF/rice_la/processed.vcf.gz

## remove int files
rm DATA/RAW/rice_la/NB-core_v4.vcf*