#!/bin/bash
#SBATCH --job-name=prunceVCFs
#SBATCH --time=0-20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

## for every processed vcf, prune
# for vcf in DATA/PROCESSED_VCF/*/processed.vcf.gz; do

#     bcftools +prune -m 0.9 -w 1000 $vcf -Oz -o $vcf"_pruned.vcf.gz"

# done

bcftools +prune -m 0.9 -w 1000 DATA/PROCESSED_VCF/rice_la/processed.vcf.gz -Oz -o DATA/PROCESSED_VCF/rice_laoybean/processed.vcf.gz_pruned.vcf.gz
