#!/bin/bash
#SBATCH --job-name=prunceVCFs
#SBATCH --time=0-10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

## for every processed vcf, prune
for vcf in DATA/PROCESSED_VCF/*/processed.vcf.gz; do

    bcftools +prune -m 0.9 -w 1000 $vcf -Oz -o $vcf"_pruned.vcf.gz"

done
