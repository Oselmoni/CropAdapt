#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait

### Compress vcf
bgzip -c DATA/RAW/corn_la/map_152k.vcf > DATA/RAW/corn_la/map_152k.vcf.gz

### Index concatenated vcf files 
bcftools index --threads 2 DATA/RAW/corn_la/map_152k.vcf.gz

### Subset vcf file, MAF filter 1% save as vcf.gz
bcftools view -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/corn_la/GENES.tsv DATA/RAW/corn_la/map_152k.vcf.gz -Oz -o DATA/PROCESSED_VCF/corn_la/processed.vcf.gz

