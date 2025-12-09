#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait

### Index bcf file 
bcftools index DATA/RAW/barley/variomeSNPs_final815samples_filteredSNPs.bcf 

### Subset bcf file and convert to vcf, only keep biallelic snps, excluded uncalled genotypes , MAF>1%, then use annotate to exclude annotations on snp effect
bcftools view --threads 6  -m2 -M2 -v snps -U -i 'MAF>=0.01&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/barley/GENES.tsv DATA/RAW/barley/variomeSNPs_final815samples_filteredSNPs.bcf  | bcftools annotate --threads 6 -x INFO/EFF -Oz -o DATA/PROCESSED_VCF/barley/processed2.vcf.gz
