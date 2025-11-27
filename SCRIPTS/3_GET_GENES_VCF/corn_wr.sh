#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-10:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

# ### index vcf files
# for vcf in DATA/RAW/corn_wr/merge_*.filter.vcf.gz; do
#     bcftools index $vcf
# done


### concatenate vcf files 
#bcftools concat --threads 6 -Oz -o DATA/RAW/corn_wr/merge.filter.vcf.gz $(ls DATA/RAW/corn_wr/merge_*.filter.vcf.gz | sort -V)

### Index concatenated vcf files 
#bcftools index --threads 6 DATA/RAW/corn_wr/merge.filter.vcf.gz

### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>1%, save as vcf.gz
bcftools view --threads 6 -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/corn_wr/GENES.tsv DATA/RAW/corn_wr/merge.filter.vcf.gz -Oz -o DATA/PROCESSED_VCF/corn_wr/processed.vcf.gz

## remove int
rm DATA/RAW/corn_wr/merge.filter.vcf.gz*