#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait

# ### convert hmp file to vcf 
run_pipeline.pl -h DATA/RAW/commonbean_wr/GBS_WildBean.hmp.txt  -export DATA/RAW/commonbean_wr/GBS_WildBean.vcf -exportType VCF 

# ### Compress vcf
bgzip -c DATA/RAW/commonbean_wr/GBS_WildBean.vcf > DATA/RAW/commonbean_wr/GBS_WildBean.vcf.gz

# ### Index vcf
bcftools index DATA/RAW/commonbean_wr/GBS_WildBean.vcf.gz

# ### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>1%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.01&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/commonbean_wr/GENES.tsv DATA/RAW/commonbean_wr/GBS_WildBean.vcf.gz -Oz -o DATA/PROCESSED_VCF/commonbean_wr/processed.vcf.gz

### remove intermediary files 
rm DATA/RAW/commonbean_wr/GBS_WildBean.vcf*
