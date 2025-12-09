#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait

### convert hmp files to vcf and zip (the tr command preliminarily removes quotes from the hmp.txt file)
 for hmp in DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter.c*.imp.hmp.txt; do
    
#         run_pipeline.pl -h <(tr -d '"' < $hmp)  -export $hmp.vcf -exportType VCF 

#         # compress vcf and index
          bgzip -c $hmp.vcf > $hmp.vcf.gz
          bcftools index $hmp.vcf.gz
          rm $hmp.vcf # remove vcf




   
 done


# ### concatenate vcf files 
bcftools concat -Oz -o DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter.imp.hmp.txt.vcf.gz $(ls DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter.c*.imp.hmp.txt.vcf.gz | sort -V)

# ### Index concatenated vcf files 
bcftools index DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter.imp.hmp.txt.vcf.gz

# ### Subset vcf file, only keep biallelic snps, excluded uncalled genotypes , MAF>5%, save as vcf.gz
bcftools view -m2 -M2 -v snps -U -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/sorghum_la/GENES.tsv DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter.imp.hmp.txt.vcf.gz -Oz -o DATA/PROCESSED_VCF/sorghum_la/processed.vcf.gz

## remove temp files
rm DATA/RAW/sorghum_la/sb_snpsDryad_sept2013_filter*.imp.hmp.txt.vcf.gz*