#!/bin/bash
#SBATCH --job-name=processvcf
#SBATCH --time=0-3:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=/data/oselmo/CropAdapt/slurmOutput/R-%x.%j.out
#SBATCH --wait


### NOTE BEFORE RUNNING!!! -> problem in two raw input files : OrRuf_reprocessed_chr9.hmp / OrRuf_reprocessed_chr10.hmp
### --> need to manually two correct the "pos" argument for two positions ORRUF10_17000000G:A and ORRUF09_12000000G:T (reported in scientific notation -> 1.2e10+7 , need convert to number)

## convert hmp files to vcf and zip (the tr command preliminarily removes quotes from the hmp.txt file)
for hmp in DATA/RAW/rice_wr/OrRuf_reprocessed_chr*.hmp.txt; do
    
         run_pipeline.pl -h <(tr -d '"' < $hmp)  -export $hmp.vcf -exportType VCF 

         # compress vcf and index
          bgzip -c $hmp.vcf > $hmp.vcf.gz
          bcftools index $hmp.vcf.gz
          rm $hmp.vcf # remove vcf
done


# ### concatenate vcf files 
bcftools concat -Oz -o DATA/RAW/rice_wr/OrRuf_reprocessed.hmp.txt.vcf.gz $(ls DATA/RAW/rice_wr/OrRuf_reprocessed_chr*.hmp.txt.vcf.gz | sort -V)

# ### Index concatenated vcf files 
bcftools index DATA/RAW/rice_wr/OrRuf_reprocessed.hmp.txt.vcf.gz

# ### Subset vcf file,  save as vcf.gz (no call filtering since low coverage, no heterozygote calls! )
bcftools view -i 'MAF>=0.015&F_MISSING<=0.5' -R DATA/PROCESSED_VCF/rice_wr/GENES.tsv DATA/RAW/rice_wr/OrRuf_reprocessed.hmp.txt.vcf.gz -Oz -o DATA/PROCESSED_VCF/rice_wr/processed.vcf.gz

# remove int files
rm DATA/RAW/rice_wr/OrRuf_reprocessed*.hmp.txt.vcf.gz*

