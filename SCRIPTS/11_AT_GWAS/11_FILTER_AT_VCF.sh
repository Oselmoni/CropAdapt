#!/bin/bash
#SBATCH --job-name=AT_vcf_filtering
#SBATCH --time=0-20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=30000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait

### filter vcfs
SnpSift filter "( CHROM = '1' )" DATA/ATHALIANA/RAW/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz | SnpSift filter -n "ANN[*].EFFECT has 'intergenic_region'" > DATA/ATHALIANA/RAW/AT_chr1_eff.vcf
SnpSift filter "( CHROM = '2' )" DATA/ATHALIANA/RAW/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz | SnpSift filter -n "ANN[*].EFFECT has 'intergenic_region'" > DATA/ATHALIANA/RAW/AT_chr2_eff.vcf
SnpSift filter "( CHROM = '3' )" DATA/ATHALIANA/RAW/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz | SnpSift filter -n "ANN[*].EFFECT has 'intergenic_region'" > DATA/ATHALIANA/RAW/AT_chr3_eff.vcf
SnpSift filter "( CHROM = '4' )" DATA/ATHALIANA/RAW/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz | SnpSift filter -n "ANN[*].EFFECT has 'intergenic_region'" > DATA/ATHALIANA/RAW/AT_chr4_eff.vcf
SnpSift filter "( CHROM = '5' )" DATA/ATHALIANA/RAW/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz | SnpSift filter -n "ANN[*].EFFECT has 'intergenic_region'" > DATA/ATHALIANA/RAW/AT_chr5_eff.vcf


### gunzip output
gunzip DATA/ATHALIANA/RAW/AT_chr1_eff.vcf.gz
gunzip DATA/ATHALIANA/RAW/AT_chr2_eff.vcf.gz
gunzip DATA/ATHALIANA/RAW/AT_chr3_eff.vcf.gz
gunzip DATA/ATHALIANA/RAW/AT_chr4_eff.vcf.gz
gunzip DATA/ATHALIANA/RAW/AT_chr5_eff.vcf.gz

### bgzip output
bgzip DATA/ATHALIANA/RAW/AT_chr1_eff.vcf
bgzip DATA/ATHALIANA/RAW/AT_chr2_eff.vcf
bgzip DATA/ATHALIANA/RAW/AT_chr3_eff.vcf
bgzip DATA/ATHALIANA/RAW/AT_chr4_eff.vcf
bgzip DATA/ATHALIANA/RAW/AT_chr5_eff.vcf

### index output
bcftools index DATA/ATHALIANA/RAW/AT_chr1_eff.vcf.gz
bcftools index DATA/ATHALIANA/RAW/AT_chr2_eff.vcf.gz
bcftools index DATA/ATHALIANA/RAW/AT_chr3_eff.vcf.gz
bcftools index DATA/ATHALIANA/RAW/AT_chr4_eff.vcf.gz
bcftools index DATA/ATHALIANA/RAW/AT_chr5_eff.vcf.gz