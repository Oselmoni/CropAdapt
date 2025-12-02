library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/rice_wr///processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/rice_wr.csv')


## Extract GT matrix
GT = vcf$gt
colnames(GT) = vcf$samples
rownames(GT) = paste0(vcf$chr,':',vcf$pos)


## Check overlap samples
ol = names(which(table(c(meta$ID, vcf$samples))==2))
print(paste0(length(ol),'/',nrow(meta),' overlapping samples'))


## Subset GT matrix
GT = GT[,ol]

## Subset meta 
rownames(meta) = meta$ID
meta = meta[ol,]
