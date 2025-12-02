library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/soybean///processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/soybean_la.csv')

# clean identifier duplicate
meta$ID[meta$ID%in%names(which(table(meta$ID)==2))][rep(c(T,F), times=12)] = NA ## set duplicate samples to NA
meta$ID[meta$ID%in%names(which(table(meta$ID)==3))][1:2] = NA ## set duplicate samples to NA
meta = meta[is.na(meta$ID)==F,]

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


