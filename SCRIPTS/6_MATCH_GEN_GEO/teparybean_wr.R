library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/teparybean////processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/teparybean_wr.csv')

## adjust sample ids
vcf$samples = unlist(lapply(strsplit(vcf$samples, ':'), function(x) {return(x[1])}))
vcf$samples[which(vcf$samples=='G40059_Menudo_Blanco_Arroz')[1]] = NA # duplicate sample, keep one only

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
