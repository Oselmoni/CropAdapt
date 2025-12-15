library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/tomato_wr/processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/tomato_wr.csv')

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


### remove duplicate SNPs

## Check which SNPs have more than one line
duplicateSNPS = names(which(table(rownames(GT))>1))

for (snp in duplicateSNPS) {
  
  # find index of duplicates
  idx = which(rownames(GT)==snp)
  
  # find which one has less missing values
  to_keep = which.min(apply(GT[idx,], 1, function(x) {sum(is.na(x))}))
  
  # get idx of duplicates to remove
  idx_to_remove = idx[1:length(idx)!=to_keep]
  
  # remove duplicate
  GT=GT[-idx_to_remove,]
}

