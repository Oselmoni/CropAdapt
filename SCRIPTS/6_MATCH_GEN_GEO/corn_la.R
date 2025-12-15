library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/corn_la///processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/corn_la.csv') # note: one coordinate per site -> to be transformed in one coordinate per sample

### adjust sample ids
vcf$samples[vcf$samples=='ZEA_3307_10_rep']  = NA # get rid of replicate sample
vcf$samples = unlist(lapply(strsplit(vcf$samples, '_'), function(x) { paste0(x[1],'_',as.numeric(x[2]),'_',x[3])}))
sites_id = unlist(lapply(strsplit(vcf$samples, '_'), function(x) { paste0(x[1],'_',as.numeric(x[2]))}))

meta$ID = unlist(lapply(strsplit(meta$ID, ' '), function(x) {paste0(x[1],'_',x[2])}))
rownames(meta) = meta$ID
meta = meta[sites_id,]
meta$ID = vcf$samples
meta = meta[is.na(meta$LON)==F,]


## Check overlap samples
ol = names(which(table(c(meta$ID, vcf$samples))==2))
print(paste0(length(ol),'/',nrow(meta),' overlapping samples'))

## Extract GT matrix
GT = vcf$gt
colnames(GT) = vcf$samples
rownames(GT) = paste0(vcf$chr,':',vcf$pos)


## Subset GT matrix
GT = GT[,ol]

## Subset meta 
rownames(meta) = meta$ID
meta = meta[ol,]


