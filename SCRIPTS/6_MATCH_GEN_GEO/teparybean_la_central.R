library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/teparybean////processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/teparybean_la.csv')

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


## Run preliminary PCoA
set.seed(0);pcoa = cmdscale(dist(t(GT[sample(1:nrow(GT), 1000),])))


## Keep only variety of interest
meta = meta[names(which(pcoa[,1]<(-10))),]
GT=GT[,names(which(pcoa[,1]<(-10)))]



## Run preliminary PCoA #2
set.seed(0);pcoa = cmdscale(dist(t(GT[sample(1:nrow(GT), 1000),])))
plot(pcoa)

## Keep only variety of interest
meta = meta[names(which(pcoa[,1]<(3))),]
GT=GT[,names(which(pcoa[,1]<(3)))]


