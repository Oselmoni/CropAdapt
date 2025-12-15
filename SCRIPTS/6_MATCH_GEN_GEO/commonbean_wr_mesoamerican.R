library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/commonbean_wr//processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/commonbean_wr.csv')

# reformat samples id
vcf$samples = unlist(lapply(strsplit(vcf$samples, ':'), function(x) {return(x[1])}))
meta$ID = gsub('G (.*)', 'G\\1',meta$ID)


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
meta = meta[names(which(pcoa[,1]>(0))),]
GT=GT[,names(which(pcoa[,1]>(0)))]


