library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/commonbean_la//processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/commonbean_la.csv')

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


### load metadata samples with info on genepools
INFO = openxlsx::read.xlsx('DATA/META_RAW//commonbean_la/41588_2019_546_MOESM3_ESM.xlsx', startRow = 2)
genepools = INFO$Genepool
names(genepools) = INFO$`ID-SEQ`

genepools = genepools[rownames(meta)]

### retain only genepool of interest
genepools = genepools[genepools=='Mesoamerican']
meta=meta[names(genepools),]
GT=GT[,rownames(meta)]

