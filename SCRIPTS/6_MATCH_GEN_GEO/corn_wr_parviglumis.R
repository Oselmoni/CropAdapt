library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/corn_wr/processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/corn_wr.csv')


## adjust sample ids
meta$ID = gsub('(.*\\d).*','\\1',meta$ID)



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


## load info on subspecies
info = openxlsx::read.xlsx('DATA/META_RAW/corn_wr/41588_2022_1184_MOESM3_ESM(1).xlsx', startRow = 2, sheet = 2)
subsp = info$Taxonomya
names(subsp) = info$ID
subsp = subsp[rownames(meta)]

## subset to species of interest
meta = meta[names(which(subsp=='Zea mays subsp. parviglumis')),]
GT = GT[,names(which(subsp=='Zea mays subsp. parviglumis'))]

