library(vcfppR)

## Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/rice_la/processed.vcf.gz_pruned.vcf.gz')

## Load meta
meta = read.csv('DATA/META/rice_la.csv')

## adjust samples ids
vcf$samples = unlist(lapply(strsplit(vcf$samples, '_'), function(x) { paste0(x[1],'_',x[2])}))

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


### find metadata on subspecies
info = openxlsx::read.xlsx('DATA/META_RAW//rice_la/41477_2020_659_MOESM4_ESM(1).xlsx', startRow = 16, sheet = 1)

# keep only 3k samples
info = info[info$SEQ=='Wang_et_al_2018',]

# find susbsp info
subsp = info$SSP
names(subsp) = info$ID
subsp = subsp[rownames(meta)]


### keep only samples of subsp of interest
GT = GT[,subsp=='JAP']
meta = meta[subsp=='JAP',]

### # some japonica samples without the LINGO annotation in metadata appear as genetic outliers, remove them
meta = meta[-which(meta$ID%in%info$ID[which(is.na(info$LINGO))]),]
GT=GT[,meta$ID]

