### Load gtf file
GTF = read.table('DATA/RAW/sorghum_la/Sbi1.4.gff3', sep='\t')


### keep only genes
GENES = GTF[GTF$V3=='mRNA',c(1,4,5)]



### Keep only genes on chromsomes
GENES = GENES[substr(GENES$V1,1,3)=='chr',]


### convert chromosome id
GENES$V1 = substr(GENES$V1,12,nchar(GENES$V1))


## remove genes outside main chromosomes 
GENES=na.omit(GENES)


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/sorghum_la/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)

