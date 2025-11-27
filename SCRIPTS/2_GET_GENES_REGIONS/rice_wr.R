### Load gtf file
GTF = read.table('DATA/RAW/rice_wr/genomic.gtf', sep='\t')


### keep only genes
GENES = GTF[GTF$V3=='transcript',c(1,4,5)]


### Create chromsome id
chrom_table = read.table('DATA/RAW/rice_wr//sequence_report.tsv', sep='\t', header=T)
chrom = as.numeric(chrom_table$Chromosome.name[1:12])
names(chrom) = chrom_table$RefSeq.seq.accession[1:12]

### convert chromosome id
GENES$V1 = chrom[GENES$V1]


## remove genes outside main chromosomes 
GENES=na.omit(GENES)


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/rice_wr/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)

