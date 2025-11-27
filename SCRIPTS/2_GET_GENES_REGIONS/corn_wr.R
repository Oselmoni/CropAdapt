### Load gtf file
GTF = read.table('DATA/RAW/corn_wr/genomic.gff', sep='\t')

### keep only genes
GENES = GTF[GTF$V3=='mRNA',c(1,4,5)]

### Create chromsome id
chrom_table = read.table('DATA/RAW/corn_wr//sequence_report.tsv', sep='\t', header=T)
chrom = 1:10
names(chrom) = chrom_table$RefSeq.seq.accession[1:10]

### convert chromosome id
GENES$V1 = chrom[GENES$V1]

## remove genes outside main chromosomes 
GENES=na.omit(GENES)


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/corn_wr/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)

