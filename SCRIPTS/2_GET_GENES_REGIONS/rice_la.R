### Load gtf file
GTF = read.table('DATA/RAW/rice_la/genomic.gtf', sep='\t')

table(GTF$V3)

### keep only genes
GENES = GTF[GTF$V3=='transcript',c(1,4,5)]

table(GENES$V1)

### Create chromsome id
chrom_table = read.table('DATA/RAW/rice_la///sequence_report.tsv', sep='\t', header=T)
chrom = 1:12
names(chrom) = chrom_table$RefSeq.seq.accession[1:12]

### convert chromosome id
GENES$V1 = chrom[GENES$V1]



### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/rice_la/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)

