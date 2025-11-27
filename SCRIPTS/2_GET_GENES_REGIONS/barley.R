### Load gtf file
GTF = read.table('DATA/RAW/barley/Hv_IBSC_PGSB_r1_HighConf.gtf')

### keep only genes
GENES = GTF[GTF$V3=='transcript',c(1,4,5)]

### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/barley/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)
