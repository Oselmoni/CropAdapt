### Load gtf file
GTF = read.table('DATA/RAW/soybean/Gmax_275_Wm82.a2.v1.gene_exons.gff3', sep='\t')


### keep only genes
GENES = GTF[GTF$V3=='mRNA',c(1,4,5)]


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/soybean/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)

