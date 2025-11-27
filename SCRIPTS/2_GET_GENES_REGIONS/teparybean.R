### Load gtf file
GTF = read.table('DATA/RAW/teparybean//Pacu.CVR.gene_models.hc.gff3', sep='\t')


### keep only genes
GENES = GTF[GTF$V3=='mRNA',c(1,4,5)]

### Keep only genes on chromosomes
GENES = GENES[substr(GENES$V1, 1,3) =='Chr',]

### Modify format of choromsomes
GENES$V1 = (substr(GENES$V1, 4,5))

table(GENES$V1)


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/teparybean//GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)


