## load datasets stats
DS_STATS = read.table('FIGURES/ds_genomic_summary/DS_STATS.txt')
DS_STATS_RAW = read.table('FIGURES/ds_genomic_summary/DS_STATS_RAW.txt')

## open datasets list
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, sep='\t')


## merge two table
OUT = DS_STATS[dslist$dataset,]
OUT = cbind(OUT, DS_STATS_RAW[dslist$annotation_file,])
OUT$Dataset = dslist$name

## sort columns
OUT = OUT[,c('Dataset','n_ind','n_loci','n_loci_aLD','N','N_snps','N_geo','N_AF_mn','N_snps_AF_mn','N_snps_AF_mn_maf','SamSites','med_D_SamSites')]

## round distance between samples
OUT$med_D_SamSites = round(OUT$med_D_SamSites)


## export table
write.table(OUT, 'FIGURES/ds_genomic_summary/DS_STATS_FINAL.txt', row.names=F, col.names=T, quote=F, sep='\t')
