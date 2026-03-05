library(vcfR)

SNPS_IND = data.frame()

#### Get stats of unfiltered vcf file
for (ds in list.dirs('DATA/PROCESSED_VCF/',  recursive = F, full.names = F)) {
  
  print(ds)
  
  ## get path of vcf file
  vcfpath  = paste0('DATA/PROCESSED_VCF/',ds,'/processed.vcf.gz')
  
  ## retrieve number of loci
  n_loci <- as.numeric(system(paste0("zcat ",vcfpath," | grep -v '^#' | wc -l"), intern = TRUE))

  ## retrieve number of individuals
  vcf = read.vcfR(vcfpath, nrows = 1)
  
  n_ind = ncol(vcf@gt)
  
  ## add to container
  SNPS_IND[ds,'n_loci'] = n_loci
  SNPS_IND[ds,'n_ind'] = n_ind
  
  ## save intermediary file
  write.table(SNPS_IND, file = 'FIGURES/ds_genomic_summary/DS_STATS_RAW.txt', sep='\t', col.names = T, row.names=T, quote=F)
  
}


#### Get stats of ld-pruned vcf file
for (ds in list.dirs('DATA/PROCESSED_VCF/',  recursive = F, full.names = F)) {
  
  ## get path of pruned vcf file
  vcfpath  = paste0('DATA/PROCESSED_VCF/',ds,'/processed.vcf.gz_pruned.vcf.gz')
  
  ## retrieve number of loci
  n_loci <- as.numeric(system(paste0("zcat ",vcfpath," | grep -v '^#' | wc -l"), intern = TRUE))
  
  ## add to container
  SNPS_IND[ds,'n_loci_aLD'] = n_loci

  ## save intermediary file
  write.table(SNPS_IND, file = 'FIGURES/ds_genomic_summary/DS_STATS_RAW.txt', sep='\t', col.names = T, row.names=T, quote=F)
  
}

