### Load phenotypes datasets
pheno_ds = list.dirs('DATA/ATHALIANA/RAW/pheno_database', full.names=F)[-1]


PHENO.META = data.frame() # dataframe to store phenotypes metadata

### For every dataset, retrieve phenotypes
for (ds in pheno_ds) {
  
  # load metadata
  meta_i = (read.csv(paste0('DATA/ATHALIANA/RAW/pheno_database/',ds,'/study_',ds,'_phenotypes.csv')))
  rownames(meta_i) = meta_i$phenotype_id
  
 
  
  
  PHENO.META = rbind(PHENO.META, meta_i)
}




save(PHENO.META, file='DATA/ATHALIANA/PHENO_META.rda')

