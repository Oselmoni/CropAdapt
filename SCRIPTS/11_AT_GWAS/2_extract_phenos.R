### Load phenotypes datasets
pheno_ds = list.dirs('DATA/ATHALIANA/RAW/pheno_database', full.names=F)[-1]


PHENO = data.frame() # dataframe to store phenotypes
PHENO.META = data.frame() # dataframe to store phenotypes metadata

### For every dataset, retrieve phenotypes
for (ds in pheno_ds) {

  # load phenotypes
  pheno_raw = read.csv(paste0('DATA/ATHALIANA/RAW/pheno_database/',ds,'/study_',ds,'_values.csv'), check.names = F)
  
  # if more than one phenotype per accession -> get average
  pheno_i = do.call(rbind, as.list(by(pheno_raw[,-c(1:2),drop=F], pheno_raw[,1], function(x) {apply(x,2,mean,na.rm=T)})))
  
  # reformat phenotypes ID 
  colnames(pheno_i) = paste0(ds,'_',colnames(pheno_raw)[-c(1:2)])
  
  # add to full table
  PHENO[rownames(pheno_i),colnames(pheno_i)] = NA
  PHENO[rownames(pheno_i),colnames(pheno_i)] = pheno_i

  # load metadata
  meta_i = (read.csv(paste0('DATA/ATHALIANA/RAW/pheno_database/',ds,'/study_',ds,'_phenotypes.csv')))
  rownames(meta_i) = paste0(ds,'_',meta_i$name)
  
  PHENO.META = rbind(PHENO.META, meta_i)
}

# remove empty spaces from names
rownames(PHENO.META) = gsub('(.*) $', '\\1', rownames(PHENO.META))
colnames(PHENO) = gsub('(.*) $', '\\1', colnames(PHENO))


# check that IDs are identical
PHENO = PHENO[,rownames(PHENO.META)]

save(PHENO, file='DATA/ATHALIANA/PHENO.rda')
save(PHENO.META, file='DATA/ATHALIANA/PHENO_META.rda')


