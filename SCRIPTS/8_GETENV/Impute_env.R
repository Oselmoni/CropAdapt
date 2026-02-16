library(missForest)

### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$dataset

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')


### Create global container of environmental variables
ENV = data.frame()

### Load all environmental matrices
for (ds in ds_list) {
  
  load(paste0('DATA/GEA_INPUT/meta_env/meta_',ds,'.rda'))
  
  meta$ds = as.factor(ds)
  
  ENV = rbind(ENV, meta[,c('LON','LAT',envVars$VariableID,'ds')])

}


### Impute missing environmental variables
set.seed(0);iENV=missForest(ENV)
iENVdf=iENV$ximp



### Write output for every dataset
for (ds in ds_list) {
  
  iENV = iENVdf[iENVdf$ds==ds,-ncol(iENVdf)]
  
  save(iENV, file=paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  
}
