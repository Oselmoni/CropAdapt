source('SCRIPTS/10_GEA/picmin.R') # load custom script to run picmin

### Load same env variable for different datasets
dslist = list.dirs('DATA/GEA_OUTPUT/WZA/', full.names = F)[-1]

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) 
registerDoParallel(cl)

### Precompute null distribution of p-values for different number of lineages
nullP = PicMinNull(linMin = 3, linMax = length(dslist))


### Create container of picmin results
PM_RES = list()

### For every environmental variable...
for (evar in envVars$VariableID) {
  
  print(evar)
  
  ####
  ####### Extract p-values for GEA vs. variable of interest
  ####
  
  PVALS = data.frame()
  
  for (ds in dslist) {
    
    load(paste0('DATA/GEA_OUTPUT/WZA/',ds,'/',evar,'.rda'))
    
    ## transform to empirical pvalues
    ep = PicMin:::EmpiricalPs(WZA[[2]]$Z_pVal)
    
    #PVALS[out$gene,ds] = out$Z_pVal
    PVALS[WZA[[2]]$gene,ds] = ep
    
  }
  
  ### Run PICMIN
  PM_RES[[evar]] = RunPicmin(PVALS, nullP = nullP)
  
  
}


### save output
save(PM_RES, file='DATA/GEA_OUTPUT/WZA/PM_RES.rda')


### Stop cluster
stopCluster(cl)

