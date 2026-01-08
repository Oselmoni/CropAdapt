source('SCRIPTS/10_GEA/picmin.R') # load custom script to run picmin

### Load same env variable for different datasets
dslist = list.dirs('DATA/GEA_OUTPUT/WZA/', full.names = F)[-1]

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')

### Create container of picmin results
PM_RES = list()

for (evar in envVars$VariableID) {
  
  print(evar)
  

  ####
  ####### Extract p-values for GEA vs. variable of interest
  ####
  
  PVALS = data.frame()

  for (ds in dslist) {
    
    load(paste0('DATA/GEA_OUTPUT/WZA/',ds,'/',evar,'.rda'))
    
    out = WZA[[2]]
    
    ## transform to empirical pvalues
    ep = PicMin:::EmpiricalPs(out$Z_pVal)
    
    #PVALS[out$gene,ds] = out$Z_pVal
    PVALS[out$gene,ds] = ep
    
  }

  ### Run PICMIN
  PM_RES[[evar]] = RunPicmin(PVALS, numReps = 20)
  
  print(sum(PM_RES[[evar]]$p<0.05))
  
}

### 

lapply(PM_RES, function(x) {
  sum(x$p<0.05)})

XXX = PM_RES$TOPO

XXX[which.min(XXX$p),]
