source('SCRIPTS/10_GEA/picmin.R') # load custom script to run picmin

### Load same env variable for different datasets
dslist = list.dirs('DATA/GEA_OUTPUT/LFMM/', full.names = F)[-1]

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
  
  for (ds in dslist[1:11]) {

    load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/LFMM_OG_',evar,'.rda'))
    
    
    
    ## transform to empirical pvalues
    ep = PicMin:::EmpiricalPs(LFMM_OG$minP)
    

    #PVALS[out$gene,ds] = out$Z_pVal
    PVALS[rownames(LFMM_OG),ds] = ep
    
  }
  
  ### Run PICMIN
  PM_RES[[evar]] = RunPicmin(PVALS, numReps = 1000)
  
  print(sum(PM_RES[[evar]]$pooled_q<0.3))

  
}

### 




lapply(PM_RES, function(x) {
  sum(x$p<0.05)})

XXX = PM_RES$TOPO

XXX[which.min(XXX$p),]
