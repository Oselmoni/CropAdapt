library(missForest)
source('SCRIPTS/custom_R_functions.R')

### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$dataset
ds_meta = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t', row.names=1)

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



#######
####### Check imputation accuracy
#######

### check which varibles were missing
apply(ENVmissing, 2, function(x) {mean(is.na(x))}) # only variables from soil grids

navars = envVars$VariableID[envVars$Source=='SoilGrids'] # create index of variables to cross-validate


# keep only sites without missing data, then scale each variable
ENVscaled = na.omit(ENV)
scENV = scale(ENVscaled[,-c(1,2,ncol(ENVscaled))], center = T, scale=T)
ENVscaled[,colnames(scENV)] = scENV


# prepare an environmental dataset where to add masked missing data
ENVmissing = ENVscaled

# container of masked sites x env
storeMissing = list()

set.seed(0);for (ds in ds_list) {
 
  # extract scaled env data
  sENV = as.matrix(ENVmissing[ENVmissing$ds==ds,-c(1,2,ncol(ENVmissing))])
  
  # create an index of sites with 10 km of distance
  sites= paste0('s_',cutree(hclust(dist(ENVmissing[ENVmissing$ds==ds,c('LON','LAT')])), h=0.01))
  
  # sample 100 variables and sites to be masked
  N = round(length(unique(sites))*0.2)*length(navars)
  na_sites = sample(unique(sites), size = N, replace = T)
  na_evar = sample(navars, size = N, replace = T)
  
  
  for (i in 1:N) {
    sENV[sites==na_sites[i],na_evar[i]] = NA
  }
  
  storeMissing[[ds]] = list(na_sites, na_evar)
  
  ### build env matrix with missing points
  ENVmissing[ENVmissing$ds==ds,-c(1,2,ncol(ENV))] = sENV
  
}


# Check that missing rate of masked sites ~ missing date of real data
apply(ENV, 2, function(x) {mean(is.na(x))}) # only variables from soil grids
apply(ENVmissing, 2, function(x) {mean(is.na(x))}) # only variables from soil grids


### Run imputation on masked sites
set.seed(0);iENVmissing=missForest(ENVmissing)$ximp



### Check how imputation accuracy
ERRORS = data.frame()

for (ds in ds_list) {

  # load original vs imputed environmental matrices
  sENV = as.matrix(ENVscaled[ENVscaled$ds==ds,-c(1,2,ncol(ENV))])
  sENVi = as.matrix(iENVmissing[ENVscaled$ds==ds,-c(1,2,ncol(ENV))])
  
  # create an index of sites with 10 km of distance
  sites= paste0('s_',cutree(hclust(dist(ENVmissing[ENVmissing$ds==ds,c('LON','LAT')])), h=0.1))
  

  # retrieve indices of masked sites and evars
  N = round(length(unique(sites))*0.2)*length(navars)
  
  na_sites = storeMissing[[ds]][[1]]
  na_evar = storeMissing[[ds]][[2]]
  
  for (i in 1:N) {
 
    ERRORS = rbind(ERRORS, data.frame(ds, 'evar'=na_evar[i], 'REAL'=mean(sENV[sites==na_sites[i],na_evar[i]]), 'IMP'=mean(sENVi[sites==na_sites[i],na_evar[i]])))
    
  }
  

}

### Plot errors by datasets
png(paste0('FIGURES/envVars/IMP_ds.png'), h=6, w=6, units = 'in', res=500)
par(mfrow=c(4,5));par(mar=c(3,3,1,1))
for (ds in sort(ds_list)) {

  s_ERRORS = na.omit(ERRORS[ERRORS$ds==ds,] )


  cplot(s_ERRORS$REAL, s_ERRORS$IMP, pch=16, cex=0.75, col=adjustcolor(1,0.5), xlab = 'real', ylab='imputed')
  title(main=ds_meta[ds,'name'], cex.main=0.75)

  legend('bottomright', paste0('R=',signif(cor(s_ERRORS$REAL, s_ERRORS$IMP, use = 'pairwise.complete.obs'),3)), cex=0.75, box.lwd=0, bg=NA)
  
  
  
}
dev.off()

### Plot errors by env vars
rownames(envVars) = envVars$VariableID

png(paste0('FIGURES/envVars/IMP_env.png'), h=6, w=6, units = 'in', res=500)
par(mfrow=c(4,3));par(mar=c(3,3,1,1))
for (evar in navars) {
  
  s_ERRORS = na.omit(ERRORS[ERRORS$evar==evar,] )
  
  cplot(s_ERRORS$REAL, s_ERRORS$IMP, pch=16, cex=0.75, col=adjustcolor(1,0.5), xlab = 'real', ylab='imputed')
  title(main=envVars[evar,'Description'], cex.main=0.75)
  
  legend('bottomright', paste0('R=',signif(cor(s_ERRORS$REAL, s_ERRORS$IMP, use = 'pairwise.complete.obs'),3)), cex=0.8, box.lwd=0, bg=NA)
  
}
dev.off()



### Make a table of missing environmental data per dataset

MissEnv = data.frame(do.call(rbind, by(ENV[,navars], ENV$ds, function(x) {apply(x,2,function(y){mean(is.na(y))})})))
rownames(MissEnv) = ds_meta[rownames(MissEnv),'name']
colnames(MissEnv) = envVars[navars,'Description']

MissEnv$TOT = apply(MissEnv,1,mean)
MissEnv['TOT',] = apply(MissEnv,2,mean)

MissEnv = signif(MissEnv,2) 

write.table(MissEnv, 'FIGURES/envVars/MissEnv.txt', col.names = T, row.names=T, quote=F)
