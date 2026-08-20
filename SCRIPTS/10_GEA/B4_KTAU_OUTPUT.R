library(SNPRelate)
library(rnaturalearth)
library(terra)
source('SCRIPTS/custom_R_functions.R')

### Load results of PICMIN for real and permuted env
load('DATA/GEA_OUTPUT/KTAU/PM_RES.rda')
load('DATA/GEA_OUTPUT/KTAU/PERM_PM_Q.rda')


### Load slist of datasets
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$dataset

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID

### Get list of selected datasets
selectedDS = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, row.names=1, sep='\t')

### Load orthogroups list
load('DATA/SNP_ANNOTATION/orthogroups.rda')
rownames(orthogroups) = orthogroups$Orthogroup
### get land for plotting
land = ne_countries()


## Set q-value significance threshold for PicMin
qt = 0.2


###
#### Calculate number of convergent signals per environmentl variable
###
conv = unlist(lapply(PM_RES, function(x) {sum(x$pooled_q<qt, na.rm=T)}))



### top variables
topenv = names(which((conv>0)))


###
#### Get detail about candidate genes
###
TOPGENES = data.frame()

for (env in topenv) {

  ### retrieve list of signficant genes
  top_genes = PM_RES[[env]][PM_RES[[env]]$pooled_q<qt,]
  
  ### add env to list
  top_genes$ENV = env
  
  ### add env description
  top_genes$ENVD = envVars[env,'Description']
  
  ### retrieve p-values of genes
  PVALS = data.frame()
  
  for (ds in dslist) {
    
    load(paste0('DATA/GEA_OUTPUT/KTAU/',ds,'/',env,'.rda'))
    
    ## transform to empirical pvalues
    ep = PicMin:::EmpiricalPs(WZA[[2]]$Z_pVal)
    
    PVALS[WZA[[2]]$gene,ds] = ep
    
  }
  
  
  
  ### for every top gene... 
  for (i in 1:nrow(top_genes)) {
    
    geneID = top_genes$locus[i] # get geneID
    
    ### retrieve p-values
    top_genes[i,paste0('P(',colnames(PVALS),')')] = PVALS[geneID,]
    
  }
  
  ### add to global container
  TOPGENES = rbind(TOPGENES, top_genes)
}



TOPGENES = TOPGENES[order(TOPGENES$pooled_q),]
rownames(TOPGENES) = 1:nrow(TOPGENES)
TOPGENES

TOPGENES$OG = orthogroups[TOPGENES$locus,'Arabidopsis_thaliana.TAIR10.pep.all']
TOPGENES

# save TOPGENES object
save(TOPGENES, file='DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')

### Output topgenes table

OUT = data.frame('orthogeneID'=TOPGENES$locus,
                 'ATgene'= unlist(lapply(strsplit(TOPGENES$OG, ', '), function(x) {
                   paste(unique(gsub('(.*)\\..*','\\1',x)), collapse=',')
                 })),
                 'q'=TOPGENES$pooled_q,
                 'env'=TOPGENES$ENVD,
                 'N_rep'=TOPGENES$n_est,
                 'N_not_NA'=apply(TOPGENES[9:22], 1, function(x) {sum(is.na(x)==F)}),
                 'crops'=apply(TOPGENES, 1, function(x) {
                   N = x['n_est']
                   crops = x[9:22]
                   crops = sort(crops)
                   crops = names(crops[1:N])
                   crops = selectedDS[substr(crops, 3, nchar(crops)-1),'name']
                   return(paste(crops, collapse=','))
                 })
)

write.table(OUT, file='FIGURES/GEAs/KTAU_picmin_topgenes.tsv', sep='\t', col.names=T, row.names = F, quote=F)

