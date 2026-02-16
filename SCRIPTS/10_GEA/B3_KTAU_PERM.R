### This script shuffles environmental data across sites, and repeats WZA GEA. This is repeated 20x per environmental variable. 
### The goal is to obtain a null p-value distribution for PicMin.

library(foreach)
library(doParallel)
library(SNPRelate)
source('SCRIPTS/custom_R_functions.R')

### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$dataset
annotation_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$annotation_file
protein_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$protein_file

names(annotation_files)=names(protein_files)=ds_list

### get list of orthologous groups
load('DATA/SNP_ANNOTATION/orthogroups.rda')

#setup parallel backend to use many processors
cores = min(c(detectCores(),120))
cl <- makeCluster(cores-2) 
registerDoParallel(cl)


### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')

## Create set of environmental variables to be permuted
set.seed(0);perm_env_vars = sample(envVars$VariableID, 100, replace=T)


### 
###### Step 1: Run WZA with permuted environmental variables
###


### for every dataset...
for (ds in ds_list[1:17]) {
  
  dir.create(paste0('DATA/GEA_OUTPUT/KTAU/',ds), showWarnings = F)
  
  gc()
  #ds=ds_list[17]
  
  print(ds)
  
  ### Load imputed GT matrix
  gti_files=list.files('DATA/GEA_INPUT/GTI/',full.names = T)
  load(gti_files[regexpr(ds, gti_files)!=-1])
  
  
  ######
  ######### LOAD AND PREPARE ANNOTATION DATA
  ######
  
  
  ### find ortologous genes
  
  # load SNP to gene table
  load(paste0('DATA/SNP_ANNOTATION/',annotation_files[ds],'.rda'))
  
  # keep only SNPs used in GEA, to speed up
  SNP_ANNOTATION = SNP_ANNOTATION[SNP_ANNOTATION$snp%in%colnames(GTI),]
  
  # keep only SNPs with a gene annotation
  SNP_ANNOTATION = SNP_ANNOTATION[is.na(SNP_ANNOTATION$geneID)==F,]
  GTI  = GTI[,colnames(GTI)%in%SNP_ANNOTATION$snp]
  
  # find ortogroup corresponding to gene annotation file
  species_OG = orthogroups[,protein_files[ds]]
  names(species_OG) = orthogroups$Orthogroup
  
  # remove empty ortogroups
  species_OG = species_OG[species_OG!='']
  
  # map genes to ortogroups
  GENE_TO_OG = rep(names(species_OG), lapply(species_OG, function(x) {length(unlist(strsplit(x, ', ')))}))
  names(GENE_TO_OG) = unlist(lapply(species_OG, strsplit, ', '))
  
  # add orthogroup column to SNP annotation
  SNP_ANNOTATION$OG = GENE_TO_OG[SNP_ANNOTATION$geneID]
  
  ### Keep only SNPs with ortogroup annotation
  SNP_ANNOTATION = SNP_ANNOTATION[is.na(SNP_ANNOTATION$OG)==F,]
  GTI  = GTI[,colnames(GTI)%in%SNP_ANNOTATION$snp]
  SNP_ANNOTATION = SNP_ANNOTATION[SNP_ANNOTATION$snp%in%colnames(GTI),]
  
  ### REFORMAT SNP ANNOTATION: one SNP per line
  SNP_ANNOTATION_OL = do.call(rbind, by(SNP_ANNOTATION, SNP_ANNOTATION$snp, function(x) {
    
    id = unique(x$snp)
    chr = unique(x$chr)
    pos = unique(x$pos)
    
    genes = paste(unique(x$geneID), collapse=';')
    Ngenes = length(unique(x$geneID))
    
    orthog = paste(unique(x$OG), collapse=';')
    Northog = length(unique(x$OG))
    
    ## pick most frequent OG
    og_count = table(x$OG)
    max_og_count = names(which(og_count==max(og_count))) # og with maximal count
    rnd_og = sample(max_og_count, 1) ## in case there are more than 1 og with maximal counts, pick one randomly
    
    return(data.frame(id,chr,pos,genes,Ngenes, orthog, Northog, rnd_og))
    
  }))
  
  rownames(SNP_ANNOTATION_OL) = SNP_ANNOTATION_OL$id
  

  ######
  ######### GET MINOR ALLELE FREQUENCY
  ######
  
  
  ### calculate minor allele frequency of every SNP
  SNP_ANNOTATION_OL = SNP_ANNOTATION_OL[colnames(GTI),]
  SNP_ANNOTATION_OL$MAF =  apply(GTI, 2, function(gt) { 
    
    af = sum(gt, na.rm=T)/(length(gt[is.na(gt)==F])*2) # calculate bi-allelic frequency
    maf= min(c(af,1-af)) # double check which allele is minor
    return(maf)             
  })
  
  
  
  
  ######
  ######### LOAD ENVIRONMENTAL DATA
  ######
  
  
  ### Load imputed env data
  load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  iENV = iENV[rownames(GTI),]
  
  
  
  ### Run LFMM for one environmental variable at the time
  pWZA =  foreach (i=1:length(perm_env_vars), .combine = cbind) %dopar%  {
    
    e = perm_env_vars[i]
    
    ### Get ID of sampling sites (within 10 km)
    coord = iENV[,1:2]
    sites = paste0('site_',cutree(hclust(dist(coord)), h = 0.05))
    names(sites) = rownames(iENV)
    
    ### get mean environmental value per sampling site
    sites_env = by(iENV[,e], sites, mean, na.rm=T)
    
    ### shuffle sites
    sh_sites = sample(unique(sites))
    names(sh_sites) = unique(sites)
    
    ### shuffle environmental data across sites
    shuffledENV = matrix(sites_env[sh_sites[sites]], ncol=1)
    
    ### Create container for WZA input 
    WZAin = SNP_ANNOTATION_OL
    
    ### Calculate correlation tau kendall for every SNP
    kTAU = -abs(cor(shuffledENV, GTI, use = 'pairwise.complete.obs', method='kendall'))
    
    ### calculate empirical p-values
    pKT =  rank(kTAU)/length(kTAU)
    
    # add pKT to WZA input
    WZAin$pKT = pKT
    
    
    ### write prepare temporary files for WZA run
    tmpIN = tempfile() # create temporary file for input
    tmpOUT = tempfile() # create temporary file for output
    
    write.table(WZAin, tmpIN, quote=F, row.names=F, sep='\t')
    
    ### run WZA via python
    system(paste0('/home/oselmo/data/conda/envs/myenv/bin/python SCRIPTS/10_GEA/WZA/general_WZA_script.py ',
                  '--correlations ',tmpIN,' ', ### input file
                  '--summary_stat pKT --window rnd_og --MAF MAF ', ### other params
                  '--output ',tmpOUT)) ### output folder
    
    
    ### read wza out
    WZAout = read.csv(tmpOUT)
    
    ## add p-value to container
    pvals = WZAout$Z_pVal
    names(pvals) = WZAout$gene
    return(pvals)
    
  }
  
  ### write output
  save(pWZA, file=paste0('DATA/GEA_OUTPUT/KTAU/',ds,'/pWZA.rda'))    
  
  
}





### 
###### Step 2: Run PicMin across species
###

### get picmin functions
source('SCRIPTS/10_GEA/picmin.R') # load custom script to run picmin


### Precompute null distribution of p-values for different number of lineages
nullP = PicMinNull(linMin = 3, linMax = length(ds_list))



### Create container of picmin results
PERM_PM_RES = list()


### for every replicate
for (i in 1:length(perm_env_vars)) {
  gc()
  print(i)
  
  ####
  ####### Extract p-values for GEA vs. variable of interest
  ####
  
  PVALS = data.frame()
  
  for (ds in ds_list) {
    
    load(paste0('DATA/GEA_OUTPUT/KTAU/',ds,'/pWZA.rda'))
    
    ## get current permutation
    wza_pval = pWZA[,i]
    
    ## transform to empirical pvalues
    ep = PicMin:::EmpiricalPs(wza_pval)
    
    PVALS[names(wza_pval),ds] = ep
    
  }
  
  ### Run PICMIN
  PERM_PM_RES[[i]] = RunPicmin(PVALS, nullP = nullP)
  
  
  ### Display number of hits in random datasets  
  hist(unlist(lapply(PERM_PM_RES, function(x) {sum(x$pooled_q<0.01)})), xlab='# hits by chance')
  
}


### create container of q-values of permuted-wza-picmin
PERM_PM_Q = data.frame()
for (i in 1:length(perm_env_vars)) {
  
  PERM_PM_Q[PERM_PM_RES[[i]]$locus,paste0(i,'_',perm_env_vars[i])] = PERM_PM_RES[[i]]$pooled_q
  
}

save(PERM_PM_Q, file='DATA/GEA_OUTPUT/KTAU/PERM_PM_Q.rda')




stopCluster(cl)
