library(lfmm)
library(SNPRelate)
library(foreach)
library(doParallel)

### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$dataset
annotation_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$annotation_file
protein_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$protein_file

names(annotation_files)=names(protein_files)=ds_list

### get list of orthologous groups
load('DATA/SNP_ANNOTATION/orthogroups.rda')

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')

#setup parallel backend to use many processors
cores = min(c(detectCores(),120))
cl <- makeCluster(cores[1]-2) 
registerDoParallel(cl)


### for every dataset...
for (ds in ds_list) {
  gc()
  #ds=ds_list[1]
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
  
  
  
  ### Load imputed env data
  load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  iENV = iENV[rownames(GTI),-c(1:2)]
  
  
  
  ### retrieve number of K from GTI filename
  K=as.numeric(gsub('.*K(.*)\\.rda','\\1', gti_files[regexpr(ds, gti_files)!=-1]))
  
  ### Create container for LFMM results
  dir.create(paste0('DATA/GEA_OUTPUT/LFMM/',ds), showWarnings = F)
  
  ### save SNP annotation
  save(SNP_ANNOTATION_OL, file=paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/SNP_ANNOTATION_OL.rda'))
  
  ### Run LFMM for one environmental variable at the time
  foreach(e=colnames(iENV)) %dopar% {
    
    library(lfmm)
    
    print(e)
    
    lfmm = lfmm_ridge(GTI, iENV[,e], K = K)
    lfmm.pv <- lfmm_test(Y=GTI, X=iENV[,e], lfmm=lfmm, calibrate="gif")
    
    ### assign p-values to every SNP
    SNP_ANNOTATION_OL$pLFMM = lfmm.pv$calibrated.pvalue[SNP_ANNOTATION_OL$id,1]
    
    ### Create container for WZA input 
    WZAin = SNP_ANNOTATION_OL
    
    
    ### write prepare temporary files for WZA run
    tmpIN = tempfile() # create temporary file for input
    tmpOUT = tempfile() # create temporary file for output
    
    write.table(WZAin, tmpIN, quote=F, row.names=F, sep='\t')
    
    ### run WZA via python
    system(paste0('/home/oselmo/data/conda/envs/myenv/bin/python SCRIPTS/10_GEA/WZA/general_WZA_script.py ',
                  '--correlations ',tmpIN,' ', ### input file
                  '--summary_stat pLFMM --window rnd_og --MAF MAF ', ### other params
                  '--output ',tmpOUT)) ### output folder
    
    
    ### read wza out
    WZAout = read.csv(tmpOUT)
    
    ### write output
    PVALS = WZAin$pLFMM;names(PVALS)=WZAin$id
    WZA=list(PVALS, WZAout)
    save(WZA, file=paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/',e,'.rda'))   
    
  }
  
  
}


stopCluster(cl)