library(lfmm)
library(missForest)
library(qvalue)

### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$dataset
annotation_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$annotation_file
protein_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T)$protein_file

names(annotation_files)=names(protein_files)=ds_list

### get list of orthologous groups
load('DATA/SNP_ANNOTATION/orthogroups.rda')



### for every dataset...
for (ds in ds_list) {
 
  #ds=ds_list[1] 

  print(ds)

  ### Load imputed GT matrix
  gti_files=list.files('DATA/GEA_INPUT/GTI/',full.names = T)
  load(gti_files[regexpr(ds, gti_files)!=-1])

  dim(GTI)
  
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

  ### Load env data
  load(paste0('DATA/GEA_INPUT/meta_env/meta_',ds,'.rda'))
  ENV = meta[,-c(1:4)]
  
  ### Impute missing environmental data
  iENV=missForest(ENV)
  iENVdf=iENV$ximp
  
  ### retrieve number of K from GTI filename
  K=as.numeric(gsub('.*K(.*)\\.rda','\\1', gti_files[regexpr(ds, gti_files)!=-1]))
  
  ### Create container for LFMM results
  dir.create(paste0('DATA/GEA_OUTPUT/LFMM/',ds), showWarnings = F)
  
  ### Run LFMM for one environmental variable at the time
  par(mfrow=c(3,3))
  for (e in colnames(iENVdf)) {
    
    print(e)
    
    lfmm = lfmm_ridge(GTI, iENVdf[,e], K = K)
    lfmm.pv <- lfmm_test(Y=GTI, X=iENVdf[,e], lfmm=lfmm, calibrate="gif")
    
    ### add lfmm results to container
    LFMM=lfmm.pv
    
    ### get lowest p-value per orthogroup
    SNP_ANNOTATION_OL$p = LFMM$calibrated.pvalue[SNP_ANNOTATION_OL$id,1]
    
    
    LFMM_OG = data.frame('meanP' = by(SNP_ANNOTATION_OL$p, SNP_ANNOTATION_OL$rnd_og, mean, na.rm=T),
                         'medP' = by(SNP_ANNOTATION_OL$p, SNP_ANNOTATION_OL$rnd_og, median, na.rm=T),
                         'minP' = by(SNP_ANNOTATION_OL$p, SNP_ANNOTATION_OL$rnd_og, min, na.rm=T))
    

    ### Store LFMM results
    save(LFMM_OG, file=paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/LFMM_OG_',e,'.rda'))
    save(LFMM, file=paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/LFMM_',e,'.rda'))
    

  }
  
 
}
  
  