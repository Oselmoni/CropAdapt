library(SNPRelate)
library(rnaturalearth)
library(terra)
source('SCRIPTS/custom_R_functions.R')

### Load results of PICMIN for real and permuted env
load('DATA/GEA_OUTPUT/LFMM/PM_RES.rda')
load('DATA/GEA_OUTPUT/LFMM/PERM_PM_Q.rda')


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
    
    load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/',env,'.rda'))
    
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
save(TOPGENES, file='DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')

### Output topgenes table

OUT = data.frame('orthogeneID'=TOPGENES$locus,
                 'ATgene'= unlist(lapply(strsplit(TOPGENES$OG, ', '), function(x) {
                   paste(unique(substr(x, 1, nchar(x)-2)), collapse=',')
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

write.table(OUT, file='FIGURES/GEAs/LFMM_picmin_topgenes.tsv', sep='\t', col.names=T, row.names = F, quote=F)
  



  
i=1
{
  
  ## Isolate gene of interest
  tg = TOPGENES[i,]
  
  ## Retrieve datasets in which gene is significant
  sortedPs = tg[1,9:length(tg)][order(as.numeric(tg[1,9:length(tg)]))]
  
  top_ds = names(sortedPs[1:tg$n_est])
  top_ds = substr(top_ds, 3,nchar(top_ds)-1)
  
  
  #####
  ######## DO MANHATTAN PLOTS OF SIGNIFICANT GENES ACROSS CROPS
  #####
  
  par(mfrow=c(2,3))
  
  for (ds in top_ds) {
    
    ## retrieve LFMM input & output
    load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/',tg$ENV,'.rda'))
    WZAout = WZA[[2]]
    rownames(WZAout) = WZAout$gene

    ## transform to empirical p-values
    WZAout$ep = PicMin:::EmpiricalPs(WZAout$Z_pVal)
    
    ## load gene annotations
    load(paste0('DATA/GENE_ANNOTATION/genes_MP_',ds,'.rda'))
    
    ## retain only genes with OG used in analysis 
    genes_MP = genes_MP[genes_MP$OG%in%WZAout$gene,]
    
    
    # find empirical p-value for every gene
    genes_MP$EP = -log(WZAout[genes_MP$OG,'ep'], 10)
    
    # load chromosome id
    metaCHR = read.csv(paste0('DATA/REF_TO_VCF/',selectedDS[ds,'annotation_file'],'.csv'))
    chIDX = metaCHR$id
    names(chIDX) = metaCHR$reference
    
    # make chromosom numerical
    genes_MP$CHRN = chIDX[genes_MP$CHR]
    
    # order genes by chromosome, then by position
    genes_MP = genes_MP[order(genes_MP$CHRN, genes_MP$POS, decreasing = F),]
    
    
    ### Do Manhattan plot
    manhattanPlot(p=genes_MP$EP, chr=genes_MP$CHRN, pos=genes_MP$POS, sig=which(genes_MP$OG==tg$locus), main=paste0(ds,'\np=',    signif(10^-genes_MP$EP[which(genes_MP$OG==tg$locus)][1],3)), chrL=metaCHR$id)
    

  }
}




## Isolate gene of interest
i=1
tg = TOPGENES[i,]
## Retrieve datasets in which gene is significant
sortedPs = tg[1,9:length(tg)][order(as.numeric(tg[1,9:length(tg)]))]

top_ds = names(sortedPs[1:tg$n_est])
top_ds = substr(top_ds, 3,nchar(top_ds)-1)


for (ds in top_ds) {
  

  ## retrieve WZA input & output
  load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/SNP_ANNOTATION_OL.rda'))
  load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/',tg$ENV,'.rda'))
  
  ##
  WZAin = SNP_ANNOTATION_OL
  WZAin$pval = WZA[[1]][WZAin$id]
  WZAout = WZA[[2]]
  rownames(WZAout) = WZAout$gene
  
  
  ## find SNPs on gene
  WZAin_gene = WZAin[WZAin$rnd_og==tg$locus,]
  
  
  ## find genotype variation
  ### Load GT matrix
  SNPS = snpgdsOpen(paste0('DATA/GEA_INPUT/GDS/gds_',ds,'.gds'), readonly = T, allow.duplicate = T)
  
  ### Load imputed GT matrix
  gti_files=list.files('DATA/GEA_INPUT/GTI/', pattern = paste0('GTI_',ds,'_K'),full.names = T)
  load(gti_files[regexpr(ds, gti_files)!=-1])
  
  
  ### Load environmental data
  load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  
  ### calculate score for every SNP in gene, using the Stouffer weighted transformation (normal deviate of pvalue, weighted by allele frequencies)
  WZAin_gene$Z = qnorm(1 - WZAin_gene$pval)*WZAin_gene$MAF*(1-WZAin_gene$MAF)
  
  ### Plot GEA
  gt = GTI[,WZAin_gene$id[which.max(WZAin_gene$Z)]]
  env = iENV[,tg$ENV]
  
  
  ### find out if GEA is positive or negative
  GEA_sign = sign(cor(gt,env, use='pairwise.complete.obs'))
  
  
  
  ### calculate  frequency of minor allele across sampling regions
  sam_reg = cutree(hclust(dist(iENV[,c('LON','LAT')])), h=1)
  #sam_reg_maf = by(gt, sam_reg, function(x) { sum(x, na.rm=T)/(length(x[is.na(x)==F])*2) })
  #sam_reg_coord = do.call(rbind, by(meta, sam_reg, function(x) { data.frame('LON'=mean(x$LON), 'LAT'=mean(x$LAT)) }))
  #sam_reg_env = by(env, sam_reg, mean)
  
  sam_reg_maf = gt
  sam_reg_coord = data.frame('LON'=iENV$LON,
                             'LAT'=iENV$LAT)
  sam_reg_env = env
  
  # ### get rid of samples too far away
  # tk = sam_reg_coord$LON>quantile(sam_reg_coord$LON, 0.25)-IQR(sam_reg_coord$LON)*2.5&
  #   sam_reg_coord$LON<quantile(sam_reg_coord$LON, 0.75)+IQR(sam_reg_coord$LON)*2.5&is.na(gt)==F
  # 
  # sam_reg_coord = sam_reg_coord[tk,]
  # sam_reg_maf = sam_reg_maf[tk]
  # sam_reg_env = sam_reg_env[tk]
  # 
  
  
  ### Prepare colorscales for plotting
  #COLMAF = colorRampPalette(c('#4B92C7','#954489','#E74D69'))(10)[cut(sam_reg_maf, breaks = seq(-0.01,1.01, length.out=10))]
  COLMAF = c('grey90','grey70','grey30')[sam_reg_maf+1]
  COLBOX = c('grey90','grey70','grey30')
  
  if (GEA_sign==-1) { 
    #COLMAF = colorRampPalette(rev(c('#4B92C7','#954489','#E74D69')))(10)[cut(sam_reg_maf, breaks = seq(-0.01,1.01, length.out=10))] 
    COLMAF = rev(c('grey90','grey70','grey30'))[sam_reg_maf+1]
    COLBOX = rev(COLBOX)
  }
  
  
  ### Plot Map
  #pdf(file = paste0('FIGURES/GEAs/',tg$locus,'_',ds,'.pdf'), width = 8, height = 8)
  
  layout(matrix(c(1,1,1,1,
                  1,1,1,1,
                  1,1,1,1,
                  2,1,1,1), nrow=4, byrow=T))
  plotGEAgeo(coord=sam_reg_coord, df=100, col=COLMAF, env=sam_reg_env, eVar=tg$ENV, main=paste0(ds, '\nP=',signif(min(WZAin_gene$pval, na.rm=T), 3)))
  
  plotGEAbp(gt = sam_reg_maf, env=sam_reg_env, envLab = tg$ENVD)
  
  
  #dev.off()
  
}




