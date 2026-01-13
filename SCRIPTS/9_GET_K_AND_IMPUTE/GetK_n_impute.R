library(SNPRelate)
library(LEA)
library(rnaturalearth)
library(terra)
### Get list of name of datasets of interest
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt')[,1]

### get country boundaries
land=ne_countries(scale='large')


## set parameters for SNMF
NCPU = 7
NREP = 5
NK = 1:20

## create empty file to store imputation accuracy results by dataset
write.table(data.frame(), file = 'DATA/GEA_INPUT/GTI/IMPACC.txt', row.names = F, col.names=F, quote = F)

### for every dataset...
for (ds in ds_list[-1]) {

    # ds=ds_list[3]
    print(ds)
  
    ### Load gds
    SNPS = snpgdsOpen(paste0('DATA/GEA_INPUT/GDS/gds_',ds,'.gds'), readonly = T, allow.duplicate = T)

    ### Get genotype matrix in .geno format
    GENO=t(snpgdsGetGeno(SNPS))
    GENO[is.na(GENO)]=9

    #####
    ####### FIND OPTIMAL K
    #####

    ### Randomly sample 1000 SNP from geno matrix (to speed up calculation)
    set.seed(0);sGENO = GENO[,sample(1:ncol(GENO), size = 1000)]
    tmpgeno = tempfile(fileext = '.geno')
    write.geno(sGENO, tmpgeno)

        
    ## Run snfm on sampled geno
    set.seed(0);project = snmf(tmpgeno,
                   K = NK,
                   entropy = TRUE,
                   repetitions = NREP, CPU=NCPU, project='new')
    
    ## Get snmf out
    snmf_out = do.call(rbind,lapply(project@runs, function(x) {
      return(data.frame('CE'=x@crossEntropy,'r'=x@run,'K'=x@K))
    }))
    
    

    ## Get coefficients optimal K: lower k after which there is not decrease
    
    # estimate K+1 change for K=1
    topK=1
    TT = t.test(CE ~ as.factor(K), data = snmf_out[snmf_out$K%in%c(topK+c(0,1)),], alternative='greater')
    PVAL = TT$p.value
    
    while(PVAL < 0.01) { ## keep incresing K, until there is no differnece between K and K+1
      topK=topK+1
      TT = t.test(CE ~ as.factor(K), data = snmf_out[snmf_out$K%in%c(topK+c(0,1)),], alternative='greater')
      PVAL = TT$p.value
    }
    

    #####
    ####### IMPUTE MISSING GT USING OPTIMAL K
    #####
    
    
    ## Now re-run SNMF with optimal K for full genotype matrix
    write.geno(GENO, tmpgeno)
    set.seed(0);project = snmf(tmpgeno,
                               K = topK,
                               entropy = TRUE,
                               repetitions = 1, CPU=NCPU, project='new')
    
    
    ### Mask random GTs to evaluate imputation accuracy
    set.seed(0);sampledGTidx = sample(1:ncell(GENO), 1000, F)
    sampledGT = GENO[sampledGTidx]
    sampledGT[sampledGT==9] = NA
    
    # keep only sampled GT that are non-missing
    nm = is.na(sampledGT)==F
    sampledGT = sampledGT[nm]
    sampledGTidx = sampledGTidx[nm]
    
    # mask sampled GT
    GENO_masked = GENO;GENO_masked[sampledGTidx] = NA
    write.geno(GENO_masked, tmpgeno)
    
    ## Run imputation
    impute(project, tmpgeno,
           method = 'mode', K = topK, run = 1)


    ## Load imputed Genotype matrix
    GTI=as.matrix(read.table(paste0( substr(tmpgeno, 1, nchar(tmpgeno)-5),'.lfmm_imputed.lfmm')))
    colnames(GTI) = read.gdsn(index.gdsn(SNPS, "snp.id"))
    rownames(GTI) = read.gdsn(index.gdsn(SNPS, "sample.id"))
    
    ## Check imputation peformance
    impGT = GTI[sampledGTidx]
    IMP.ACC = mean(sampledGT==impGT, na.rm=T) # imputation accuracy
    
    write.table(data.frame(ds, IMP.ACC), file = 'DATA/GEA_INPUT/GTI/IMPACC.txt', append = T, row.names = F, col.names=F, quote = F)
  
    ## restore real GT in imputed matrix
    GTI[sampledGTidx] = sampledGT

    ## save imputed Genotype matrix
    save(GTI, file=paste0('DATA/GEA_INPUT/GTI/GTI_',ds,'_K',topK,'.rda'))
 
    ## load metadata
    load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))

    ## Get admixture coefficients top run
    ADCO = Q(project, topK, 1)

    ## get ancestral populations for every sample
    ANC.POP = rep(NA, length=nrow(meta))
    for (i in 1:ncol(ADCO)) {
        ANC.POP[ADCO[,i]>0.5]=i
    }

    ## Plot K decision
    pdf(paste0('FIGURES/snmf_k/',ds,'.pdf'), h=3, w=5)
    layout(matrix(c(1,2,2), nrow=1))
    ## Plot cross-entropy by number of ancestral populations
    boxplot(snmf_out$CE~snmf_out$K, xlab='K', ylab='cross-entropy', main=ds)
    
    plot(meta[,2:3], pch=16, col=NA, main=paste0('K=',topK))
    plot(land, add=T, col='grey90', border='grey90')
    points(meta[,2:3], pch=16, col=ANC.POP)
    box()
    dev.off()    
}    
