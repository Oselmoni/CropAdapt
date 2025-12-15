library(SNPRelate)
library(rnaturalearth)
source('SCRIPTS/custom_R_functions.R')

lf = list.files('SCRIPTS/6_MATCH_GEN_GEO_ENV/')[-1]

land=ne_countries()

## create dataframe on datasets stats
DS_STATS=data.frame()

for (f in lf) {

  print(f)
  d=substr(f, 1 ,nchar(f)-2)
  
  ### load vcf and metadata
  source(paste0('SCRIPTS/6_MATCH_GEN_GEO_ENV/',f))
  
  ### enter dataset stats
  DS_STATS[d,'N'] = nrow(meta)
  DS_STATS[d,'N_snps'] = nrow(GT)
  
  ### check missing coord
  meta = meta[is.na(meta$LON)==F&is.na(meta$LAT)==F,]
  GT=GT[,rownames(meta)]
  
  DS_STATS[d,'N_geo'] = nrow(meta)
  
  ### Create GDS file 
  tmp=tempfile() # create temporary file
  
  ### create gds object 
  snpgdsCreateGeno(tmp, GT, sample.id = colnames(GT), snp.id = rownames(GT), snp.chromosome= unlist(lapply(strsplit(rownames(GT), ':'), function(x) {return(x[1])})), snp.position = unlist(lapply(strsplit(rownames(GT), ':'), function(x) {return(x[1])})))
  
  ### open gds object
  SNPS <- snpgdsOpen(tmp)
  
  
  
  ### Filter missingness by SNP
  MN_snp = snpgdsSNPRateFreq(SNPS)$MissingRate
  snps_mn = snpgdsSelectSNP(gdsobj = SNPS,    maf = 0, missing.rate = 0.1, autosome.only=F, remove.monosnp=T)
  
  ## Filter missigness by IND
  MN_ind = snpgdsSampMissRate(SNPS, snp.id=snps_mn)
  ind_mn = colnames(GT)[which(MN_ind<0.1)]
  
  ## Filter MAF 
  MAF_snp = snpgdsSNPRateFreq(SNPS, sample.id = ind_mn , snp.id=snps_mn)$MinorFreq
  snps_mn_maf = snpgdsSelectSNP(gdsobj = SNPS, maf = 0.05, sample.id = ind_mn , snp.id=snps_mn, autosome.only=F, remove.monosnp=T)
  
  
  ### enter dataset stats
  DS_STATS[d,'N_AF_mn'] = length(ind_mn)
  DS_STATS[d,'N_snps_AF_mn'] = length(snps_mn)
  DS_STATS[d,'N_snps_AF_mn_maf'] = length(snps_mn_maf)
  
  
  ### Calculate PCA 
  PCA = snpgdsPCA(SNPS, sample.id=ind_mn, snp.id=snps_mn_maf, autosome.only=F)
  
  ### calculate % var explained
  PVE = PCA$eigenval/sum(PCA$eigenval, na.rm=T)*100
  
  ## Create colorscale for first two PCOA axes
  col = rgb(red =   (scale(PCA$eigenvect[,3], center=min(PCA$eigenvect[,3]), scale=diff(range(PCA$eigenvect[,3])))/1)+0, 
            blue =  (scale(PCA$eigenvect[,2], center=min(PCA$eigenvect[,2]), scale=diff(range(PCA$eigenvect[,2])))/2)+0,
            green = (scale(PCA$eigenvect[,1], center=min(PCA$eigenvect[,1]), scale=diff(range(PCA$eigenvect[,1])))/1)+0,
            alpha =1)
  
  {
    pdf(paste0('FIGURES/ds_genomic_summary/',d,'.pdf'), h=7, w=6)
    ## plot summary figure
    layout(matrix(c(1,2,3,
                    4,5,6,
                    7,8,9,
                    10,10,10,
                    10,10,10), nrow=5, byrow = T))
    
    ## plot stats
    par(mar=c(3,3,2,1))
    hist(MN_ind, breaks=100, main='\nA-Miss. rate by ind.', xlim=c(0,1), las=2)
    rect(0.1, 0, 1, par('usr')[4], border=NA, col=adjustcolor('red',0.1))
    hist(MN_snp, breaks=100, main=paste0(d,'\nB-Miss. rate by SNP'), xlim=c(0,1), las=2)
    rect(0.1, 0, 1, par('usr')[4], border=NA, col=adjustcolor('red',0.1))
    hist(MAF_snp, breaks=100, main='\nC-Minor Allele Freq.', xlim=c(0,1), las=2)
    rect(0, 0, 0.05, par('usr')[4], border=NA, col=adjustcolor('red',0.1))
    
    ## plot PCA
    par(mar=c(3,3,1,1))
    cplot(1:length(PVE[is.na(PVE)==F]), PVE[is.na(PVE)==F], main='D - Pr. Comp. Analysis', ylab='PVE', xlab='PC#')
    cplot(PCA$eigenvect[,1], PCA$eigenvect[,2], col=col, xlab=paste0('PC1 (PVE=',signif(PVE[1],2),'%)'), ylab=paste0('PC2 (PVE=',signif(PVE[2],2),'%)'))
    cplot(PCA$eigenvect[,3], PCA$eigenvect[,4], col=col, xlab=paste0('PC3 (PVE=',signif(PVE[3],2),'%)'), ylab=paste0('PC4 (PVE=',signif(PVE[4],2),'%)'))
    cplot(PCA$eigenvect[,5], PCA$eigenvect[,6], col=col, xlab=paste0('PC5 (PVE=',signif(PVE[5],2),'%)'), ylab=paste0('PC6 (PVE=',signif(PVE[6],2),'%)'))
    cplot(PCA$eigenvect[,7], PCA$eigenvect[,8], col=col, xlab=paste0('PC7 (PVE=',signif(PVE[7],2),'%)'), ylab=paste0('PC8 (PVE=',signif(PVE[8],2),'%)'))
    cplot(PCA$eigenvect[,9], PCA$eigenvect[,10], col=col, xlab=paste0('PC9 (PVE=',signif(PVE[9],2),'%)'), ylab=paste0('PC10 (PVE=',signif(PVE[10],2),'%)'))
    
    
    ## Plot map
    par(mar=c(1,1,0,1))
    plot(meta[ind_mn,c('LON','LAT')], col=NA, axes=F)
    plot(land, add=T, col='grey90', border=NA)
    points(meta[ind_mn,c('LON','LAT')], col=col, pch=16)
    box()
    dev.off()
  }
  
  #### save filter SNP and metadata table
  snpgdsCreateGenoSet(src.fn=tmp, dest.fn=paste0('DATA/GEA_INPUT/gds_',d,'.gds'), snp.id=snps_mn_maf, sample.id=ind_mn)
  
  meta = meta[ind_mn,]
  save(meta, file=paste0('DATA/GEA_INPUT/meta_',d,'.rda'))
  
  closefn.gds(SNPS)
  gc()
}

write.table(DS_STATS, 'FIGURES/ds_genomic_summary/DS_STATS.txt', col.names=T, row.names=T, quote=F, sep='\t')
