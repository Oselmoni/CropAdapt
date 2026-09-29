library(SNPRelate)
library(rnaturalearth)
library(terra)
library(sf)
source('SCRIPTS/custom_R_functions.R')


### Load list of datasets
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$dataset

### Get list of selected datasets
selectedDS = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, row.names=1, sep='\t')

### Load orthogroups list
load('DATA/SNP_ANNOTATION/orthogroups.rda')
rownames(orthogroups) = orthogroups$Orthogroup
### get land for plotting
land = ne_countries(scale='large')

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID
envVars$ShortName =  gsub("\\\\n", "\n", envVars$ShortName) # maintain the line break

### Set type of gea
GEATs = c('LFMM','KTAU') 


for (GEAT in GEATs) {

### Load TOPGENES
load(paste0('DATA/GEA_OUTPUT/',GEAT,'/TOPGENES.rda'))


for (i in 1:nrow(TOPGENES)) {

## Isolate gene of interest
tg = TOPGENES[i,-ncol(TOPGENES)]
tg


####
###### Show Output Plots
####


## Retrieve datasets in which gene is significant
sortedPs = tg[1,9:length(tg)][order(as.numeric(tg[1,9:length(tg)]))]

top_ds = names(sortedPs[1:tg$n_est])
top_ds = substr(top_ds, 3,nchar(top_ds)-1)

{
png(paste0('FIGURES/GEAs/topgene_',GEAT,'_',tg$locus,'_',tg$ENV,'.png'), width = 6, height = 1.5*length(top_ds), units = 'in', res = 300)
par(mfrow=c(length(top_ds),4))


### Set Colors
xCOLCS = colorRampPalette(c('#8194E9','#DD7D7D'))(10)
xCOLMAF = colorRampPalette(c('#8194E9','#DAA1DA','#DD7D7D'))(10)
xCOLBOX =c('#8194E9','#DAA1DA','#DD7D7D')


for (ds in top_ds) {

  
  ## retrieve WZA input & output
  load(paste0('DATA/GEA_OUTPUT/',GEAT,'/',ds,'/SNP_ANNOTATION_OL.rda'))
  load(paste0('DATA/GEA_OUTPUT/',GEAT,'/',ds,'/',tg$ENV,'.rda'))
  
  WZAin = SNP_ANNOTATION_OL
  WZAin$pval = WZA[[1]][WZAin$id]
  WZAout = WZA[[2]]
  rownames(WZAout) = WZAout$gene
  
  
  ## find SNPs on gene
  WZAin_gene = WZAin[WZAin$rnd_og==tg$locus,]
  

  ## transform GEA p-values to empirical p-values
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

  # use hierarchical clustering to split populations by geographical distancw
  DIST = as.dist(st_distance(st_as_sf(iENV, coords = c("LON","LAT"), crs = 4326))/200)
  sam_reg = cutree(hclust(DIST), h=50) # group sites up to 50 km apart
 
  ### Calculate mean coordnates, genotype and environmental variable by region 
  GTENV_reg = do.call(rbind, by(iENV[,c('LON','LAT')], sam_reg, function(x) { data.frame('LON'=mean(x$LON), 'LAT'=mean(x$LAT)) }))
  
  # add maf and env
  GTENV_reg$MAF = as.numeric(by(gt, sam_reg, function(x) { sum(x, na.rm=T)/(length(x[is.na(x)==F])*2) }))
  GTENV_reg$ENV = as.numeric(by(env, sam_reg, mean))
  
  #### For display purposes, remove points that are geographical outliers
  meanD = apply(as.matrix(dist(GTENV_reg[,1:2])), 1, mean)
  print(paste0('removed:', sum(meanD>quantile(meanD,0.75)+IQR(meanD)*5)))
  tk = (meanD<quantile(meanD,0.75)+IQR(meanD)*5)
  GTENV_reg = GTENV_reg[tk,]
  
  ### Set extent of plotting map area
  coord = GTENV_reg[,c('LON','LAT')]
  dX = diff(range(coord[,1]))
  dY = diff(range(coord[,2]))
  

  ### Set boundaries of plotted area, so that overall plot is a square
  if (dX<dY) {
    minY = min(coord[,2])-dY*0.05
    maxY = max(coord[,2])+dY*0.05
    
    delta = ((dY*1.1)-dX)/2
    minX = min(coord[,1])-delta
    maxX = max(coord[,1])+delta
    
  } else {
    minX = min(coord[,1])-dX*0.05
    maxX = max(coord[,1])+dX*0.05
    
    delta = ((dX*1.1)-dY)/2
    minY = min(coord[,2])-delta
    maxY = max(coord[,2])+delta
  }
  
  # Set colors
  COLCS = xCOLCS
  COLBOX = xCOLBOX
  COLMAF = xCOLMAF
  
  # If negative association, invert GEA colors
  if (GEA_sign==-1) {COLMAF=rev(COLMAF);COLBOX=rev(COLBOX)}
  
  
  ####
  ###### Make plots
  ####
  #par(mfrow=c(1,4))
  
  
  ## Plot Manhattan plot
  par(mar=c(3,3,2,1))
  manhattanPlot(p=genes_MP$EP, chr=genes_MP$CHRN, pos=genes_MP$POS, sig=which(genes_MP$OG==tg$locus), main='', chrL=gsub('chr(.*)','\\1',metaCHR$id), col='purple')
  legend('topright', '',bg=NA, box.lwd = 0,title=paste0('P=',signif(10^-genes_MP$EP[which(genes_MP$OG==tg$locus)][1],3)), cex=1, title.col = 'purple')

  pos.title = par('usr')[1]-(par('usr')[2]-par('usr')[1])*0.4
  mtext(paste0(LETTERS[which(top_ds==ds)],') ', selectedDS[ds,'name']), line=0.75, font=2, at=pos.title, cex=0.8, adj=0)
  


  
  

  ## Plot GEA
  
  # boxplot
  par(mar=c(3,3,2,1))
  plotGEAbp(gt = gt, env=env, envLab = envVars[tg$ENV,'ShortName'], COLCS = COLCS)
  
  ## Set scalebar size, based on extent of study area
  if (dX > 50 ) {scale=NA}
  if (dX < 50 & dX > 20 ) {scale= 1000}
  if (dX < 20 & dX > 10 ) {scale= 500}
  if (dX < 10 & dX > 5 ) {scale= 250}
  if (dX < 5 ) {scale=100}
  
  # plot ENV on map
  par(mar=c(1,1,2,1))
  plot(NA, xlim=c(minX,maxX), ylim=c(minY,maxY), axes=F)
  plot(land, add=T, col='grey90', border=NA)
  points(GTENV_reg$LON, GTENV_reg$LAT, col=COLCS[cut(GTENV_reg$ENV, 10)], pch=16 , add=T, border=NA, cex=0.75)
  if (ds==top_ds[1]) {title(main=paste0(' ',envVars[tg$ENV,'ShortName']), line=0.2, cex.main=0.85)}
  legend('bottomleft', legend = c(signif(min(GTENV_reg$ENV),3), '', '', signif(max(GTENV_reg$ENV),3)), col=COLCS[c(1,4,7,10)], pch=15, box.lwd = 0, bg=adjustcolor('white',0.4), pt.cex = 2, title.adj = 0, cex=0.7)
  scalebar(mean(par('usr')[1])+diff(par('usr')[1:2])*0.5, mean(par('usr')[3])+diff(par('usr')[3:4])*0.15, length_km = scale)
  box(lwd=0.5)
  
         
  # plot GT on map
  par(mar=c(1,1,2,1))
  plot(NA, xlim=c(minX,maxX), ylim=c(minY,maxY), axes=F)
  plot(land, add=T, col='grey90', border=NA)
  points(GTENV_reg$LON, GTENV_reg$LAT, col=COLMAF[cut(GTENV_reg$MAF, 10)], pch=16 , add=T, border=NA, cex=0.75)
  if (ds==top_ds[1]) {title(main='Allele Frequency', line=0.2, cex.main=0.85)}
  legend('bottomleft', legend = c(signif(min(GTENV_reg$MAF),3), '','', signif(max(GTENV_reg$MAF),3)), col=COLMAF[c(1,4,7,10)], pch=15, box.lwd = 0, bg=adjustcolor('white',0.4), pt.cex = 2, title.adj = 0, cex=0.7)
  scalebar(mean(par('usr')[1])+diff(par('usr')[1:2])*0.5, mean(par('usr')[3])+diff(par('usr')[3:4])*0.15, length_km = scale)
  box(lwd=0.5)
  
  
  

}
dev.off()
} 
}
} 


