library(qvalue)
library(Rmisc)

### Load dataset info
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, sep='\t')
rownames(dslist) = dslist$dataset

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars) = envVars$VariableID


### Create table to store p-vals (kTAU)
TAB.KTAU = data.frame() # to store fraction
TAB.KTAU.abs = data.frame() # to store abs number

for (ds in dslist$dataset) {
for (env in envVars$VariableID) {
      
      
    load(paste0('DATA/GEA_OUTPUT/KTAU/',ds,'/',env,'.rda'))

    WZAout = WZA[[2]]
    
    QS = qvalue(WZAout$Z_pVal)$qvalue
    
    TAB.KTAU[ds,env] = mean(QS<0.05)
    TAB.KTAU.abs[ds,env] = sum(QS<0.05)
    
}
}
 

### Create table to store p-vals (LFMM)
TAB.LFMM = data.frame()
TAB.LFMM.abs = data.frame() # to store abs number


for (ds in dslist$dataset) {
  for (env in envVars$VariableID) {
    
    
    load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/',env,'.rda'))
    
    WZAout = WZA[[2]]
    
    QS = qvalue(WZAout$Z_pVal)$qvalue
    
    TAB.LFMM[ds,env] = mean(QS<0.05)
    TAB.LFMM.abs[ds,env] = sum(QS<0.05)
    
    
  }
} 

max(TAB.KTAU)*100

### Function to transform table to figure
tab_plot = function(relative, absolute, maxval) {


  ## set colorscale
  CCC = adjustcolor(colorRampPalette(c('white','lightblue','pink'))(10), 0.8)

  
  ### Set the composed plot layout
  layout(matrix(c(3,3,3,3,2,
                  1,1,1,1,4,
                  1,1,1,1,4,
                  1,1,1,1,4), byrow=T, ncol=5))
  
 
  ## plot a canvas
  par(mar=c(18,15,0,0))
  plot(NA, xlim=c(1,ncol(relative)), ylim=c(1,nrow(relative)), axes=F, xlab='', ylab='')

  ### add grid on bg
  abline(v=1:ncol(relative), col='grey90')
  abline(h=1:nrow(relative), col='grey90')
  
  ### store x-y range
  xyrange = par('usr')
  
  # add value for every cell
  for (x in 1:ncol(relative)) {
    for (y in 1:nrow(relative)) {
      
      points(x,y, bg=CCC[cut(relative[y,x], breaks = seq(-0.000001,maxval,length.out=10))], pch=22, col=NA, cex=3)
      
      if (absolute[y,x]>0) {text(x,y, label=absolute[y,x], cex=0.75)}
      
      
    }
  }
  
  # add labels
  axis(2, 1:nrow(relative), labels=dslist[rownames(relative),'name'], las=2, tick = F)  
  axis(1, 1:ncol(relative), labels=envVars[colnames(relative),'Description'], las=2, tick = F)  

  
  ## add legend
  par(mar=c(0,0,0,0))
  plot(1,1, col=0, axes=F, xlab='', ylab='') # canvas for legend
  legend('bottomleft', pt.bg=CCC, pt.cex = 3, col=0, legend = signif(seq(-0.000001,maxval,length.out=10),3), pch=22, box.lwd = 0, title='% of significant \ngenes (q<0.05)', cex=0.8)
  
  ## Add historgram for env
  par(mar=c(1,15,1,0))
  barplot(apply(relative, 2, mean), width=0.5, space=1, axes=T, names='', las=2, border=0, xlim=c(1,43)-0.25)

 
  
  ## Add historgram for ds
  par(mar=c(18,0,0,1))
  barplot(apply(relative, 1, mean),width=0.5, space=1,horiz=T, axes=T, names='', las=2, border=0,  ylim=c(1,14)-0.25)

}



### Create plots DS by ENV 
png(filename = 'FIGURES/GEAs/KTAU_GEA.png', width = 9, height = 5.5, units = 'in', res = 300)
tab_plot(relative = TAB.KTAU*100, absolute = TAB.KTAU.abs, maxval = 0.56)
dev.off()

png(filename = 'FIGURES/GEAs/LFMM_GEA.png', width = 9, height = 5.5, units = 'in', res = 300)
tab_plot(relative = TAB.LFMM*100, absolute = TAB.LFMM.abs, maxval = 0.56)
dev.off()


### Create plot comparing landraces vs. method
TAB.KTAU.WR = TAB.KTAU[lapply(strsplit(rownames(TAB.KTAU), '_'), function(x) {return(x[2])})=='wr',]
TAB.KTAU.LA = TAB.KTAU[lapply(strsplit(rownames(TAB.KTAU), '_'), function(x) {return(x[2])})=='la',]
TAB.LFMM.WR = TAB.LFMM[lapply(strsplit(rownames(TAB.LFMM), '_'), function(x) {return(x[2])})=='wr',]
TAB.LFMM.LA = TAB.LFMM[lapply(strsplit(rownames(TAB.LFMM), '_'), function(x) {return(x[2])})=='la',]

png(filename = 'FIGURES/GEAs/GEA_types.png', width = 7, height = 5, units = 'in', res = 300)
par(mar=c(10,5,1,1))
boxplot( unlist(as.vector(TAB.KTAU.WR))*100,
         unlist(as.vector(TAB.KTAU.LA))*100,
         unlist(as.vector(TAB.LFMM.WR))*100,
         unlist(as.vector(TAB.LFMM.LA))*100,
         names=c('Kendall Τ, CWR', 'Kendall T, landraces', 'LFMM, CWR', 'LFMM, landraces'), las=2)
         title(ylab ='% of significant genes (q<0.05)',line=4)
dev.off()        
 

### Get values for reporting

# difference between tau and lfmm
CI(as.matrix(TAB.KTAU))*100
CI(as.matrix(TAB.LFMM))*100

# difference between variable
sort(apply(rbind(TAB.KTAU, TAB.LFMM), 2, mean))*100
CI(rbind(TAB.KTAU, TAB.LFMM)[,'DRYDA'])*100
CI(rbind(TAB.KTAU, TAB.LFMM)[,'ARIDITY'])*100
CI(rbind(TAB.KTAU, TAB.LFMM)[,'OCSTO'])*100

# difference crop types
CI(as.matrix(rbind(TAB.KTAU.LA, TAB.LFMM.LA)))*100
CI(as.matrix(rbind(TAB.KTAU.WR, TAB.LFMM.WR)))*100

# differences between datasets
sort(apply(cbind(TAB.KTAU, TAB.LFMM), 1, mean))*100
CI(as.matrix(cbind(TAB.KTAU, TAB.LFMM)['corn_wr',]))*100
CI(as.matrix(cbind(TAB.KTAU, TAB.LFMM)['tomato_wr',]))*100








#### Evaluate number of orthogroups shared between datasets





# LOAD AND PREPARE ANNOTATION DATA
load('DATA/SNP_ANNOTATION/orthogroups.rda')
annotation_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$annotation_file
protein_files = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$protein_file
names(protein_files) = names(annotation_files) = dslist$dataset


# create a tab to store number of genes per dataset
OG.TAB = data.frame() 
OG = list() # create list to store orthogene from every dataset


# for every dataset: count genes with SNP & orthogenes with SNPs
for (ds in dslist$dataset) {

  ### Load imputed GT matrix
  gti_files=list.files('DATA/GEA_INPUT/GTI/', pattern = paste0('GTI_',ds,'_K'),full.names = T)
  load(gti_files[regexpr(ds, gti_files)!=-1])
  
  # load SNP to gene table
  load(paste0('DATA/SNP_ANNOTATION/',annotation_files[ds],'.rda'))
  
  # keep only SNPs used in GEA, to speed up
  SNP_ANNOTATION = SNP_ANNOTATION[which(SNP_ANNOTATION$snp%in%colnames(GTI)),]
  
  # load SNP to orthogene table
  load(paste0('DATA/GEA_OUTPUT/LFMM/',ds,'/SNP_ANNOTATION_OL.rda'))
  
  # find genes with SNPs
  GENE_w_SNP = unique(SNP_ANNOTATION$geneID)
  
  # find orthogene with SNPs 
  OGENE_w_SNP = unique(SNP_ANNOTATION_OL$rnd_og)
  
  
  # add counts table
  OG.TAB[ds,'Genes'] = length(GENE_w_SNP)
  OG.TAB[ds,'OrthoGenes'] = length(OGENE_w_SNP)
  OG.TAB[ds,'SNPs'] = nrow(SNP_ANNOTATION_OL)
  
  # keep list of OrthoGenes
  OG[[ds]] = OGENE_w_SNP  
}

### Now check how many orthogroups are shared between pairs of datasets
# get pairs
PAIRS = combn(dslist$dataset, 2)

for (i in 1:ncol(PAIRS)) {
  
  ds1 = PAIRS[1,i]
  ds2 = PAIRS[2,i]
  
  OG1 = OG[[ds1]]
  OG2 = OG[[ds2]]
  
  shared =  sum(OG1%in%OG2)
  
  OG.TAB[ds1,ds2] = shared
}

rownames(OG.TAB) = dslist[rownames(OG.TAB),'name']
colnames(OG.TAB)[4:16] = dslist[colnames(OG.TAB)[3:15],'name']

write.csv(OG.TAB, 'FIGURES/GEAs/OG_per_dataset.csv', row.names = T, col.names  = T, quote=F)




