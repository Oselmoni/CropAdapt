### Figure to summarize annotation of GO and phenotype for every topgene

### Load list of top candidate genes from LFMM picmin
load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
TOPGENES.LFMM = TOPGENES
TOPGENES.LFMM$GEA = 'LFMM'
### Load list of top candidate genes from KTAU picmin
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
TOPGENES.KTAU = TOPGENES
TOPGENES.KTAU$GEA = 'KTAU'

### Create unique table of top genes
TOPGENES = rbind(TOPGENES.LFMM, TOPGENES.KTAU)

### Keep only best env per top gene
TOPGENES = do.call(rbind, by(TOPGENES, TOPGENES$locus, function(x) {return(x[which.min(x$pooled_q),])}))

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID
envVars$ShortName =  gsub("\\\\n", "\n", envVars$ShortName) # maintain the line break


### load GO annotations
load('DATA/ATHALIANA/GOtab.rda')

### Create manually curated groups of ontologies for plotting, assign clusters from GOtab
GO_sem = list('Response to amino acid'=c(1), 
           'Gibberellin biosynthesis'=c(2),
           'Response to water deprivation / salt stress'=c(5),
           'Calcium ion homoestasis / transport / signaling'=c(3,4,6),
           'Response to cold'=c(7),
           'Response to light'=c(8))


### Add gene names without GO annotation
AT = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

GO.tab[,AT[AT%in%colnames(GO.tab)==F]] = as.character(NA)


### for every topGenes, check if they have GO annotation 
TOPGENES.GO = data.frame()

### 
for (og in unique(TOPGENES$locus)) {
  

  ## subset topgenes, find best env
  TOPGENES_og = TOPGENES[TOPGENES$locus==og,]
  TOPGENES_og = TOPGENES_og[which.min(TOPGENES_og$pooled_q),]
  
  ## find top env
  TOPGENES.GO[og,'env'] = TOPGENES_og$ENVD
  
  ## find at orthologs
  at = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES_og$OG, ', '))))
  
  if (length(at)==0) {
    TOPGENES.GO[og,names(GO_sem)] = ''
  } else {
    
      for (sem in names(GO_sem)) {
          
        annotations = GO.tab[GO.tab$group==GO_sem[[sem]],at] # get annotations for semantic group
        
        if (sum(annotations=='x', na.rm=T)>0) { # if at least one hit: mark as present
              TOPGENES.GO[og,sem] = 'x'
        } else {
          TOPGENES.GO[og,sem] = ''
        }
        
      }
    
  }
  
}



### load phenotypes
load('DATA/ATHALIANA/PROTANN.rda')
PROT.ANN


### for every topGenes, check if they have GO annotation 
TOPGENES.PHENO = data.frame()

for (og in unique(TOPGENES$locus)) {
  
  tg = TOPGENES[TOPGENES$locus==og,]
  
  ## get phenotype with strongest association with gene (if multiple paralogs, pick the most significant)
  PROT.ANN_og = PROT.ANN[PROT.ANN$OrthogeneID==og,]
  
  PHs = PROT.ANN_og$TO_name[which.max(as.numeric(PROT.ANN_og$top_pheno_score))] 
  
  # get info if phenotype made it through multiple testing
  Bonferroni = PROT.ANN_og$top_pheno_bonferroni[which.max(as.numeric(PROT.ANN_og$top_pheno_score))]
  
  # score
  if (length(Bonferroni)==1) {
  if (Bonferroni=='False') { TOPGENES.PHENO[og,PHs] = 'f'}
  if (Bonferroni=='True') { TOPGENES.PHENO[og,PHs] = 't'}
  }
  
}

## Sort by env and by q
TOPGENES.GO = TOPGENES.GO[order(TOPGENES.GO$env, TOPGENES$pooled_q),]
TOPGENES.PHENO = TOPGENES.PHENO[rownames(TOPGENES.GO),]
TOPGENES.PHENO[is.na(TOPGENES.PHENO)] = ''

## set order of phenos
TOPGENES.PHENO = TOPGENES.PHENO[c('days to flowering trait', 'yield trait',
                 'root morphology trait','root branching','relative root length',
                 'bacterial disease resistance','anthocyanin content','zinc concentration')]

## Set order of GO terms
GO.COUNTS = sort(apply(TOPGENES.GO[,-1], 2,function(x) {sum(x=='x', na.rm=T)}), decreasing = T)

TOPGENES.GO = TOPGENES.GO[,c('env',names(GO.COUNTS))]

#### Create summary figure


pdf('FIGURES/athaliana/GO_PHENO.pdf', width = 5, height = 6)
{
layout(matrix(c(2,2,2,4,4,
                    1,1,1,3,3,
                    1,1,1,3,3,
                    1,1,1,3,3), nrow=4, byrow=T))


## Plot goterms by gene
par(mar=c(10,12,1,1))
plot(NA, xlim=c(2,ncol(TOPGENES.GO)), ylim=c(-nrow(TOPGENES.GO),-1), axes=F, xlab='', ylab='')
abline(v=2:ncol(TOPGENES.GO), col='pink')
abline(h=-1:-nrow(TOPGENES.GO), col='pink')
for (i in 1:nrow(TOPGENES.GO)) {
  for (l in 2:(ncol(TOPGENES.GO))) {
      if (TOPGENES.GO[i,l]=='x') {points(l,-i, col='pink', pch=16, cex=2)}
  }
}

axis(1, at = 2:ncol(TOPGENES.GO), labels = rep('', ncol(TOPGENES.GO)-1), tick =T)
text(x = 2:ncol(TOPGENES.GO), y = par("usr")[3]-1,  # Position labels slightly below axis
     labels = colnames(TOPGENES.GO)[2:ncol(TOPGENES.GO)], srt = 45, xpd = TRUE, pos=2, cex=0.75)

axis(2, at = -1:-nrow(TOPGENES.GO), labels = paste0(rownames(TOPGENES.GO)), las=2)


## Plot goterms counts
par(mar=c(0,12,1,1))
barplot(GO.COUNTS, width=0.5, space=1, border=NA, col='pink', axes=F, names='', xlim=c(1,length(GO.COUNTS))-0.25)
axis(2, las=2)


## Plot phenotypes by gene
par(mar=c(10,0,1,3))
plot(NA, xlim=c(1,ncol(TOPGENES.PHENO)), ylim=c(-nrow(TOPGENES.PHENO),-1), axes=F, xlab='', ylab='')
abline(v=1:ncol(TOPGENES.PHENO), col='lightblue')
abline(h=-1:-nrow(TOPGENES.PHENO), col='lightblue')
for (i in 1:nrow(TOPGENES.PHENO)) {
  for (l in 1:(ncol(TOPGENES.PHENO))) {
    if (TOPGENES.PHENO[i,l]=='f') {points(l,-i, col=adjustcolor('lightblue',.5), pch=16, cex=2)}
    if (TOPGENES.PHENO[i,l]=='t') {points(l,-i, col='lightblue', pch=16, cex=2)}
    
  }
}
axis(1, at = 1:ncol(TOPGENES.PHENO), labels = rep('', ncol(TOPGENES.PHENO)), tick =T)
text(x = 1:ncol(TOPGENES.PHENO), y = par("usr")[3]-1,  # Position labels slightly below axis
     labels = colnames(TOPGENES.PHENO), srt = 45, xpd = TRUE, pos=2, cex=0.75)



## Plot phenotype counts
COUNTS = apply(TOPGENES.PHENO,2,function(x) {sum(x!='', na.rm=T)})
COUNTS_sig = apply(TOPGENES.PHENO,2,function(x) {sum(x=='t', na.rm=T)})
par(mar=c(0,0,1,3))
barplot(COUNTS, width=0.5, space=1, names='',axes=F, border=NA, col=adjustcolor('lightblue',0.5), xlim=c(1,length(COUNTS))-0.25)
barplot(COUNTS_sig, width=0.5, space=1, names='',axes=F, border=NA, col=adjustcolor('lightblue',1), xlim=c(1,length(COUNTS))-0.25, add=T)
axis(4, las=2)
}
dev.off()

