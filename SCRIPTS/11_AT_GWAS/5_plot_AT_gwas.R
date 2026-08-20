library(vcfppR)
library(openxlsx)
library(terra)
source('SCRIPTS/custom_R_functions.R')
library(randomForest)
library(foreach)
library(doParallel)
library(ontologyIndex)
 
### load gene lists
load('DATA/ATHALIANA/GENE.rda')
GENE$geneID = gsub('ID=gene:(.*);Name.*','\\1',GENE$V9)

# load exons list
load('DATA/ATHALIANA/EXONS.rda')
EXONS$transcriptID = gsub('Parent=transcript:(.*);Name.*','\\1',EXONS$V9)
EXONS$geneID = substr(EXONS$transcriptID, start = 1, stop = nchar(EXONS$transcriptID)-2)
EXONS$transcriptV = gsub('.*\\.(.*)','\\1', EXONS$transcriptID)

### If necessary, downlaod the plant trait ontology index
if (file.exists('DATA/ATHALIANA/to.obo')==F) {
  download.file("http://purl.obolibrary.org/obo/to.obo",
                destfile = "DATA/ATHALIANA/to.obo")}


# Build ontology index of plant traits
to = get_ontology("DATA/ATHALIANA/to.obo")



### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID


### Load env data
load('DATA/ATHALIANA/meta.rda')


####
###### PREPARE PHENOTYPES
####
### load Phenotypes
load('DATA/ATHALIANA/PHENO_META.rda')

# remove "climate of origin" phenotypes
PHENO.META = PHENO.META[PHENO.META$to_name!='other miscellaneous trait',]
PHENO.META = PHENO.META[PHENO.META$name!='clim-bio1',]


### Load araGWAS results
load('DATA/ATHALIANA/araGWAS.rda')




####
###### PREPARE FUNCTIONS
####

### Match definition of variants
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='INTRON'] = 'intron_variant'
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='NON_SYNONYMOUS_CODING'] = 'missense_variant'
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='STOP_GAINED'] = 'stop_gained'
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='SYNONYMOUS_CODING'] = 'synonymous_variant'
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='UTR_3_PRIME'] = '3_prime_UTR_variant'
araGWAS$snp.annotations.0.effect[araGWAS$snp.annotations.0.effect=='UTR_5_PRIME'] = '5_prime_UTR_variant'



## define symbols for mutation types
MUT_PCH = do.call(rbind, list(
  "downstream_gene_variant" = c(22,0, 'I')         ,                                                                                     
  "upstream_gene_variant" = c(23, 90, 'I')                   ,                                                                               
  "intron_variant"  =  c(21, 0, 'I')                    ,                                                                                           
  "5_prime_UTR_variant" = c(23, 90, 'U'),
  "3_prime_UTR_variant" = c(22,0, 'U')             ,                                                                               
  "synonymous_variant" = c(21, 0, 'S')                     ,                                                                                
  "5_prime_UTR_premature_start_codon_gain_variant"   = c(23, 90, 'S') ,                                                                           
  "splice_region_variant"  = c(24, 0, 'S')              , 
  "missense_variant" =  c(24, 0, 'M')                 ,                                                                                           
  "stop_gained" = c(22, 0, 'S')))
colnames(MUT_PCH) = c('pch','ang','txt')



#### Check Geno-Pheno association for locus of interest

#gene='AT1G74960'
plot_gene = function(gene) {

## get info on mutations in gene
load(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_',gene,'.rda'))
GT=SNP_AT[[1]] # get genotypes in gene
MUT=SNP_AT[[2]] # get mutations in gene

## Subest aragwas results
gwas = araGWAS[araGWAS$snp.gene_name==gene,]
if (nrow(gwas)>0) {
gwas = do.call(rbind, by(gwas, gwas$snp.position, function(x) {
  
  return(x[which.max(x$score),])
  
}))
}
## get gene exons
s_EXONS = EXONS[EXONS$geneID==gene,]


  
########
############ PLOT
########


### Plot the gene context, with mutation types and importance of every mutation 

{
par(mar=c(3,4,1,1))

### Plot as a background the transcripts exons
transcript_versions = sort(unique(s_EXONS$transcriptV))

# plot background, set y based on how many transcripts
plot(NA, xlim=range(min(s_EXONS$V4), max(s_EXONS$V5)), ylim=c(-1,max(c(length(transcript_versions), 3))+1.5), axes=F, xlab='', ylab='', main=gene)  
for (i in sort(unique(s_EXONS$transcriptV))) {
    
    exons = s_EXONS[s_EXONS$transcriptV==i,]
    
    for (e in 1:nrow(exons)) {
      rect(exons$V4[e], as.numeric(i), exons$V5[e], as.numeric(i)+0.5, col=adjustcolor('grey',0.3), border=NA)
    }
    
    # add direction arrow and transcript version
    if (exons$V7[1]=='+') {     arrows(min(exons$V4), as.numeric(i)+0.25, max(exons$V5), as.numeric(i)+0.25, col='grey');text(min(exons$V4),as.numeric(i)+0.25, paste0('.',i), pos=2, col='grey')   }
    if (exons$V7[1]=='-') {     arrows(max(exons$V5), as.numeric(i)+0.25, min(exons$V4), as.numeric(i)+0.25, col='grey');text(max(exons$V5),as.numeric(i)+0.25, paste0('.',i), pos=4, col='grey')    }
}


### Add mutations without phenotypic association
mut_jitters = c(-0.2, -0.4, -0.6, -0.8) # vertical positions of mutations on plot
mut_y = 1

for (mut in names(MUT)) {

  pos.mut = as.numeric(gsub('.*\\:(.*)','\\1',mut))
  mut_types = MUT[[mut]]
  
  
  for (mt in mut_types[mut_types%in%rownames(MUT_PCH)])  {

    if (pos.mut%in%gwas$snp.position==F) {

    ## add vertical line to mark SNP position
    lines(c(pos.mut,pos.mut), c(mut_jitters[mut_y],par('usr')[4]), col='grey', lwd=0.25 )
    
    ## add mutation symbol
    points(pos.mut, mut_jitters[mut_y], pch=as.numeric(MUT_PCH[mt,'pch']), bg='grey80', cex=1.5, lwd=0.5) 
    text(pos.mut, mut_jitters[mut_y], MUT_PCH[mt,'txt'], cex=0.5)
    
    ## update vertical position to avoid overlaps
    mut_y=mut_y+1;mut_y[mut_y>4] = 1
    
    }
  }
  }

### Plot barplots indicating SNP-importance in  phenotype associations

## set colorscale per phenotype type
set.seed(0);traits = unique(PHENO.META[as.character(gwas$study.id),'to_name'])
traits.cs = rainbow(length(traits));names(traits.cs)=traits


if (nrow(gwas)>0) {
for (pos in unique(gwas$snp.position)) {

  s_gwas = gwas[gwas$snp.position==pos,] # subset gwas results for SNP
  s_gwas = s_gwas[which.max(s_gwas$score),] # find phenotype with top score
  pto_trait =   PHENO.META[as.character(s_gwas$study.id),'to_name']
  
  colSNP = traits.cs[pto_trait] # find color code for phenotype type
  
  

  ## add vertical line to mark SNP position
  lines(c(pos,pos), c(mut_jitters[mut_y],par('usr')[4]), col='grey', lwd=0.25 )
  
  ## add mutation symbol
  points(pos, mut_jitters[mut_y], pch=as.numeric(MUT_PCH[s_gwas$snp.annotations.0.effect,'pch']), bg=adjustcolor(colSNP, 0.5), cex=1.5, lwd=0.5) 
  text(pos, mut_jitters[mut_y], MUT_PCH[s_gwas$snp.annotations.0.effect,'txt'], cex=0.5)
  
  ## update vertical position to avoid overlaps
  mut_y=mut_y+1;mut_y[mut_y>4] = 1
  
  
  
  # find score of top phenotype for a given snp
  score = s_gwas$score
  
 
  # set relative bar height
  r.imp = score/max(gwas$score)
  bar_height = r.imp*round(par('usr')[4]) # put relative importance in the scale of the plot
  
  
  ## add importance of every SNP
  rect(pos-10, 0, pos+10, bar_height, border=NA, col=adjustcolor(colSNP,0.5))
  
}
}

### Add genome line
lines(c(min(s_EXONS$V4), max(s_EXONS$V5)), c(0,0), lwd=3)

### Add labels on axes
axis(1, lwd=0, tick = T, lwd.ticks = 1)
if (nrow(gwas)>0) {axis(2, at=seq(0,par('usr')[4], length.out=3), labels = signif(seq(0, max(gwas$score), length.out=3),2 ), las=1); title(line=3, ylab='                      -log(P)')}
title(xlab=paste0('Chr ',s_EXONS$V1[1]), line=2)

}


}

## function to get a legend for a list of genes
getLegend = function(gene) {
  
  
## get info on mutations in gene
load(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_',gene,'.rda'))
MUT=SNP_AT[[2]] # get mutations in gene

## Subest aragwas results
gwas = araGWAS[araGWAS$snp.gene_name==gene,]
if (nrow(gwas)>0) {
gwas = do.call(rbind, by(gwas, gwas$snp.position, function(x) {
  
  return(x[which.max(x$score),])
  
}))}


# SUBSET legend of mutation types
mut_types_in = unique(c(gwas$snp.annotations.0.effect, unlist(MUT)))
mut_legend = MUT_PCH[rownames(MUT_PCH)%in%mut_types_in,]

## set colorscale per phenotype type
set.seed(0);traits = unique(PHENO.META[as.character(gwas$study.id),'to_name'])
traits.cs = rainbow(length(traits));names(traits.cs)=traits
# 

par(mar=c(0,0,0,0))
plot(1,1,col=0, axes=F, xlab='', ylab='')
if (length(traits.cs)>0) {
legend('left', legend = traits, lty=1, col=traits.cs, lwd=2, title='Trait associated', box.lwd=0, bg=NA, cex=,1, title.adj = 0, title.font = 2)}
par(mar=c(0,0,0,0))
plot(1,1,col=0, axes=F, xlab='', ylab='')
if (nrow(mut_legend)>0) {
legend('left', legend = rownames(mut_legend), pch =as.numeric(mut_legend[,1]), bg = NA, pt.cex=1.5, box.lwd = 0, title = 'Mutation type', cex=1, title.adj = 0, title.font = 2)
legend('left', legend = rownames(mut_legend), pch =mut_legend[,3], bg = NA, pt.cex=0.5, box.lwd=0, title='Mutation type', cex=1, title.adj = 0, title.font = 2)
}

}



### Plot genes for every orthogroup

for (geat in c('LFMM','KTAU')) {

  load(paste0('DATA/GEA_OUTPUT/',geat,'/TOPGENES.rda'))
  
for (og in unique(TOPGENES$locus)) {
  
  tgene = TOPGENES[TOPGENES$locus==og,]
  tgene = tgene[which.min(tgene$pooled_q),] # keep only best association per gene

  ## find all a-thaliana orthologs
  ATs = unique(gsub('(.*)\\..*','\\1', unique(unlist(strsplit(tgene$OG,', ')))))
  

  for (at in ATs) {
    
    ## plot every gene
    png(filename = paste0('FIGURES/athaliana/',geat,'_',at,'.png'), width = 7, height = 3, res = 300 ,units = 'in')
    
    layout(matrix(c(1,1,2,
                    1,1,3), ncol=3, byrow=T))
    
    plot_gene(gene=at)
    getLegend(gene = at)
    
    dev.off()
  }
  }
}


