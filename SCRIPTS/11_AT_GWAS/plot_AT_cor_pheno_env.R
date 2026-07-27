library(vcfppR)
library(openxlsx)
library(terra)
source('SCRIPTS/custom_R_functions.R')
library(randomForest)
library(foreach)
library(doParallel)

### load gene lists
load('DATA/ATHALIANA/GENE.rda')
GENE$geneID = gsub('ID=gene:(.*);Name.*','\\1',GENE$V9)
EXONS = GFF[GFF$V3=='exon',]
EXONS$transcriptID = gsub('Parent=transcript:(.*);Name.*','\\1',EXONS$V9)
EXONS$geneID = substr(EXONS$transcriptID, start = 1, stop = nchar(EXONS$transcriptID)-2)
EXONS$transcriptV = substr(EXONS$transcriptID, start = nchar(EXONS$transcriptID), stop = nchar(EXONS$transcriptID))

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID


### Load env data
load('DATA/ATHALIANA/meta.rda')



####
###### PREPARE POPULATION STRUCTURE
####

## load neutral SNP matrix 
load('DATA/ATHALIANA/SNPS.rda')

## calculate neutral population structure using a PCoA
PCoA = cmdscale(dist(SNPS), k = 20, eig = T)
PS = PCoA$points[,1:9]
colnames(PS) = paste0('PS',1:ncol(PS))




####
###### PREPARE PHENOTYPES
####
### load Phenotypes
load('DATA/ATHALIANA/PHENO.rda')
load('DATA/ATHALIANA/PHENO_META.rda')

# remove "climate of origin" phenotypes
PHENO.META = PHENO.META[PHENO.META$to_name!='other miscellaneous trait',]
PHENO.META = PHENO.META[PHENO.META$name!='clim-bio1',]
PHENO = PHENO[,rownames(PHENO.META)]

## Keep only accessions that are genotyped
PHENO = PHENO[rownames(PHENO)%in%rownames(SNPS),]

# Check how many sample per phenotype
hist(apply(PHENO, 2, function(x) {sum(is.na(x)==F)}), breaks=100)
quantile(apply(PHENO, 2, function(x) {sum(is.na(x)==F)}), breaks=100)

## Keep only phenotype with at least 100 N
PHENO = PHENO[,apply(PHENO, 2, function(x) {sum(is.na(x)==F)})>50]

## Keep only phenotype with at least 10 different values (avoid categorical phenotype)
PHENO = PHENO[,apply(PHENO, 2, function(x) {length(unique(x))})>10]




####
###### PREPARE FUNCTIONS
####

# function to extract SNPs from a given gene
getGENO = function(gene, MN.T=0.2, MAF.T=0.05) {
  
  ### Find start stop of gene
  s_GENE = GENE[GENE$geneID==gene,]
  
  chr = s_GENE$V1
  sta = s_GENE$V4
  end = s_GENE$V5
  
  if (chr %in% as.character(1:5)==F) { return() } 
  
  
  ### Extract SNPs
  vcf = vcftable(paste0('DATA/ATHALIANA/AT_chr',chr,'_eff.vcf.gz'), region = paste0(chr,':',sta,'-',end), vartype = 'snps')
  
  ### Extract GT matrix, add id for SNPs and samples
  GT = t(vcf$gt)
  
  if (ncol(GT)<=1|nrow(GT)<=1) {return()}
  rownames(GT) = vcf$samples
  GT = GT[rownames(meta),]
  colnames(GT) = paste0(vcf$chr,':',vcf$pos)
  
  ### Get info about mutation types
  MUT = gsub('.*EFF=(.*)', '\\1',vcf$info)
  MUT = strsplit(MUT, ',')
  MUT = lapply(MUT, function(x) {return(x[regexpr(gene, x)!=-1])}) ## keep only mutations with an effect gene of interest
  MUT = lapply(MUT, function(x) {unique(gsub('(.*)\\(.*','\\1',x))}) ## simplify mutation description
  names(MUT) = colnames(GT)
  
  
  ## Filter for missingness and MAF
  GT = GT[,apply(GT,2, function(x) {mean(is.na(x))})<MN.T,drop=F]
  MAF = apply(GT, 2, mean, na.rm=T)/2
  GT = GT[,MAF>MAF.T,drop=F]
  MUT = MUT[colnames(GT)]
  
  ## return GT and MUT table
  return(list(GT,MUT))
  
}


# function to get correlation of SNPs with an environmental variable
getGEAcor = function(gt, env) {
  
COR = cor(gt, env, use='pairwise.complete.obs')
return(COR)
}


# function to get RF importance for genotype vs. an environmental variable
GENOPHENO_IMP=function(GT, PH, PS, vect=F, N=50, rep=100) {
  
  ## Remove individuals without genotype or phenotypes
  PH = na.omit(PH)
  GT = na.omit(GT)
  
  ## container of r2
  r2=c()
  IMP = c()
  
  
  for (i in 1:rep) {
    
    ## Keep only individuals with phenotype and genotype
    shared_accessions = names(which(table(c(rownames(PH), rownames(GT)))==2))
    
    if (length(shared_accessions)<N) {return()} # if there are no sufficient N with pheno and geno data: abort...
    
    ### Restrict N: stratified sampling across the phenotype range
    Y = PH[shared_accessions,]
    cls = cutree(hclust(dist(Y)), k=N)
    names(Y) = names(cls) = shared_accessions
    shared_accessions = c()
    for (i in unique(cls)) {shared_accessions = c(shared_accessions,sample(names(Y[cls==i,drop=F]), 1))}
    
    ### Preapre data matrix: DATA1: pop structure only, DATA2: pop structure + GT
    DATA = data.frame(PS[shared_accessions,], GT[shared_accessions,], 'Y'=scale(PH[shared_accessions,]))
    
    RF.FULL=randomForest(Y ~ . , data=DATA, ntree=N*5)
    
    R2i = c()
    
    for (var in colnames(DATA)[-ncol(DATA)]) {
      
      RF.I = randomForest(Y ~ . , data=DATA[,colnames(DATA)!=var], ntree=N*5)
      
      RF.I.r2 = mean(RF.I$rsq)
      RF.I.r2[RF.I.r2<0] = 0
      R2i[var] = RF.I.r2
    }
    
    ## 
    R2.full = mean(RF.FULL$rsq)
    
    R2i = R2i - R2.full
    
    
    ## Get importance of every explanatory variable in model 2
    IMP = cbind(IMP, R2i)
    
  }
  
  
  return(-IMP)
  
  
}


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

gene = 'AT1G74960' # GRODD
ph_id = '22_GR63 warm'

sEXONS = EXONS[substr(EXONS$geneID,1,nchar(EXONS$geneID)-2)==gene,]
load(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_',gene,'.rda'))
GT=SNP_AT[[1]]
MUT=SNP_AT[[2]]

PH = na.omit(PHENO[,ph_id,drop=F])

## plot SNPs
IMP = GENOPHENO_IMP(GT, PH, PS)

IMP.SNPS = apply(IMP,1,mean)[10:nrow(IMP)]
POS = as.numeric(substr(names(MUT), 3, nchar(names(MUT))))

s_EXONS = EXONS[EXONS$geneID==gene,]



{
par(mfrow=c(2,1));par(mar=c(3,3,0,0))

plot(NA, xlim=range(min(s_EXONS$V4), max(s_EXONS$V5)), ylim=c(0,max(IMP.SNPS)), axes=F)

for (i in 1:length(IMP.SNPS)) {
  
  pos = POS[i]
  imp = IMP.SNPS[i]

  ## add importance of every SNP
  rect(pos-10, 0, pos+10, imp, border=NA, col='grey30')
    
}

axis(1)

for (mut in rownames(MUT_PCH)) {
  
  pos.mut = POS[unlist(lapply(MUT, function(x) {mut%in%x}))]
  
  if (length(pos.mut)>0) {
  
  points(pos.mut, rep(0, length(pos.mut)), pch=as.numeric(MUT_PCH[mut,'pch']), bg='yellow', cex=1.5)
  text(pos.mut, rep(0, length(pos.mut)), MUT_PCH[mut,'txt'], cex=0.5)
  }
}


par(mar=c(0,3,0,0))
plot(NA, xlim=range(min(s_EXONS$V4), max(s_EXONS$V5)), ylim=c(-length(unique(s_EXONS$transcriptV)),-0.5), axes=F)

for (i in unique(s_EXONS$transcriptV)) {
  
  exons = s_EXONS[s_EXONS$transcriptV==i,]
  
  for (e in 1:nrow(exons)) {
    rect(exons$V4[e], -as.numeric(i), exons$V5[e], -as.numeric(i)+0.5)
    rect(exons$V4[e], -as.numeric(i), exons$V5[e], -as.numeric(i)+0.5)
  }
  
  
}
}



# try to make a plot where, 
# --> in the middle, genome line, with position of SNPs. Symbols for different type of mutations. 
# --> on top: barplot with effect of every SNP
# --> below: exons as boxes,on multiple lines if multiple versions
# --> for one SNP, add an asterisk or something, then show a boxplot showing GwP on the right.
# --> add 2 p-values for GxP

