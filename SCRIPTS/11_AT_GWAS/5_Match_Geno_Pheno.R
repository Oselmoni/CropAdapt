library(vcfppR)
library(openxlsx)
library(terra)
source('SCRIPTS/custom_R_functions.R')
library(randomForest)
library(foreach)
library(doParallel)


### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID


### Load env data
load('DATA/ATHALIANA/meta.rda')

### Load Genes info on AT
load('DATA/ATHALIANA/GENE.rda')

#setup parallel backend to use many processors
cores = 6
cl <- makeCluster(cores-1) 
registerDoParallel(cl)


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

## Keep only phenotype with at least 50 N
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
  vcf = vcftable(paste0('DATA/ATHALIANA/RAW/AT_chr',chr,'_eff.vcf.gz'), region = paste0(chr,':',sta,'-',end), vartype = 'snps')
  
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


# function to build geno-pheno association
GENOPHENO=function(gene, PHENO, PS, vect=F, N=50) {
  
  ##### Get genotype data for gene
  load(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_',gene,'.rda'))
  GT = SNP_AT[[1]] # genotype table

  
  ## If no SNPS left on gene, abort
  if (is.null(GT)) {return()}
  if (ncol(GT)<=1|nrow(GT)<=1) {return()}
  
  ## Prepare SNP matrix: discard indiiduals with missing GT
  GT = na.omit(GT)
  
  ## Correct SNP identifiers to be compatible with RF
  colnames(GT) = paste0('snps', 1:ncol(GT))
  
  ## screen all gwas studies
  phenoscore = do.call(rbind, foreach(ph=colnames(PHENO)) %dopar% {

    library(randomForest)
    library(missForest)
   
    ## Keep only individuals with phenotype of interest
    PHENO_i = na.omit(PHENO[,ph,drop=F])
    
    ## Keep only individuals with phenotype and genotype
    shared_accessions = names(which(table(c(rownames(PHENO_i), rownames(GT)))==2))

    if (length(shared_accessions)<N) {return()} # if there are no sufficient N with pheno and geno data: abort...
    
    ### Restrict N: stratified sampling across the phenotype range
    Y = PHENO_i[shared_accessions,]
    cls = cutree(hclust(dist(Y)), k=N)
    names(Y) = names(cls) = shared_accessions
    shared_accessions = c()
    for (i in unique(cls)) {shared_accessions = c(shared_accessions,sample(names(Y[cls==i,drop=F]), 1))}

    ### Preapre data matrix: DATA1: pop structure only, DATA2: pop structure + GT
    DATA1 = data.frame(PS[shared_accessions,], 'Y'=scale(PHENO_i[shared_accessions,]))
    DATA2 = data.frame(PS[shared_accessions,], GT[shared_accessions,], 'Y'=scale(PHENO_i[shared_accessions,]))

    
    ## Run Random Forest
    RF1=randomForest(Y ~ . , data=DATA1, ntree=N*5)
    RF2=randomForest(Y ~ . , data=DATA2, ntree=N*5)
    
    ## Extract pseudo-R2, id < 0 -> 0
    RF1r2  = mean(RF1$rsq)
    RF2r2 =  mean(RF2$rsq)
    RF1r2[RF1r2<0] = 0
    RF2r2[RF2r2<0] = 0
    
    ## Calcualte ∆R2 and MSE
    r2 = RF2r2 - RF1r2
    mse = mean(RF2$mse)-mean(RF1$mse)

    ## Prepare output
    out = data.frame()
    out[1,'gene'] = gene
    out[1,'ph'] = ph
    out[1,'dR2'] = r2
    out[1,'dMSE'] = mse
    out[1,'R2_PS'] = RF1r2
    out[1,'R2_PS_GT'] = RF2r2
    out[1,'nind'] = nrow(DATA2)
    out[1,'nsnps'] = ncol(GT)
    
    return(out)
  })
  
  if (is.null(phenoscore)) {return()} # if no phenoscore .. abort
  
  ## sort phenotypes based on dR2
  phenoscore = phenoscore[order(phenoscore$dR2, decreasing = T),]

  ## if vector format: return all phenotypes...
  if (vect==T) {return(phenoscore)}
  
  ## otherwise, return just best phenotype 
  phenoscore = phenoscore[which.max(phenoscore$dR2),]
  return(phenoscore)
  
}


# function to test geno-pheno association and test randomness. 
GENOPHENOTEST = function(gene) {
  print(gene)
  ### Step 1: Find best phenotype association for  gene of interest
  
  print('Looking for best phenotype association for gene...')
  ## Prepare container of different runs for a given gene
  TOP.OUT = data.frame() 
  set.seed(123);for (i in 1:100) { # Run 100 replicates for a given gene

    top=GENOPHENO(gene=gene, PHENO=PHENO, PS=PS, vect=T) # check phenotype associations 
    
    if (is.null(top)) {break}
    
    TOP.OUT[top$ph,i] = top$dR2 # fill in output table with dR2 of all phenotpes
    
  }
  
  if (nrow(TOP.OUT)==0) {print('no SNP found');return()}
  
  #### Step 2: check: How likely it is to obtain this results with a random phenotype? 
  
  print('Checking gene association with random phenotypes...')
  ### Create container of permutations runs
  topgene_shufflePheno = data.frame()
  
  ### Continue until 1000 permutations are reached
  while (nrow(topgene_shufflePheno)<100) {
    
    ## calculate permutated phenotypes left to compute
    pheno_missing = 100-nrow(topgene_shufflePheno)
    
    ### Create random phenotypes
    sh_phenos = sample(colnames(PHENO), pheno_missing, T)
    SHF_PHENO = c()
    for (shp in sh_phenos) {
      
      x= PHENO[,shp]  
      
      notNAs = is.na(x)==F
      
      x[notNAs] = sample(x[notNAs]) # shuffle phenotypes randomly
      
      SHF_PHENO = cbind(SHF_PHENO, x)
    }
    rownames(SHF_PHENO) = rownames(PHENO)
    colnames(SHF_PHENO) = paste0(1:pheno_missing,'_',sh_phenos)
    
    # Re run: same locus, 1000 random phenotypes
    topgene_shufflePheno = rbind(topgene_shufflePheno, GENOPHENO(gene=gene, PHENO=SHF_PHENO, PS=PS, vect=T))
  }
  

  ### Step3, check: how likely it is to obtain same association for same phenotype random genes? 
  print('Checking gene association with random genotypes...')
  
  # Create container
  rndgene_samePheno=data.frame()
  
  while (nrow(rndgene_samePheno)<=100) {
    
    rndgene = sample(GENE$geneID,1)
    SNP_AT = getGENO(rndgene)
    save(SNP_AT, file=paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_temp.rda')) # create a temporary file with SNPs of a random gene
    rndgene_samePheno = rbind(rndgene_samePheno, GENOPHENO(gene='temp', PHENO=PHENO[,names(which.max(apply(TOP.OUT,1,mean))),drop=F], PS=PS))
    file.remove(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_temp.rda')) # remove temporary file
  }

  ### Store output
  OUT = list('GxP'=TOP.OUT, 'RND.PH'=topgene_shufflePheno, 'RND.GT'=rndgene_samePheno)
  
  save(OUT, file=paste0('DATA/ATHALIANA/GwP/GwP_',gene,'.rda'))
  
}


GENOPHENOTEST2 = function(gene) {
  print(gene)
  ### Step 1: Find best phenotype association for  gene of interest
  
  if (file.exists(paste0('DATA/ATHALIANA/GwP/GwP_',gene,'.rda'))) {load(paste0('DATA/ATHALIANA/GwP/GwP_',gene,'.rda'))} else {break}
  TOP.OUT = OUT$GxP
  
  if (nrow(TOP.OUT)==0) {print('no SNP found');return()}
  
  #### Step 2: check: How likely it is to obtain this results with a random phenotype? 
  
  print('Checking gene association with random phenotypes...')
  ### Create container of permutations runs
  topgene_shufflePheno = data.frame()
  
  ### Continue until 1000 permutations are reached
  while (nrow(topgene_shufflePheno)<100) {
    
    ## calculate permutated phenotypes left to compute
    pheno_missing = 100-nrow(topgene_shufflePheno)
    
    ### Create random phenotypes
    sh_phenos = sample(colnames(PHENO), pheno_missing, T)
    SHF_PHENO = c()
    for (shp in sh_phenos) {
      
      x= PHENO[,shp]  
      
      notNAs = is.na(x)==F
      
      x[notNAs] = sample(x[notNAs]) # shuffle phenotypes randomly
      
      SHF_PHENO = cbind(SHF_PHENO, x)
    }
    rownames(SHF_PHENO) = rownames(PHENO)
    colnames(SHF_PHENO) = paste0(1:pheno_missing,'_',sh_phenos)
    
    # Re run: same locus, 1000 random phenotypes
    topgene_shufflePheno = rbind(topgene_shufflePheno, GENOPHENO(gene=gene, PHENO=SHF_PHENO, PS=PS, vect=T))
  }
  
  
  ### Step3, check: how likely it is to obtain same association for same phenotype random genes? 
  print('Checking gene association with random genotypes...')
  
  # Create container
  rndgene_samePheno=data.frame()
  
  while (nrow(rndgene_samePheno)<=100) {
    
    rndgene = sample(GENE$geneID,1)
    SNP_AT = getGENO(rndgene)
    save(SNP_AT, file=paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_temp.rda')) # create a temporary file with SNPs of a random gene
    rndgene_samePheno = rbind(rndgene_samePheno, GENOPHENO(gene='temp', PHENO=PHENO[,names(which.max(apply(TOP.OUT,1,mean))),drop=F], PS=PS))
    file.remove(paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_temp.rda')) # remove temporary file
  }
  
  ### Store output
  OUT = list('GxP'=TOP.OUT, 'RND.PH'=topgene_shufflePheno, 'RND.GT'=rndgene_samePheno)
  
  save(OUT, file=paste0('DATA/ATHALIANA/GwP/GwP_',gene,'.rda'))
  
}


# set genes of interest
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
genes_of_interest_KTAU = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
genes_of_interest_LFMM = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

genes_of_interest = unique(c(genes_of_interest_KTAU, genes_of_interest_LFMM))

# Run pheno test for every gene
for (gene in genes_of_interest) { GENOPHENOTEST(gene)}



####
####
####. >>>> below HERE. - - to be checked or deleted
####
####
#### Check Geno-Pheno association for locus of interest
gene = 'AT1G74960' # GRODD
gene = 'AT2G40830' # VSOILW

## Prepare container of different runs for a given gene
TOP.OUT = data.frame() 
set.seed(123);for (i in 1:100) { # Run 100 replicates for a given gene
  print(i)
  top=GENOPHENO(gene=gene, PHENO=PHENO, PS=PS, vect=T) # check phenotype associations 
  
  TOP.OUT[top$ph,i] = top$dR2 # fill in output table with dR2 of all phenotpes
  
}


# show boxplot of importantce of phenotypes
par(mar=c(15,3,1,1));boxplot(t(TOP.OUT), las=2) 

sort(apply(TOP.OUT,1,mean), decreasing=T)[1:5] # get top phenotypes
# for AT1G74960 ->  22_GR63 warm, dR2 ~ 0.11
# for AT2G40830 -> 36_LRL125, dR ~0.055

TOP.OUT.S = data.frame(t(TOP.OUT)[,names(sort(apply(TOP.OUT,1,mean), decreasing=T))[1:5]])
TOP.OUT.S$others =  apply(t(TOP.OUT)[,names(sort(apply(TOP.OUT,1,mean), decreasing=T))[6:nrow(TOP.OUT)]],1,mean)
boxplot(TOP.OUT.S, las=2)


###############################################################
###############################################################

# Visualize score of best real phenotypes vs 1000 random phenotypes
hist(topgene_shufflePheno$dR2, breaks=10, main='')
abline(v=max(apply(TOP.OUT,1,mean), na.rm=T))
mean(max(apply(TOP.OUT,1,mean), na.rm=T)<topgene_shufflePheno$dR2) # empirical p


###############################################################
###############################################################
#### 2nd check: How likely it is to obtain this results with a random gene? 



hist(rndgene_samePheno$dR2, breaks=100)
abline(v=max(apply(TOP.OUT,1,mean), na.rm=T))
mean(max(apply(TOP.OUT,1,mean), na.rm=T)<rndgene_samePheno$dR2) # empirical p




## Stop Cluster
stopCluster(cl)





