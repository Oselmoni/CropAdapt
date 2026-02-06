library(vcfppR)
library(seqinr)
library(foreach)
library(doParallel)


### Load gtf file
GTF = read.table('DATA/RAW/barley/Hv_IBSC_PGSB_r1_HighConf.gtf')
# keep only annotation of interest
GTF=GTF[GTF$V3=='CDS',]


### load gene names
fasta = read.fasta('DATA/RAW/PROTEINS/barley_Hv_IBSC_PGSB_r1_proteins_HighConf.fa')
nameGenes=names(fasta)
rm(fasta);gc()

# Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/barley///processed.vcf.gz_pruned.vcf.gz')

snpid = paste0(vcf$chr,':',vcf$pos)
chr = vcf$chr
pos = as.numeric(vcf$pos)
rm(vcf);gc()

names(pos)=names(chr)=snpid


## add column with genes ID
GTF$GENEID = GTF$V13


## double check that ~all gene identifiers match fasta headers
table(GTF$GENEID%in%nameGenes)
table(nameGenes%in%GTF$GENEID)


##
SNP_ANNOTATION = data.frame()

for (CH in unique(chr)) {
  
  print(CH)
  ## find SNPs on chromosome
  pos_CH = pos[chr==CH]
  
  ## find GENES on chromosome
  GTF_CH = GTF[GTF$V1==CH,]
  
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-2) #not to overload your computer
  registerDoParallel(cl)
  
  snp_annotation_chr <- foreach(snp=names(pos_CH), .combine=rbind) %dopar% {
    
    # get SNP position
    ppp = as.numeric(pos_CH[snp])
    
    # get geneID of overlapping genes
    matches = unique(GTF_CH[which(ppp>=GTF_CH$V4&ppp<=GTF_CH$V5),'GENEID'])
    
    if (length(matches)>0) {
      return(data.frame(snp, 'chr'=CH, 'pos'=ppp, 'geneID'=matches))
    } else {
       return(data.frame(snp, 'chr'=CH, 'pos'=ppp, 'geneID'=NA))
      
    }
  }
  
  stopCluster(cl)  
  
  SNP_ANNOTATION = rbind(SNP_ANNOTATION, snp_annotation_chr)
  
}

save(SNP_ANNOTATION, file='DATA/SNP_ANNOTATION/barley.rda', compression_level = 9)

### Save genes metadata
load('DATA/SNP_ANNOTATION/barley.rda')
CDS = GTF[,c('V1','V4','V5','V7','GENEID')] 
colnames(CDS) = c('CHR','STA','END','STR','GENEID')
CDS = CDS[CDS$GENEID%in%SNP_ANNOTATION$geneID,]
save(CDS, file='DATA/GENE_ANNOTATION/barley.rda')
  
