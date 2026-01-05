library(vcfppR)
library(seqinr)
library(foreach)
library(doParallel)


### Load gtf file
GTF = read.table('DATA/RAW/sorghum_la//Sbi1.4.gff3', sep='\t')

### load gene names
fasta = read.fasta('DATA/RAW/PROTEINS/sorghum_Sbi1.4.pep.fa')
nameGenes=names(fasta)
rm(fasta);gc()

# Load VCF
vcf = vcftable('DATA/PROCESSED_VCF/sorghum_la/processed.vcf.gz_pruned.vcf.gz')

snpid = paste0(vcf$chr,':',vcf$pos)
chr = vcf$chr
pos = as.numeric(vcf$pos)
rm(vcf);gc()

names(pos)=names(chr)=snpid

#### Get protein to exon annotation

# keep only annotation of interest
GTF=GTF[GTF$V3=='CDS',]


GTF$GENEID = gsub('.*Parent=(.*)$', '\\1', GTF$V9)



## double check that ~all gene identifiers match fasta headers
table(GTF$GENEID%in%nameGenes)
table(nameGenes%in%GTF$GENEID)


## match chromsome IDs
chrID = as.character(1:10)
names(chrID) = paste0('chromosome_',1:10)

GTF$V1 = chrID[GTF$V1]



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


## random check 
randomSNPS = SNP_ANNOTATION[sample(which(is.na(SNP_ANNOTATION$geneID)==F), 1),]
randomSNPS
GTF[GTF$GENEID%in%randomSNPS$geneID,]

#save(SNP_ANNOTATION, file='DATA/SNP_ANNOTATION/sorghum_la.rda', compression_level = 9)

