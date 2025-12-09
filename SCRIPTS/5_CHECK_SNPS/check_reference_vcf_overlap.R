library(vcfppR)
library(Rsamtools)

## load list of datasets and path to reference genomes
ref_lists = read.csv('DATA/REF_TO_VCF/paths_to_references.csv', header=F, row.names=1)


## prepare output
OUTCHECK = data.frame()

## for every dataset...
for (d in rownames(ref_lists)) {

print(d)
## Load VCF
vcf = vcftable(paste0('DATA/PROCESSED_VCF/',d,'///processed.vcf.gz_pruned.vcf.gz'))

## load fasta of reference genome
fa <- FaFile(ref_lists[d,1])
open(fa) # create index


## load table linking chromsome names from reference to vcf
ref_to_vcf = read.csv(paste0('DATA/REF_TO_VCF/',d,'.csv'), colClasses = c('character','character'))

## check chromosomes names
CHRS_FA = ref_to_vcf$reference
CHRS_VCF = ref_to_vcf$vcf



### for every chromosome: pick 10 SNPs randomly, and check if they exist on reference
for (ci in 1:length(CHRS_VCF)) {
 
  chr = CHRS_FA[ci]
  
  ## randomly pick 10 SNPs on chromosome
  idx = sample(which(vcf$chr==CHRS_VCF[ci]), 20)
  
  for (i in idx) {

    pos = vcf$pos[i]
    ref = vcf$ref[i]
    alt = vcf$alt[i]
    
    nucl = c(substr(vcf$id[i], nchar(vcf$id[i])-2, nchar(vcf$id[i])-2),
    substr(vcf$id[i], nchar(vcf$id[i]), nchar(vcf$id[i])))
  
    
    # find nucleotide in fasta
    ref_fa <- scanFa(fa, param = GRanges(chr, IRanges(pos, pos)))
    
    # check match
    match = as.character(ref_fa)==ref|as.character(ref_fa)%in%unlist(strsplit(alt,','))

    # add to output
    OUTCHECK = rbind(OUTCHECK, data.frame(d, ci, pos, match))
    
  } 
  
}


OUTCHECK = OUTCHECK[order(OUTCHECK$ci, OUTCHECK$pos),]
print(table(OUTCHECK$match, OUTCHECK$d))


gc()

}




### check that all SNPs matched with reference
table(OUTCHECK$match, OUTCHECK$d)




