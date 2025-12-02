library(vcfppR)
library(Rsamtools)

## load list of datasets and path to reference genomes
ref_lists = read.csv('DATA/REF_TO_VCF/paths_to_references.csv', header=F, row.names=1)

## prepare output
OUTCHECK = data.frame()

## for every dataset...
#for (d in rownames(ref_lists)) {

d=rownames(ref_lists)[6]
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



##########
chrs_fa = seqlevels(fa)
CHRS_FA = chrs_fa[1:length(CHRS_VCF)]
chrs_vcf = sort(unique(vcf$chr))
x=as.matrix(cbind(chrs_vcf, chrs_fa[1:length(chrs_vcf)]))
########


### for every chromosome: pick 10 SNPs randomly, and check if they exist on reference
for (ci in 1:length(CHRS_VCF)) {
  
  chr = CHRS_FA[ci]
  
  ## randomly pick 10 SNPs on chromosome
  idx = sample(which(vcf$chr==CHRS_VCF[ci]), 10)
  
  for (i in idx) {

    pos = vcf$pos[i]
    ref = vcf$ref[i]
    alt = vcf$alt[i]
    
    # find nucleotide in fasta
    ref_fa <- scanFa(fa, param = GRanges(chr, IRanges(pos, pos)))
    
    # check match
    match = as.character(ref_fa)==ref|as.character(ref_fa)==alt
    
    # add to output
    OUTCHECK = rbind(OUTCHECK, data.frame(d, ci, pos, match))
    
  } 
  
}

print(table(OUTCHECK$match, OUTCHECK$ci))


OUTCHECK = OUTCHECK[order(OUTCHECK$ci, OUTCHECK$pos),]

plot(OUTCHECK$match, col=as.factor(OUTCHECK$ci), pch=16, cex=0.5)

boxplot(OUTCHECK$pos~OUTCHECK$match)

gc()

#}

### check that all SNPs matched with reference
table(OUTCHECK$match, OUTCHECK$d)




