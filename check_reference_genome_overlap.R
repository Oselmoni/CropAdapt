library(vcfppR)
library(vcfR)
library(Rsamtools)

## get list of  datsets
DS = list.dirs('DATA/PROCESSED_VCF/', full.names = F)[-1]

## for every dataset
for (ds in DS) {
  
  ### load a vcf
#  SNPS = vcftable('DATA/PROCESSED_VCF/barley/processed.vcf.gz', region=region )
 
  VCF = read.vcfR('DATA/PROCESSED_VCF/tomato_wr/processed.vcf.gz', nrows = 100)
  
  
  ### load fasta reference
  refPath = "DATA/RAW/tomato_wr/GCF_000188115.4_SL3.0_genomic.fna"

  fa <- FaFile(refPath)
  open(fa)
  
  
  ### scan SNPs
  meta_snps = VCF@fix
  
  for (i in sample(nrow(meta_snps), 5)) {

    chrom = VCF@fix[i,1]
    pos = as.numeric(VCF@fix[i,2])
    nu =  VCF@fix[i,4]
    
    seq <- scanFa(fa, param = GRanges('NC_015438.3', IRanges(pos, pos)))
    

    print(paste0(nu,'   vs.   ',seq))
    
  }
  
  



reference = read.fasta('DATA/RAW/barley/barley_morex_pseudomolecules.fasta')

fa <- FaFile("DATA/RAW/barley/barley_morex_pseudomolecules.fasta")
open(fa)



SNPS$ref[20]

unique(SNPS$pos)[1]






