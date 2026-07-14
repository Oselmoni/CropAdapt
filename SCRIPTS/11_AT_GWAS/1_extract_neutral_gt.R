library(vcfppR)

### Get info on chromosome length
GFF3 = read.table('DATA/ATHALIANA/RAW/Arabidopsis_thaliana.TAIR10.62.gff3', sep='\t')
CHROM = GFF3[GFF3$V3=='chromosome',]


### Sample 10,000 SNPs
n=10000
CHROM$REL = CHROM$V5/sum(CHROM$V5) # get relative length of every chromsoome
CHROM$N = round(CHROM$REL*n) # get number of SNPs per chromosome
CHROM$cumN = cumsum(CHROM$N) # get cumulative number of SNPs 

### For every chromosome... 
SNPS = data.frame(nrow=1135, col=0)

for (chr in 1:5) {

  print(chr)
  ## set target number of SNPs
  target = CHROM[CHROM$V1==chr,'cumN']
  
  ## as long as the target number of SNPs per chromosome is not reached... 
  while (ncol(SNPS)<target) {
    
    # sample a position on the chromosome
    pos = sample(1000:(CHROM[CHROM$V1==chr,'V5']-1000), 1)
    sta = pos-1000
    end = pos+1000
    
    # get VCF around the region
    vcf = vcftable(paste0('DATA/ATHALIANA/RAW/AT_chr',chr,'_eff.vcf.gz'), region = paste0(chr,':',sta,'-',end), vartype = 'snps')
    GT = t(vcf$gt)
    
    if (ncol(GT)>=1&nrow(GT)>1) {

    # set individuals and snp ids
    colnames(GT) = paste0(vcf$chr,':',vcf$pos)
    rownames(GT) = vcf$samples
    
    
    ## Filter for missingness and MAF
    GT = GT[,apply(GT,2, function(x) {mean(is.na(x))})<0.2,drop=F]
    MAF = apply(GT, 2, mean, na.rm=T)/2
    GT = GT[,MAF>0.01,drop=F]

    if (ncol(GT)>=1) {
    
    # sample one SNP in the region
    SNPS = cbind(SNPS, GT[,sample(colnames(GT), 1),drop=F])
    
    }}
  }
}

### Save
save(SNPS, file='DATA/ATHALIANA/SNPS.rda', compress = T)
