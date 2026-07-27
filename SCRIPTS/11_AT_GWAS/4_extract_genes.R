### load gene lists
GFF = read.table('DATA/ATHALIANA/RAW/Arabidopsis_thaliana.TAIR10.62.gff3', sep='\t', quote='')
GENE = GFF[GFF$V3=='gene',]
GENE$geneID = gsub('ID=gene:(.*);Name.*','\\1',GENE$V9)

save(GENE, file='DATA/ATHALIANA/GENE.rda')

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


# set genes of interest
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
genes_of_interest_KTAU = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
genes_of_interest_LFMM = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

genes_of_interest = unique(c(genes_of_interest_KTAU, genes_of_interest_LFMM))

# get SNP information for genes of interest and save
for (gene in genes_of_interest) {
  
  SNP_AT = getGENO(gene)
  
  save(SNP_AT, file=paste0('DATA/ATHALIANA/SNP_AT/SNP_AT_',gene,'.rda'))
  
  
}

