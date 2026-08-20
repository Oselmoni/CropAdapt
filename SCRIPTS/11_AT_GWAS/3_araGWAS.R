
###  Load set genes of interest
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
genes_of_interest_KTAU = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
genes_of_interest_LFMM = unique(gsub('(.*)\\..*','\\1',unlist(strsplit(TOPGENES$OG, ', '))))

genes_of_interest = unique(c(genes_of_interest_KTAU, genes_of_interest_LFMM))



# Get araGWAS hits
araGWAS_reports = list.files('DATA/ATHALIANA/RAW/araGWAS/', full.names = T) # get list of csv files with araGWAS results for specific genes

araGWAS = data.frame()
for (gene in genes_of_interest) { 

  ## import aragwas report for gene of interest
  agwas = read.csv(  araGWAS_reports[regexpr(gene, araGWAS_reports)!=-1]) 
  
  ## filter to keep only snp on gene
  agwas = agwas[agwas$snp.gene_name==gene,]

  ## add to container
  araGWAS = rbind(araGWAS, agwas[,c('snp.gene_name','snp.chr','snp.position','score','mac','maf','study.id','study.name','study.phenotype.name','overPermutation','overBonferroni','snp.annotations.0.effect')])
  
}
  
# get rid of "environmental phenotypes"
araGWAS = araGWAS[substr(araGWAS$study.phenotype.name,1,4)!='clim',]


# save table of GWAS hits
save(araGWAS, file='DATA/ATHALIANA/araGWAS.rda')


