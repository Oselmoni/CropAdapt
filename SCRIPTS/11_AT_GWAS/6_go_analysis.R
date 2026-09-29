library(xml2)
### Get list of topgenes for gene set enrichment analysis with Panther

### Load list of top candidate genes from LFMM picmin
load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
TOPGENES.LFMM = TOPGENES
TOPGENES.LFMM$GEA = 'LFMM'
### Load list of top candidate genes from KTAU picmin
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
TOPGENES.KTAU = TOPGENES
TOPGENES.KTAU$GEA = 'KTAU'

### Create unique table of top genes
TOPGENES = rbind(TOPGENES.LFMM, TOPGENES.KTAU)


### Get list of A. thaliana orthologs
AT = unique(gsub('(.*)\\..*', '\\1',unlist(strsplit(TOPGENES$OG, ', '))))

### Save list of top genes
write.table(AT, 'DATA/ATHALIANA/topgenes_gsea.txt', quote=F, row.names=F, col.names=F)

###
###
### --> Run the analysis on https://pantherdb.org, using the A. thaliana genome as refeence
###
###

### read output from Panther
panther_xml = as_list(read_xml('DATA/ATHALIANA/panther_analysis.xml', encoding = "UTF-8"))

### load list of GO terms to no annotate (too generic)
GO_exclude = read.table('DATA/ATHALIANA/gocheck_do_not_annotate.tsv', header=T)
GO_exclude$ID = paste0('GO:',gsub('.*GO_(.*)>', '\\1',GO_exclude$X.x))

# create table
GO.tab = data.frame()

## browse panther groups 
for (i in which(names(panther_xml$overrepresentation)=='group')) {
  
  group = panther_xml$overrepresentation[[i]]
  
  ## browse every GOterm
  for (go in 1:length(group)) {
    
    id = group[[go]]$term$id[[1]]
    
    if (id%in%GO_exclude$ID==F) {
    
    # create a first column merging label and level
    level = group[[go]]$term$level
    label = group[[go]]$term$label
    
    if (level==0) {name=label 
    } else {
    name = paste0(paste(rep(' ', as.numeric(level)), collapse = ''),'->',label)
    }
    
    GO.tab[id,'name'] = name
    GO.tab[id,'id'] = id
    GO.tab[id,'level'] = level
    
    GO.tab[id,'group'] = i-min(which(names(panther_xml$overrepresentation)=='group'))+1
    
    GO.tab[id,'N_reference'] = group[[go]]$number_in_reference
    GO.tab[id,'N_significant'] = group[[go]]$input_list$number_in_list
    GO.tab[id,'Expected'] = group[[go]]$input_list$expected
    GO.tab[id,'pval'] = group[[go]]$input_list$pValue
    GO.tab[id,'fold'] = group[[go]]$input_list$fold_enrichment
    GO.tab[id,'FDR'] = group[[go]]$input_list$fdr
    
    GO.tab[id,unlist(group[[go]]$input_list$mapped_id_list)] = 'x'
    }
  }
  
}


## save as table
write.table(GO.tab, file = 'FIGURES/athaliana/GOterms_panther.txt',row.names=F, col.names=T, quote=F, sep='\t') 

## save as rddata
save(GO.tab, file='DATA/ATHALIANA/GOtab.rda')


