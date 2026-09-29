### Load slist of datasets
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')
rownames(dslist) = dslist$dataset

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

## correct names of columns with dataset ids
col.ds = colnames(TOPGENES)[substr(colnames(TOPGENES),1,2)=='P(']
colnames(TOPGENES)[colnames(TOPGENES)%in%col.ds]=paste0('P(',dslist[substr(col.ds, 3, nchar(col.ds)-1),'name'],')')

## add pheno metadata
load('DATA/ATHALIANA/PHENO_META.rda')
rownames(PHENO.META) = PHENO.META$phenotype_id


## Write output
write.table(TOPGENES,'FIGURES/GEAs/TOPGENES.tsv', sep='\t', quote=F, row.names = F, col.names=T)






### Load UNIPROT table of AT annotations
UNIPROT = read.table('DATA/ATHALIANA/RAW/uniprotkb_taxonomy_id_3702_2026_07_15.tsv', sep='\t', quote="", header=T, fill = T)
UNIPROT = UNIPROT[order(UNIPROT$Reviewed),] # set reviewed prot on top



### Load araGWAS results
load('DATA/ATHALIANA/araGWAS.rda')



### Retrieve protein annotations
PROT.ANN = data.frame()
       
for (og in unique(TOPGENES$locus)) {

    AT = unique(TOPGENES[TOPGENES$locus==og,'OG'])
    AT = unlist(strsplit(AT, ', '))
    AT = unique(gsub('(.*)\\..*','\\1',AT))
    
    if (length(AT)>0) { 
  
    
    for (at in AT) {

      ### Find uniprot annotations
      UP_name =   UNIPROT$Protein.names[UNIPROT$TAIR==paste0(at,';')][1]
      UP_ID =   UNIPROT$Entry.Name[UNIPROT$TAIR==paste0(at,';')][1]
      MF =   UNIPROT$Gene.Ontology..molecular.function.[UNIPROT$TAIR==paste0(at,';')][1]
      BP =   UNIPROT$Gene.Ontology..biological.process.[UNIPROT$TAIR==paste0(at,';')][1]
      CC =   UNIPROT$Gene.Ontology..cellular.component.[UNIPROT$TAIR==paste0(at,';')][1]
      
      if (length(UP_name)==0) {UP_name=''}
      if (length(UP_ID)==0) {UP_ID=''}
      if (length(MF)==0) {MF=''}
      if (length(BP)==0) {BP=''}
      if (length(CC)==0) {CC=''}
      
      
      ### Find GTxPhenotype associations
      gwas = araGWAS[araGWAS$snp.gene_name==at,]
      
      if (nrow(gwas)>0) {
      
      gwas$pos = paste0(gwas$snp.chr,':',gwas$snp.position)
      
      ## get info of top phenotype
      top_pheno_pos =   gwas$pos[which.max(gwas$score)]
      top_pheno_type =   gwas$snp.annotations.0.effect[which.max(gwas$score)]
    
      top_pheno = gwas$study.phenotype.name[which.max(gwas$score)]
      top_pheno_id = gwas$study.id[which.max(gwas$score)]
      top_pheno_score = max(gwas$score)
      top_pheno_bonferroni = gwas$overBonferroni[which.max(gwas$score)]
      top_pheno_permutation = gwas$overPermutation[which.max(gwas$score)]
      
      ## get summary of type of all phenotypes associated
      all_pheno = paste(unique(PHENO.META[as.character(gwas$study.id),'to_name']), collapse = ';')
      
      } else {
        top_pheno_pos = top_pheno_type = top_pheno = top_pheno_score = top_pheno_bonferroni = top_pheno_permutation = top_pheno_id = all_pheno = ''
      }
    
      
      PROT.ANN = rbind(PROT.ANN, data.frame('OrthogeneID'=og, 'A. thaliana OG'=at, 'UP_ID'=UP_ID, 'UP_Name'=UP_name, 'GO.MF'=MF, 'GO.BP'=BP, 'GO.CC'=CC, top_pheno_pos, top_pheno_type, top_pheno_id, top_pheno, top_pheno_score, top_pheno_bonferroni, top_pheno_permutation, all_pheno))
      
    }
      
    }

}

head(PROT.ANN)

### add info on phenotypes
PROT.ANN$Pheno = PHENO.META[PROT.ANN$top_pheno_id,'name']
PROT.ANN$Study = PHENO.META[PROT.ANN$top_pheno_id,'study']
PROT.ANN$Scoring = PHENO.META[PROT.ANN$top_pheno_id,'scoring']
PROT.ANN$TO_name = PHENO.META[PROT.ANN$top_pheno_id,'to_name']
PROT.ANN$TO_definition = PHENO.META[PROT.ANN$top_pheno_id,'to_definition']

PROT.ANN$Growth.Conditions = PHENO.META[PROT.ANN$top_pheno_id,'growth_conditions']
PROT.ANN$DOI = PHENO.META[PROT.ANN$top_pheno_id,'doi']

# Write output tables
write.table(PROT.ANN,'FIGURES/GEAs/topgenes_prot.tsv', sep='\t', quote=F, row.names = F, col.names=T)
save(PROT.ANN, file='DATA/ATHALIANA/PROTANN.rda')


