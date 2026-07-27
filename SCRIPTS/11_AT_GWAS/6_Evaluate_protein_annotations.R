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


## Write output
write.table(TOPGENES,'FIGURES/GEAs/TOPGENES.tsv', sep='\t', quote=F, row.names = F, col.names=T)






### Load UNIPROT table of AT annotations
UNIPROT = read.table('DATA/ATHALIANA/RAW/uniprotkb_taxonomy_id_3702_2026_07_15.tsv', sep='\t', quote="", header=T, fill = T)
UNIPROT = UNIPROT[order(UNIPROT$Reviewed),] # set reviewed prot on top







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
      
      if (file.exists(paste0('DATA/ATHALIANA/GwP/GwP_',at,'.rda'))) {
        load(paste0('DATA/ATHALIANA/GwP/GwP_',at,'.rda'))
        
        toph =  names(sort(apply(OUT$GxP,1,mean), decreasing = T)[1])
        tophR2 =  sort(apply(OUT$GxP,1,mean), decreasing = T)[1]
        check1 = mean(OUT$RND.PH$dR2>tophR2)
        check2 = mean(OUT$RND.GT$dR2>tophR2)
      
    } else { toph = tophR2 = check1 = check2 = ''}
    
      PROT.ANN = rbind(PROT.ANN, data.frame('OrthogeneID'=og, 'A. thaliana OG'=at, toph, 'R2'=tophR2, 'p1'=check1 , 'p2'=check2, 'UP_ID'=UP_ID, 'UP_Name'=UP_name, 'GO.MF'=MF, 'GO.BP'=BP, 'GO.CC'=CC))
      
   
   
  
}}
}

### add info on phenotypes
load('DATA/ATHALIANA/PHENO_META.rda')

PROT.ANN$Pheno = PHENO.META[PROT.ANN$toph,'name']
PROT.ANN$Study = PHENO.META[PROT.ANN$toph,'study']
PROT.ANN$Scoring = PHENO.META[PROT.ANN$toph,'scoring']
PROT.ANN$Growth.Conditions = PHENO.META[PROT.ANN$toph,'growth_conditions']
PROT.ANN$DOI = PHENO.META[PROT.ANN$toph,'doi']

# Write output tables
write.table(PROT.ANN,'FIGURES/GEAs/topgenes_prot.tsv', sep='\t', quote=F, row.names = F, col.names=T)
