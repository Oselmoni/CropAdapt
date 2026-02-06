### Load orthogroups and store in R
orthogroups = read.table('DATA/ORTHOFINDER/Orthogroups/Orthogroups.tsv', sep='\t', header=T)

### save output
save(orthogroups, file='DATA/SNP_ANNOTATION/orthogroups.rda')


### For every dataset -> create dataset of position of orthogroups along genome

## load dataset ID
selectedDS = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, row.names=1)


for (ds in rownames(selectedDS)) {
  print(ds)
  
  ### load CDS
  load(paste0('DATA/GENE_ANNOTATION/',selectedDS[ds,"annotation_file"],'.rda'))
  
  # find ortogroup corresponding to gene annotation file
  species_OG = orthogroups[,selectedDS[ds,"protein_file"]]
  names(species_OG) = orthogroups$Orthogroup
  
  # remove empty ortogroups
  species_OG = species_OG[species_OG!='']
  
  # map genes to ortogroups
  GENE_TO_OG = rep(names(species_OG), lapply(species_OG, function(x) {length(unlist(strsplit(x, ', ')))}))
  names(GENE_TO_OG) = unlist(lapply(species_OG, strsplit, ', '))
  
  # get midpoint position of every gene
  genes_MP = do.call(rbind , by(CDS, CDS$GENEID, function(x) {

    mi = min(x$STA)
    ma = max(x$END)
    mid = mi+(ma-mi)/2
    
    return(data.frame('CHR'=unique(x$CHR), 'POS'=mid))
      }))
  
  # find orthologs
  genes_MP$OG = GENE_TO_OG[rownames(genes_MP)]
  

  save(genes_MP, file=paste0('DATA/GENE_ANNOTATION/genes_MP_',ds,'.rda'))  
  
}
