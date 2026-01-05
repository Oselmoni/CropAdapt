### Load orthogroups and store in R
orthogroups = read.table('DATA/RAW/PROTEINS/OrthoFinder/Results_Dec17/Orthogroups/Orthogroups.tsv', sep='\t', header=T)

### save output
save(orthogroups, file='DATA/SNP_ANNOTATION/orthogroups.rda')
