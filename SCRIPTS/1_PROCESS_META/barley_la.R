library(xlsx)

### load metadata
meta = xlsx::read.xlsx('DATA/META_RAW//barley/tpj15908-sup-0015-tabless1-s7.xlsx', sheetIndex = 1, startRow = 3)


### get samples of interest
meta = meta[meta$Domestication.status=='Landrace',]


### keep column of interest
meta = meta[,c('Accession.identifier.current.assembly','Longitude','Latitude')]
colnames(meta) = c('ID','LON','LAT')

### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/barley_la.csv', row.names = F)

