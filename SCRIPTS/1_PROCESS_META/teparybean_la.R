library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW/teparybean//tpg220363-sup-0002-tabless1-s10.xlsx', startRow = 4)
meta$Type
### get samples of interest
meta = meta[meta$Type=='Landrace',]


### keep column of interest
meta = meta[,c('Accession.ID','Longitude','Latitude')]
colnames(meta) = c('ID','LON','LAT')

### fix coordinates
meta$LON = as.numeric(meta$LON)
meta$LAT = as.numeric(meta$LAT)

### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/teparybean_la.csv', row.names = F)

