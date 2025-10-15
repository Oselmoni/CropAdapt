### load metadata
meta = xlsx::read.xlsx('DATA/META_RAW//rice_wr/TableS2.xls', sheetIndex = 1, header = T, startRow = 2)




### keep column of interest
meta = meta[,c('Accession.ID','Longitude','Latitude')]
colnames(meta) = c('ID','LON','LAT')



### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/rice_wr.csv', row.names = F)

