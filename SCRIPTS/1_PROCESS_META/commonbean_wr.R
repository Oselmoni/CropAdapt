library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW/commonbean_wr/Table 1.xlsx', startRow = 1)

### keep column of interest
meta = meta[,c('Accesion','Longitude.[decimal]','Latitude.[decimal]','Date')]
colnames(meta) = c('ID','LON','LAT','YEAR')

### fix coordinates
meta$LON = as.numeric(meta$LON)
meta$LAT = as.numeric(meta$LAT)

### fix date
meta$YEAR = as.numeric(substr(meta$YEAR, nchar(meta$YEAR)-3, nchar(meta$YEAR)) )


### save output
write.csv(meta, 'DATA/META/commonbean_wr.csv', row.names = F)

