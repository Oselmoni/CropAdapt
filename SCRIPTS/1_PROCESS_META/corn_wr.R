library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW/corn_wr/41588_2022_1184_MOESM3_ESM(1).xlsx', startRow = 2, sheet = 2)


### keep only species of interest
meta = meta[meta$New.taxonomyb%in%c('Zea mays subsp. mexicana'     ,  'Zea mays subsp. parviglumis'),]

### keep column of interest
meta = meta[,c('Accessiona','Longitudea','Latitudea')]
colnames(meta) = c('ID','LON','LAT')

meta$LON = as.numeric(meta$LON)
meta$LAT = as.numeric(meta$LAT)


### fix date
meta$YEAR = NA


### save output
write.csv(meta, 'DATA/META/corn_wr.csv', row.names = F)

