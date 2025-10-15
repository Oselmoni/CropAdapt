library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//corn_la/EMLP passport information.xlsx', startRow = 1, sheet = 2)


### keep column of interest
meta = meta[,c('Accession_no','Longitude.(oE)','Latitude.(oN)','Collection_date')]
colnames(meta) = c('ID','LON','LAT','YEAR')




### save output
write.csv(meta, 'DATA/META/corn_la.csv', row.names = F)

