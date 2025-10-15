library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW/millet_la/41467_2020_19066_MOESM4_ESM.xlsx', startRow = 2)


colnames(meta)
### keep column of interest
meta = meta[,c('Landraces.accession.number','Longitude','Latitude','Year.of.collect')]
colnames(meta) = c('ID','LON','LAT','YEAR')




### save output
write.csv(meta, 'DATA/META/millet_la.csv', row.names = F)

