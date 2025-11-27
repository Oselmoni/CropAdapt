library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//tomato_wr/abiotic_environmental_data.xlsx', startRow = 1)

# subset metadata
meta = meta[,c(1,3,4)]

### keep column of interest
colnames(meta) = c('ID','LON','LAT')



### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/tomato_wr.csv', row.names = F)

