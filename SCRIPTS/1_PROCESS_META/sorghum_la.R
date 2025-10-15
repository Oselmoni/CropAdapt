library(openxlsx)

### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//sorghum_la/Tables S1 to S21.xlsx', sheet = 1, startRow = 3)


## remove empty first line
meta = meta[-1,]

### keep column of interest
meta = meta[,c('is_no','Longitude','Latitude')]
colnames(meta) = c('ID','LON','LAT')

### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/sorghum_la.csv', row.names = F)

