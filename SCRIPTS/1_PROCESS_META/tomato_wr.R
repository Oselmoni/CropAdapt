library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//tomato_wr/mec15477-sup-0002-TableS1.xlsx', startRow = 1)



### keep column of interest
colnames(meta) = c('ID','LON','LAT')



### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/tomato_wr.csv', row.names = F)

