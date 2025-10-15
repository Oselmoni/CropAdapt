library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//rice_la/41477_2020_659_MOESM4_ESM(1).xlsx', startRow = 16, sheet = 1)


### keep only 3k samples
meta = meta[meta$SEQ=='Wang_et_al_2018',]
head(meta)
### keep column of interest
meta = meta[,c('ID','LON','LAT')]



### save output
write.csv(meta, 'DATA/META/rice_la.csv', row.names = F)

