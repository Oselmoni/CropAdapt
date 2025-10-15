library(openxlsx)



### load metadata
meta = openxlsx::read.xlsx('DATA/META_RAW//commonbean_la/41588_2019_546_MOESM3_ESM.xlsx', startRow = 2)

### get samples of interest
meta = meta[meta$Type=='L',]



### keep column of interest
meta = meta[,c('ID-SEQ','Longtud','Latitude')]
colnames(meta) = c('ID','LON','LAT')


### Fix longitude
HE = substr(meta$LON,1,1) # get hemisphere info
DEG = as.numeric(substr(meta$LON,2,nchar(meta$LON))) # get decimal degree

DEG[HE=='S'] = -DEG[HE=='S'] # change sign in hemisphere south
meta$LON = DEG



### Fix latitude
HE = substr(meta$LAT,1,1) # get hemisphere info
DEG = as.numeric(substr(meta$LAT,2,nchar(meta$LAT))) # get decimal degree

DEG[HE=='W'] = -DEG[HE=='W'] # change sign in hemisphere west
meta$LAT = DEG


### add info on year
meta$YEAR = NA

### save output
write.csv(meta, 'DATA/META/commonbean_la.csv', row.names = F)

