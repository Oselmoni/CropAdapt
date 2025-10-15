library(xlsx)

### load metadata
meta = xlsx::read.xlsx('DATA/META_RAW//soybean_la/12864_2022_8326_MOESM1_ESM(1).xlsx', sheetIndex = 1, startRow = 2)


### load id converter: three identifiers GID (the one used in the SNP matrix), PID (accession code), common name
idtable  = read.table('DATA/META_RAW/soybean_la/AnLab_1.5K.SampleIDs.txt', sep='\t', header=T, quote="", fill=T)

## link GID to PID
GID.PID = idtable$GID
names(GID.PID) = idtable$PI.Number

## link GID to common name
GID.CN = idtable$GID
names(GID.CN) = idtable$Common.Name



### get samples of interest
meta = meta[which(meta$Type=='landrace'),]

### replace IDS
meta$ID = GID.PID[meta$PI]
meta$ID[is.na(meta$ID)] = GID.CN[meta$Common.Name[is.na(meta$ID)]]


### keep column of interest
meta = meta[,c('ID','longitude','latitude')]
colnames(meta) = c('ID','LON','LAT')


meta=na.omit(meta)



### add info on year
meta$YEAR = NA



### save output
write.csv(meta, 'DATA/META/soybean_la.csv', row.names = F)
