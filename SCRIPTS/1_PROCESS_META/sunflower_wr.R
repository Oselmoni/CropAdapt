library(xlsx)

### load metadata
meta = xlsx::read.xlsx('DATA/META_RAW//sunflower_wr/41586_2020_2467_MOESM3_ESM(2).xlsx', sheetIndex = 1, startRow = 1)


### get samples of interest
meta = meta[which(meta$Taxon=='Helianthus annuus'),]



### keep column of interest, expand to sample identifiers
meta = do.call(rbind, by(meta, meta$Population.ID, function(x) {
  
  ID = paste0('ANN',sprintf('%04d', as.numeric(substr(x$Individuals,4,7)):as.numeric(substr(x$Individuals,9,12))))

  return(data.frame(ID, 'LON'=x$Longitude, 'LAT'=x$Latitude, 'YEAR'=x$Collected))  
  
}))


### add info on year
meta$YEAR = as.numeric(substr(meta$YEAR, 1,4))



### save output
write.csv(meta, 'DATA/META/sunflower_wr.csv', row.names = F)
