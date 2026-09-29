### load libraries
library(terra)
library(rnaturalearth)
library(sf)

### Load data
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, sep='\t')
rownames(ds_list)=ds_list$dataset
land = ne_countries( scale = 'small')

## Convert land to EE projection
land_ee <- st_transform(land, crs = 8857)


### set ds colors 
ds_cols = c('tomato_wr'='#b83727', 
            'barley_la'='#b86327', 'barley_wr'='#b86327',
            'commonbean_la'='lightblue','commonbean_wr'='lightblue',
            'corn_la'='#ffdf3d','corn_wr'='#ffdf3d',
            'rice_la'='#696760','rice_wr'='#696760',
            'sorghum_la'='#dbaf79',
            'soybean_la'='#3f914f','soybean_wr'='#3f914f',
            'teparybean_la'='#451c67', 'teparybean_wr'='#451c67'
            )

## split landraces from wild relatives
ds_la = ds_list$dataset[regexpr('_la', ds_list$dataset)!=-1]
ds_wr = ds_list$dataset[regexpr('_wr', ds_list$dataset)!=-1]

## generate shared metadata file for landraces and cwr
meta_la = data.frame()
meta_wr = data.frame()

for (ds in ds_la) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))

  # add id
  meta$DS = ds
  
  ## add to shared container
  meta_la = rbind(meta_la, meta)
  
}

for (ds in ds_wr) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))
  
  # add id
  meta$DS = ds
  
  ## add to shared container
  meta_wr = rbind(meta_wr, meta)
  
}

### Create container of sample sizes
Ns = c()

### Run distribution plots
pdf('FIGURES/samples_map.pdf', width = 6, height = 6)
{
par(mfrow=c(2,1))

### plot landraces
par(mar=c(0,0,1,0))
plot(NA, xlim=ext(land_ee)[1:2], ylim=ext(land_ee)[3:4], axes=F)
plot(land_ee, add=T, col='grey95', border=NA)

for (ds in names(sort(table(meta_la$DS), decreasing = T))) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))

  # Convert to sf object (WGS84)
  sf_df <- st_as_sf(meta, coords = c("LON", "LAT"), crs = 4326)
  
  # Transform to Equal Earth (EPSG:8857)
  sf_df_equal_earth <- st_transform(sf_df, crs = 8857)
  
  # Retrive coord
  coord_ee = st_coordinates(sf_df_equal_earth)
  
  # plot
  points(coord_ee[,1], coord_ee[,2], pch=16, cex=0.25, col=ds_cols[ds])
  
  # store sample size
  Ns[ds] = nrow(meta)
    
}

box()
title(main='A) Landraces datasets', adj=0, cex.main=0.75)
legend('bottomleft', paste0(gsub('(.*) landraces(.*)', '\\1 \\2', ds_list[ds_la,'name']),' (',Ns[ds_la],')'), col=ds_cols[ds_la], cex=0.5, pch=16, box.lty=0, bg=NA)


### Create container of sample sizes
Ns = c()

### plot cwr
par(mar=c(0,0,1,0))
plot(NA, xlim=ext(land_ee)[1:2], ylim=ext(land_ee)[3:4], axes=F)
plot(land_ee, add=T, col='grey95', border=NA)

for (ds in names(sort(table(meta_wr$DS), decreasing = T))) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))
  
  # Convert to sf object (WGS84)
  sf_df <- st_as_sf(meta, coords = c("LON", "LAT"), crs = 4326)
  
  # Transform to Equal Earth (EPSG:8857)
  sf_df_equal_earth <- st_transform(sf_df, crs = 8857)
  
  # Retrive coord
  coord_ee = st_coordinates(sf_df_equal_earth)
  
  # plot
  points(coord_ee[,1], coord_ee[,2], pch=16, cex=0.25, col=ds_cols[ds])
  

  # store sample size
  Ns[ds] = nrow(meta)
  
}
box()
title(main='B) Crop Wild Relatives datasets', adj=0, cex.main=0.75)
legend('bottomleft', paste0(gsub('(.*) wild relatives(.*)', '\\1 \\2', ds_list[ds_wr,'name']),' (',Ns[ds_wr],')'), col=ds_cols[ds_wr], cex=0.5, pch=16, box.lty=0, bg=NA)
}
dev.off()

