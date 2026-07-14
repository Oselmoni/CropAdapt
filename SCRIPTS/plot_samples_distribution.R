### load libraries
library(terra)
library(rnaturalearth)

### Load data
ds_list = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, sep='\t')
rownames(ds_list)=ds_list$dataset
land = ne_countries( scale = 'small')

### set ds colors 
ds_cols = c('tomato_wr'='#b83727', 
            'barley_la'='#dbaf79', 'barley_wr'='#dbaf79',
            'commonbean_la'='#de7eab','commonbean_wr'='#b55783',
            'corn_la'='#ffdf3d','corn_wr'='#ffdf3d',
            'rice_la'='#696760','rice_wr'='#696760',
            'sorghum_la'='#b86327',
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
par(mfrow=c(2,1))

### plot landraces
par(mar=c(0,0,1,0))
plot(NA, xlim=c(-120,140), ylim=c(-40,70), axes=F)
plot(land, add=T, col='grey90', border=NA)

for (ds in names(sort(table(meta_la$DS), decreasing = T))) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))

  # plot
  points(meta$LON, meta$LAT, pch=16, cex=0.25, col=ds_cols[ds])
  
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
plot(NA, xlim=c(-130,140), ylim=c(-40,70), axes=F)
plot(land, add=T, col='grey90', border=NA)

for (ds in names(sort(table(meta_wr$DS), decreasing = T))) {
  
  # load coordinates
  load(paste0('DATA/GEA_INPUT/meta/meta_',ds,'.rda'))
  
  # plot
  points(meta$LON, meta$LAT, pch=16, cex=0.25, col=ds_cols[ds])
  
  # store sample size
  Ns[ds] = nrow(meta)
  
}
box()
title(main='B) Crop Wild Relatives datasets', adj=0, cex.main=0.75)
legend('bottomleft', paste0(gsub('(.*) wild relatives(.*)', '\\1 \\2', ds_list[ds_wr,'name']),' (',Ns[ds_wr],')'), col=ds_cols[ds_wr], cex=0.5, pch=16, box.lty=0, bg=NA)

