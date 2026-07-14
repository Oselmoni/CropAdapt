library(terra)
library(ecmwfr)
library(rnaturalearth)
library(elevatr)
library(openxlsx)
library(missForest)


### load sampels metadata
META = read.xlsx('DATA/ATHALIANA/RAW/SupData_from_ExpositoAlonso_2019.xlsx', sheet = 2, startRow = 3)
rownames(META) = META$id

land=ne_countries(scale='large')

### load list of env dataset ids
envdata = read.csv('DATA/ENV/envlist.csv', row.names=3)


meta = META[,c('id','longitude','latitude')]

  
  ### 
  #### Get point extraction of every variable
  ###
  
  for (evar in rownames(envdata)) {
    
    ### load raster of variable of interest
    if (envdata[evar,'Source'] %in% c('SoilGrids','EarthEnv')) { # if input variable is fro SoilGrids or EarthEnv... need to use processed versions
      RAST = rast(paste0(envdata[evar,'folder'],'_Processed/',envdata[evar,'variableRaw']))
    } else {
      RAST = rast(paste0(envdata[evar,'folder'],'/',envdata[evar,'variableRaw']))
    }
    
    ### extract point values
    meta[,evar] = extract(RAST, meta[,2:3])[,2]
    
  }
  
### impute missing data
imeta = missForest(meta)

### imputed metadata matrix
meta = imeta$ximp

save(meta, file='DATA/ATHALIANA/meta.rda')
  


