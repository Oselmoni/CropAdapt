library(terra)
library(ecmwfr)
library(rnaturalearth)
library(elevatr)

lf = list.files('DATA/GEA_INPUT/meta/', pattern='.rda')

land=ne_countries(scale='large')

### load list of env dataset ids
envdata = read.csv('DATA/ENV/envlist.csv', row.names=3)

for (f in lf) {

  #f=lf[1]
  print(f)
  
  d=substr(f, 1 ,nchar(f)-4)
  
  ### load vcf and metadata
  load(paste0('DATA/GEA_INPUT/meta/',f))
  

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
  
  ### save extracted env data
  save(meta, file=paste0('DATA/GEA_INPUT/meta_env//',f))
    
}
 

