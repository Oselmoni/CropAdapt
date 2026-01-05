library(terra)


#### Reproject SoilGrids data
lf = list.files('DATA/ENV/SoilGrids/')

for (f in lf) {
  
  print(f)
  
  input = paste0('DATA/ENV/SoilGrids/',f)
  output = paste0('DATA/ENV/SoilGrids_Processed/',f)
  
  project(rast(input), y='epsg:4326', filename=output)

}



#### Reproject SoilGrids data
lf = list.files('DATA/ENV/SoilGrids/')

# for every soilgrid layer...
for (f in lf) {
  
  print(f)
  
  input = paste0('DATA/ENV/SoilGrids/',f)
  output = paste0('DATA/ENV/SoilGrids_Processed/',f)
  
  ## reproject to 4326
  project(rast(input), y='epsg:4326', filename=output)
  
}



#### Resample Landcover classes of interest
lf = list.files('DATA//ENV/CONSENSUS_LC//')

## for every landcover layer
for (f in lf) {
  
  print(f)
  
  input = paste0('DATA/ENV/CONSENSUS_LC//',f)
  output = paste0('DATA/ENV/CONSENSUS_LC_Processed//',f)
  
  RI = rast(input)
  
  # create canvas with downsampled resolution
  CANVAS = RI;res(CANVAS) = c(1,1)

  ## resample to 100 km
  resample(RI, y=CANVAS, filename=output, overwrite=T)

}
