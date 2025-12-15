library(terra)
library(ecmwfr)
library(rnaturalearth)

lf = list.files('DATA/GEA_INPUT/', pattern='.rda')

land=ne_countries()


for (f in lf) {
  
  print(f)
  d=substr(f, 1 ,nchar(f)-4)
  
  ### load vcf and metadata
  load(paste0('DATA/GEA_INPUT/',f))
  
  ### 
  #### Get copernicus data
  ###
  wf_set_key(read.table('CopernicusAPIkey')[1,1], user=read.table('CopernicusAPIkey')[2,1])
  
  wf_get_key()
  
  # translated using the RStudio IDE addin)
  request <- list(
    dataset_short_name = "sis-biodiversity-era5-global",
    data_format = "netcdf_legacy",
    download_format = "unarchived",
    area = c(60, -20, 33, 20),
    target = "test.nc"
  )
  
  # initiate a download
  file <- wf_request(
    request,
    path = "./" ,
    transfer = TRUE
  )
  
  
  ?wf_set_key
  # Read in data using the `terra`
  # geospatial library
  library(terra)
  r <- rast(file)
  
  