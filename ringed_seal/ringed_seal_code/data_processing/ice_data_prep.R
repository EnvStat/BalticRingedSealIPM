

# global parameters
source_directory <- paste(data_directory, 'source_data/raw_data/Ice/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

years <- 1988:2023 # years included in the study

# function for extracting ice data from netcdf rasters
extract.ice.cover <- function(year, lat_min) {
  ice <- nc_open(paste(source_directory, year, '_ice_apr.nc', sep = ''))
  lon <- ncvar_get(ice, "lon")
  lat <- ncvar_get(ice, "lat")
  ice_cover <- mean(ncvar_get(ice, "kFmiIceConcentration")[,lat > lat_min] > 0, na.rm=T)
  return(ice_cover)
}

# compute ice cover as (% ice cover)x(area of Bothnian Bay in 000's)
ice.cover <- as.matrix(37*sapply(years, extract.ice.cover, lat_min=63))
rownames(ice.cover) <- years
colnames(ice.cover) <- 'ice_cover'

# save result
write.csv(ice.cover, file = paste(output_directory, 'ice_cover.csv', sep=''))


