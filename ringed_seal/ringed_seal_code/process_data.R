


# set directories for data, model code and model outputs
data_directory <- 'ringed_seal_data/'
code_directory <- 'ringed_seal_code/'
output_directory <- 'ringed_seal_results/'


# uncomment this line to process raw ice raster files
#library(ncdf4) # package for reading netcdf files (for ice data)

# Process raw data
source(paste0(code_directory, 'data_processing/SW_data_prep.R'), local=T)
source(paste0(code_directory, 'data_processing/FI_data_prep.R'), local=T)

# Prepare model input data
source(paste0(code_directory, 'data_processing/SW_hunt_data_prep.R'), local=T)
source(paste0(code_directory, 'data_processing/FI_hunt_data_prep.R'), local=T)
source(paste0(code_directory, 'data_processing/bycatch_data_prep.R'), local=T)
source(paste0(code_directory, 'data_processing/reproduction_data_prep.R'), local=T)
source(paste0(code_directory, 'data_processing/sampling_bias_data_prep.R'), local=T)

# uncomment this line to process raw ice raster files
#source(paste(code_directory, 'data_processing/ice_data_prep.R', sep=''))


