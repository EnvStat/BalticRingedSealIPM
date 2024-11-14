
# main model
rm(list=ls())
source('ringed_seal_code/process_data.R', local=T)
rm(list=ls())
source('ringed_seal_code/fit_model.R', local=T)
rm(list=ls())
source('ringed_seal_code/get_results.R', local=T)

# alternative models
rm(list=ls())
source('ringed_seal_code/sensitivity_analyses/fit_data_sensitivities.R', local=T)
rm(list=ls())
source('ringed_seal_code/sensitivity_analyses/fit_model_sensitivities.R', local=T)
rm(list=ls())
source('ringed_seal_code/sensitivity_analyses/get_sensitivity_results.R', local=T)