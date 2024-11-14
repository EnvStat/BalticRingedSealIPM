

# set directories for data, model code and model outputs
data_directory <- 'ringed_seal_data/'
code_directory <- 'ringed_seal_code/'
model_data_directory <- paste0(data_directory, 'model_data/')

# load packages and functions
library(rstan)
library(LDATS)
library(ggplot2)
library(gghalves)
source(paste0(code_directory, 'functions/load_data_function.R'), local=T)
source(paste0(code_directory, 'functions/analysis_functions.R'), local=T)
source(paste0(code_directory, 'functions/predictive_functions.R'), local=T)
source(paste0(code_directory, 'functions/plotting_functions.R'), local=T)

# global parameters
quantiles <- c(0.025, 0.5, 0.975)
font_size <- 25

# load data
data <- load_data(model_data_directory)
indices <- load_indices(data)
data$y_survey[data$y_survey == 0] <- NA

# create sensitivity datasets
source(paste0(code_directory, 'functions/create_sensitivity_datasets.R'), local=T)

# prior & posterior comparisons
source(paste0(code_directory, 'analysis/prior_vs_posterior.R'), local=T)

######## simplified model ###############
output_directory <- paste0('ringed_seal_results/sensitivities/m1/')

# load posterior samples
model <- readRDS(paste0(output_directory, 'ringed_seal_post_sens_m1.rds'))
samples <- rstan::extract(model)

compute_loo <- TRUE # whether to output LOO_CV results

source(paste0(code_directory, 'analysis/generate_outputs.R'), local=T) # extract main estimates and plots
source(paste0(code_directory, 'analysis/post_pred_checks.R'), local=T) # generate plots for posterior predictive checks
source(paste0(code_directory, 'analysis/PSIS.R'), local=T) # Pareto-smoothed Importance Sampling

########## prior sensitivities ###########
for(i in 1:2) {
  
  output_directory <- paste0('ringed_seal_results/sensitivities/m',i+1,'/')
  
  # load posterior samples
  model <- readRDS(paste0(output_directory, 'ringed_seal_post_sens_m',i+1,'.rds'))
  samples <- rstan::extract(model)
  
  source(paste0(code_directory, 'analysis/generate_outputs.R'), local=T) # extract main estimates and plots
  source(paste0(code_directory, 'analysis/post_pred_checks.R'), local=T) # generate plots for posterior predictive checks
  source(paste0(code_directory, 'analysis/PSIS.R'), local=T) # Pareto-smoothed Importance Sampling
  
}

########## data sensitivities ###########

for(i in 1:length(sensitivity_data)) {

  output_directory <- paste0('ringed_seal_results/sensitivities/',i,'/')
  
  data <- sensitivity_data[[i]]
  data$y_survey[data$y_survey == 0] <- NA
  
  # load posterior samples
  model <- readRDS(paste0(output_directory, 'ringed_seal_post_sens',i,'.rds'))
  samples <- rstan::extract(model)
  
  # remove chain 1 from sensitivity 7 due to poor mixing
  if(i == 6) {
    samples <- lapply(samples, function(x) {
      if(length(dim(x))==1) {return(x[(dim(x)[1]/4):dim(x)[1]])}
      if(length(dim(x))==2) {return(x[(dim(x)[1]/4):dim(x)[1],])}
      if(length(dim(x))==3) {return(x[(dim(x)[1]/4):dim(x)[1],,])}
      if(length(dim(x))==4) {return(x[(dim(x)[1]/4):dim(x)[1],,,])}
    })
  }
  
  source(paste0(code_directory, 'analysis/generate_outputs.R'), local=T) # extract main estimates and plots
  source(paste0(code_directory, 'analysis/post_pred_checks.R'), local=T) # generate plots for posterior predictive checks
  source(paste0(code_directory, 'analysis/PSIS.R'), local=T) # Pareto-smoothed Importance Sampling

}

