

# set directories for data, model code and model outputs
data_directory <- 'ringed_seal_data/'
code_directory <- 'ringed_seal_code/'
output_directory <- 'ringed_seal_results/'


model_data_directory <- paste0(data_directory, 'model_data/')

# load packages and functions
library(rstan)
library(ggplot2)
library(gghalves)
source(paste0(code_directory, 'functions/load_data_function.R'), local=T)
source(paste0(code_directory, 'functions/analysis_functions.R'), local=T)
source(paste0(code_directory, 'functions/predictive_functions.R'), local=T)
source(paste0(code_directory, 'functions/plotting_functions.R'), local=T)


# load data
data <- load_data(model_data_directory)
data$y_survey[data$y_survey == 0] <- NA
indices <- load_indices(data)

# load posterior samples
model <- readRDS(paste0(output_directory, 'ringed_seal_post.rds'))
samples <- rstan::extract(model)

# global parameters
quantiles <- c(0.025, 0.5, 0.975)
font_size = 25

compute_loo <- TRUE

source(paste0(code_directory, 'analysis/generate_outputs.R'), local=T) # extract main estimates and plots
source(paste0(code_directory, 'analysis/post_pred_checks.R'), local=T) # generate plots for posterior predictive checks
source(paste0(code_directory, 'analysis/PSIS.R'), local=T) # Pareto-smoothed Importance Sampling
source(paste0(code_directory, 'analysis/critical_hunting.R'), local=T) # calculate critical hunting level and output plot
source(paste0(code_directory, 'analysis/scenario_predict.R'), local=T) # run scenario-based simulations and output plot


