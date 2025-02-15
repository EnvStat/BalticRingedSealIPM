
# set directories for data and model code
data_directory <- 'ringed_seal_data/'
code_directory <- 'ringed_seal_code/'


#################################################################################################
##########################     SETUP      #######################################################
#################################################################################################

######################## PARAMETER INITIALIZATION ###########################

# proper initialization for efficient sampling
# initialization may fail if parameters result in population going extinct
generate.inits <- function() {
  return( list(phi_f0_scale = runif(1, 0.5, 1), 
               phi_f5_unscaled = runif(1, 0.8, 1),
               n0_unscaled = rnorm(1, 2.5/0.65, 0.5), 
               u_pup = rnorm(35, 0, 0.1), 
               u_sex = rnorm(35, 0, 0.1), 
               u0 = matrix(rnorm(12*27, 0, 0.1), 12, 27), 
               u1 = matrix(rnorm(36*1, 0, 0.1), 36, 1), 
               u2 = matrix(rnorm(48*7, 0, 0.1), 48, 7), 
               b_max = runif(1, 0.6, 1), 
               b_scale = runif(1, 0.5, 1),
               sqrt_r_inv = 1/sqrt(abs(rnorm(1, 40, 5))),
               pq_fi = rbeta(1, 1, 3), 
               pq_sw_spring = rbeta(1, 1, 3), 
               pq_sw_fall = rbeta(1, 1, 3),
               sigma_E_fi = runif(1, 0.1, 2), 
               sigma_E_sw_spring = runif(1, 0.1, 2),
               sigma_E_sw_fall = runif(1, 0.1, 2),
               g_fi = rnorm(11, 0, 1), 
               g_sw_spring = rnorm(11, 0, 1), 
               g_sw_fall = rnorm(11, 0, 1), 
               g_bycatch = rnorm(11, 0, 1))
  )}


###################### PREPARE DATA AND MODEL INPUTS ############################

# load packages and functions
library(rstan)
source(paste0(code_directory, 'functions/load_data_function.R'), local=T)
source(paste0(code_directory, 'functions/analysis_functions.R'), local=T)

data <- load_data('ringed_seal_data/model_data/')
indices <- load_indices(data)

# sensitivity of prior bounds on w
w_bounds <- matrix(NA, 2, 2)
w_bounds[1,] <- c(0.3, 0.8) # shift 10% left
w_bounds[2,] <- c(0.5, 1) # shift 10% right


###############################################################################################################
###############################     RUN MODEL    ##############################################################
###############################################################################################################

seed <- 1
set.seed(seed)

# optimize settings for Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

iter <- 4000
warmup <- iter/2

##################################################
####### Fit alternative pregnancy loss model #####
##################################################

model <- stan_model(paste0(code_directory, 'sensitivity_analyses/ringed_seal_ploss_decr.stan'))

output_directory <- paste0('ringed_seal_results/sensitivities/m0/')

stan_data <- c(data, indices)

# run Stan model
fit <- sampling(model, data = stan_data, seed=seed,
                iter = iter, warmup=warmup, chains = 4, 
                init = list(generate.inits(), generate.inits(),
                            generate.inits(), generate.inits()),
                chain_id = c(1,2,3,4), control = list(adapt_delta=0.97, max_treedepth = 11))

####### Save Stan output ###########

fit@stanmodel@dso <- new("cxxdso")
saveRDS(fit, file = paste0(output_directory, 'ringed_seal_post_sens_m0.rds'))

##################################################
######### Fit simplified haul-out model ##########
##################################################

model <- stan_model(paste0(code_directory, 'sensitivity_analyses/ringed_seal_simplified.stan'))

output_directory <- paste0('ringed_seal_results/sensitivities/m1/')

stan_data <- c(data, indices)

# run Stan model
fit <- sampling(model, data = stan_data, seed=seed,
                iter = iter, warmup=warmup, chains = 4, 
                init = list(generate.inits(), generate.inits(),
                            generate.inits(), generate.inits()),
                chain_id = c(1,2,3,4), control = list(adapt_delta=0.97, max_treedepth = 11))

####### Save Stan output ###########

fit@stanmodel@dso <- new("cxxdso")
saveRDS(fit, file = paste0(output_directory, 'ringed_seal_post_sens_m1.rds'))

##################################################
######### Fit prior sensitivities ################
##################################################

model <- stan_model(paste0(code_directory, 'ringed_seal.stan'))

for(i in 1:nrow(w_bounds)) {
  
  output_directory <- paste0('ringed_seal_results/sensitivities/m',i+1,'/')
  
  # data input for Stan model
  data$w_bounds <- w_bounds[i,]
  stan_data <- c(data, indices)
  
  # run Stan model
  fit <- sampling(model, data = stan_data, seed=seed,
                  iter = iter, warmup=warmup, chains = 4, 
                  init = list(generate.inits(), generate.inits(),
                              generate.inits(), generate.inits()),
                  chain_id = c(1,2,3,4), control = list(adapt_delta=0.97, max_treedepth = 11))
  
  ####### Save Stan output ###########
  
  fit@stanmodel@dso <- new("cxxdso")
  saveRDS(fit, file = paste0(output_directory, 'ringed_seal_post_sens_m', i+1,'.rds'))

}







