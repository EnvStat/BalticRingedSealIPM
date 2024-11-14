
load_data <- function(data_directory) {
  
  # load aerial survey data
  y_survey_raw <- read.csv(file = paste(data_directory, 'aerial_survey.csv', sep=''), header = T)
  y_survey <- y_survey_raw$survey_estimate
  y_survey[is.na(y_survey)] <- 0

  years <- 1988:2022 # years included in the model
  years_std <- (years - mean(years)) / sd(years)
  
  # load ice data
  ice_raw <- read.csv(file = paste(data_directory, 'ice_cover.csv', sep=''), header = T)
  ice <- ice_raw$ice_cover

  # load Finnish hunting bag and quotas
  y_hb_fi_raw <- read.csv(file = paste(data_directory, 'hunting_bag_fi.csv', sep=''), header = T)
  y_hb_fi <- y_hb_fi_raw$bag
  Q_fi <- y_hb_fi_raw$quota
  
  # load Swedish hunting bag and quotas
  y_hb_sw_raw <- read.csv(file = paste(data_directory, 'hunting_bag_sw.csv', sep=''), header = T)
  y_hb_sw_spring <- y_hb_sw_raw$spring_bag
  y_hb_sw_fall <- y_hb_sw_raw$fall_bag
  y_hb_sw <- y_hb_sw_spring + y_hb_sw_fall
  Q_sw <- y_hb_sw_raw$quota
  
  # load Finnish hunting seal samples
  y_hs_fi_raw <- read.csv(file = paste(data_directory, 'hunting_samples_fi.csv', sep=''), header = T)
  y_hs_fi_ageNA_raw <- read.csv(file = paste(data_directory, 'hunting_samples_fi_ageNA.csv', sep=''), header = T)
  y_hs_fi_sexNA_raw <- read.csv(file = paste(data_directory, 'hunting_samples_fi_sexNA.csv', sep=''), header = T)
  y_hs_fi <- t(y_hs_fi_raw[,-1])
  y_hs_fi_ageNA = t(y_hs_fi_ageNA_raw[,-1]) 
  y_hs_fi_sexNA = t(y_hs_fi_sexNA_raw[,-1])
  
  # load Swedish hunting samples from spring season
  y_hs_sw_spring_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_spring.csv', sep=''), header = T)
  y_hs_sw_ageNA_spring_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_spring_ageNA.csv', sep=''), header = T)
  y_hs_sw_sexNA_spring_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_spring_sexNA.csv', sep=''), header = T)
  y_hs_sw_spring <- t(y_hs_sw_spring_raw[,-1])
  y_hs_sw_ageNA_spring = t(y_hs_sw_ageNA_spring_raw[,-1]) 
  y_hs_sw_sexNA_spring = t(y_hs_sw_sexNA_spring_raw[,-1])
  
  # load Swedish hunting samples from fall season
  y_hs_sw_fall_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_fall.csv', sep=''), header = T)
  y_hs_sw_ageNA_fall_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_fall_ageNA.csv', sep=''), header = T)
  y_hs_sw_sexNA_fall_raw <- read.csv(file = paste(data_directory, 'hunting_samples_sw_fall_sexNA.csv', sep=''), header = T)
  y_hs_sw_fall <- t(y_hs_sw_fall_raw[,-1])
  y_hs_sw_ageNA_fall = t(y_hs_sw_ageNA_fall_raw[,-1]) 
  y_hs_sw_sexNA_fall = t(y_hs_sw_sexNA_fall_raw[,-1])
  
  # load sampling bias data for Swedish hunting in spring
  l <- read.csv(file = paste(data_directory, 'length_by_class.csv', sep=''), header = T)
  l <- as.matrix(l[,-1])
  o_16_19 <- read.csv(file = paste(data_directory, 'sampling_by_length_16_19.csv', sep=''), header = T)
  o_16_19 <- as.matrix(o_16_19[,-1])
  o_20_21 <- read.csv(file = paste(data_directory, 'sampling_by_length_20_21.csv', sep=''), header = T)
  o_20_21 <- as.matrix(o_20_21[,-1])
  
  # load mean occurrence times of hunting
  tau_fi_raw <- read.csv(file = paste(data_directory, 'hunting_times_fi.csv', sep=''), header = T)
  tau_sw_raw <- read.csv(file = paste(data_directory, 'hunting_times_sw.csv', sep=''), header = T)
  tau_fi <- c(tau_fi_raw[,-1])
  tau_sw_spring <- c(tau_sw_raw$spring)
  tau_sw_fall <- c(tau_sw_raw$fall)

  # load bycaught samples
  y_bc_raw <- read.csv(file = paste(data_directory, 'bycatch_samples.csv', sep=''), header = T)
  y_bc_ageNA_raw <- read.csv(file = paste(data_directory, 'bycatch_samples_ageNA.csv', sep=''), header = T)
  y_bc_sexNA_raw <- read.csv(file = paste(data_directory, 'bycatch_samples_sexNA.csv', sep=''), header = T)
  y_bc <- t(y_bc_raw[,-1])
  y_bc_ageNA = t(y_bc_ageNA_raw[,-1]) 
  y_bc_sexNA = t(y_bc_sexNA_raw[,-1])
  
  # load pregnancy data
  y_pregnant_raw <- read.csv(file = paste(data_directory, 'pregnancies.csv', sep=''), header = T)
  pregnant_ss <- y_pregnant_raw$sample_size
  y_pregnant <- y_pregnant_raw$pregnant
  
  # load joint assessments of placental scars and CA
  z_raw <- read.csv(file = paste(data_directory, 'reproductive_assessments.csv', sep=''), header = T)
  z <- t(z_raw[,-1])
  
  # load placental scar data
  y_scar_raw <- read.csv(file = paste(data_directory, 'scars.csv', sep=''), header = T)
  scar_ss <- y_scar_raw$sample_size
  y_scar <- y_scar_raw$scars
  
  # load CA data
  y_CA_raw <- read.csv(file = paste(data_directory, 'CA.csv', sep=''), header = T)
  CA_ss <- y_CA_raw$sample_size
  y_CA <- y_CA_raw$CA
  
  # covariance matrix for logistic-Gaussian bias priors (see Appendix  on prior distributions)
  D <- matrix(1e100, 12, 12)
  D[2:5, 2:5] <- D[8:11, 8:11] <- as.matrix(dist(1:4))
  D[2:5, 8:11] <- D[8:11, 2:5] <- 1 + as.matrix(dist(1:4))
  D[1, 7] <- D[7, 1] <- D[6, 12] <- D[12, 6] <- 1
  D <- D*(matrix(1, 12, 12) - diag(1, 12))
  rho <- 0.9887^D
  
  # compute conditional covariance matrix using adult males as the reference level
  sigma <- 5
  K <- diag(sigma,12)%*%rho%*%diag(sigma,12)
  K <- K[1:11,1:11] - (K[1:11,12]%*%t(K[12,1:11]))/sigma^2
  L_bias <- t(chol(K))
  
  # prior bounds on baseline haul-out %
  w_bounds = c(0.4, 0.9)
  
  # initial guess for demographic structure in 1988 (used to initialize eigenvector computation)
  phi <- rep(c(0.65, 0.89, 0.89, 0.89, 0.89, 0.95), 2) # guess for survival probabilities
  N0 <- stable.structure.eq(phi=phi/lambda(phi, b=0.28)) # guess for initial age structure
  
  # datapoint exclusions for sensitivity
  exclude_survey <- rep(0, length(y_survey))
  exclude_rep <- rep(0, length(years))
  exclude_hb <- rep(0, length(years))
  exclude_samples <- rep(0, length(years))

  return(list(years=years, 
              years_std = years_std, 
              ice = ice,
              
              y_survey = y_survey,
              
              Q_fi = Q_fi, 
              Q_sw = Q_sw, 
              
              y_hb_fi = y_hb_fi, 
              y_hb_sw_spring = y_hb_sw_spring, 
              y_hb_sw_fall = y_hb_sw_fall,
              y_hb_sw = y_hb_sw,
              
              tau_fi = tau_fi, 
              tau_sw_spring = tau_sw_spring, 
              tau_sw_fall = tau_sw_fall,
              
              y_hs_fi = y_hs_fi, 
              y_hs_fi_sexNA = y_hs_fi_sexNA, 
              y_hs_fi_ageNA = y_hs_fi_ageNA, 
              
              y_hs_sw_spring = y_hs_sw_spring, 
              y_hs_sw_fall = y_hs_sw_fall,
              y_hs_sw_ageNA_spring = y_hs_sw_ageNA_spring, 
              y_hs_sw_ageNA_fall = y_hs_sw_ageNA_fall,
              y_hs_sw_sexNA_spring = y_hs_sw_sexNA_spring, 
              y_hs_sw_sexNA_fall = y_hs_sw_sexNA_fall,
              l = l,
              o_16_19 = o_16_19,
              o_20_21 = o_20_21,

              y_bc = y_bc,
              y_bc_sexNA = y_bc_sexNA, 
              y_bc_ageNA = y_bc_ageNA,

              pregnant_ss = pregnant_ss, 
              y_pregnant = y_pregnant,
              z = z,
              CA_ss = CA_ss, 
              y_CA = y_CA, 
              scar_ss = scar_ss, 
              y_scar = y_scar,
              
              L_bias = L_bias,
              w_bounds = w_bounds,
              
              N0 = N0,
              
              exclude_survey = exclude_survey,
              exclude_rep = exclude_rep,
              exclude_hb = exclude_hb,
              exclude_samples = exclude_samples)
         )
}

load_indices <- function(data) {
  tau_basking <- 3/12 # length of the nursing+basking period (1 March - 1 June)
  tau_foraging <- 7/12 # length of the foraging period (1 June - 1 January)
  tau_subnivean <- 1 - tau_basking - tau_foraging # length of the subnivean period (1 January - 1 March)
  tau_0_fi <- 1.5/12 # time from pupping to onset of FI hunting (1 March - 14 April)
  tau_0_sw_spring <- 2/12 # time from pupping to onset of SW spring hunting season (March 1 - May 1)
  tau_0_sw_fall <- 6/12 # time from pupping to onset of SW spring hunting season (March 1 - Sept 1)
  tau_p <- 1/2 # average time between mating (early April) and sampling of pregnant females (mid-September)
  
  return(list(a = 6, # number of age classes
              k = 4, # number of mortality sources (FI hunting + SW hunting + natural)
              t = length(data$years), # total number of years in the study period
              t_h_sw = 2022-2014, # number of years with legal SW hunting
              t_h_fi = 2022-2015, # number of years with legal FI hunting
              tau_basking = tau_basking,
              tau_foraging = tau_foraging,
              tau_subnivean = tau_subnivean,
              tau_0_fi = tau_0_fi,
              tau_0_sw_spring = tau_0_sw_spring,
              tau_0_sw_fall = tau_0_sw_fall,
              tau_p = tau_p
              ))
}




