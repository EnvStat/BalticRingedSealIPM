
library(ggplot2)

# set directories for data, model code and model outputs
data_directory <- 'ringed_seal_data/'
code_directory <- 'ringed_seal_code/'
output_directory <- 'ringed_seal_results/'

# import required functions
source(paste0(code_directory, 'functions/analysis_functions.R'), local=T)
source(paste0(code_directory, 'functions/load_data_function.R'), local=T)

# load data
data <- load_data(paste0(data_directory, 'model_data/'))
indices <- load_indices(data)
hunting_times_fi <- read.csv(paste0(data_directory,'source_data/hunting_times_samples_fi.csv'))[,'x']
hunting_times_fi <- hunting_times_fi[which(hunting_times_fi > 1.5/12 & hunting_times_fi < 3/12)]
hunting_times_sw <- read.csv(paste0(data_directory,'source_data/hunting_times_samples_sw.csv'))[,'x']

# load posterior samples
model <- readRDS(paste0(output_directory, 'ringed_seal_post.rds'))
samples <- rstan::extract(model)


# compute hunting rate throughout the hunting season
hunting_bag <- function(t, N0, mu, E, w, ice) {
  n_eff = w*N0/ice
  dH = exp(-E%*%(n_eff*(1-exp(-mu*t))/mu))*(E%*%(n_eff*(mu*exp(-mu*t))/mu))
  return(dH)
}


hdist_fi <- hist(hunting_times_fi, breaks=15, plot = F)
dh_fi <- hdist_fi$counts / sum(hdist_fi$counts)
t_fi <- hdist_fi$mids
dt_fi <- t_fi[2] - t_fi[1]

hdist_sw <- hist(hunting_times_sw, breaks=35, plot = F)
dh_sw <- hdist_sw$counts / sum(hdist_sw$counts)
t_sw <- hdist_sw$mids
t_sw_spring <- t_sw[which(t_sw > 2/12 & t_sw < 4/12)]
t_sw_fall <- t_sw[which(t_sw > 6/12 & t_sw < 8/12)]
dt_sw <- t_sw[2] - t_sw[1]

ice <- mean(data$ice[which(1988:2023 >= 2015)]) # use mean ice extent since 2015

# compute model predicted hunting distributions

dh_fi_pred <- matrix(NA, nrow(samples$n0), length(t_fi))
dh_sw_spring_pred <- matrix(NA, nrow(samples$n0), length(t_sw_spring))
dh_sw_fall_pred <- matrix(NA, nrow(samples$n0), length(t_sw_fall))
for(i in 1:nrow(samples$n0)) {
  # Finnish hunting
  dh_fi_pred[i,] <- sapply(t_fi-t_fi[1], hunting_bag, 
                           N0 = rowMeans(samples$N[i,,which(1988:2023 >= 2015)]), 
                           mu = samples$mu[i,],
                           E = samples$E_fi[i,],
                           w = samples$w_vec[i,],
                           ice = ice)*dt_fi
  dh_fi_pred[i,] <- dh_fi_pred[i,] / sum(dh_fi_pred[i,])
  
  # Swedish hunting (spring + fall)
  dh_sw_spring_pred[i,] <- sapply(t_sw_spring-t_sw_spring[1], hunting_bag,
                           N0 = rowMeans(samples$N[i,,which(1988:2023 >= 2015)]),
                           mu = samples$mu[i,],
                           E = samples$E_sw_spring[i,],
                           w = 1,
                           ice = 1)*dt_sw
  dh_sw_fall_pred[i,] <- sapply(t_sw_fall-t_sw_fall[1], hunting_bag,
                                  N0 = rowMeans(samples$N[i,,which(1988:2023 >= 2015)]),
                                  mu = samples$mu[i,],
                                  E = samples$E_sw_fall[i,],
                                  w = 1,
                                  ice = 1)*(1-sum(dh_sw_spring_pred[i,]))*dt_sw
  h_sw_tot <- sum(c(dh_sw_spring_pred[i,], dh_sw_fall_pred[i,]))
  dh_sw_spring_pred[i,] <- dh_sw_spring_pred[i,] / h_sw_tot
  dh_sw_fall_pred[i,] <- dh_sw_fall_pred[i,] / h_sw_tot
}


############ plot results ###########

font_size <- 25

ggplot() +
  geom_ribbon(aes(x = t_fi, 
                  ymin = apply(dh_fi_pred, 2, quantile, 0.025),
                  ymax = apply(dh_fi_pred, 2, quantile, 0.975)), alpha = 0.25, fill='blue') +
  geom_line(aes(x = t_fi, y = apply(dh_fi_pred, 2, median))) +
  geom_col(aes(x = t_fi, y = dh_fi), alpha=0.65) +
  labs(x = 'Time', y='Proportion of hunting in Finland') +
  theme_classic() +
  theme(text = element_text(size=font_size))
ggsave(filename='hunting_dist_fi.png', width=8, height=8, path=paste0(output_directory, 'plots/'))

ggplot() +
  geom_ribbon(aes(x = t_sw_spring, 
                  ymin = apply(dh_sw_spring_pred, 2, quantile, 0.025),
                  ymax = apply(dh_sw_spring_pred, 2, quantile, 0.975)), alpha = 0.25, fill='blue') +
  geom_line(aes(x = t_sw_spring, y = apply(dh_sw_spring_pred, 2, median))) +
  geom_ribbon(aes(x = t_sw_fall, 
                  ymin = apply(dh_sw_fall_pred, 2, quantile, 0.025),
                  ymax = apply(dh_sw_fall_pred, 2, quantile, 0.975)), alpha = 0.25, fill='blue') +
  geom_line(aes(x = t_sw_fall, y = apply(dh_sw_fall_pred, 2, median))) +
  geom_col(aes(x = t_sw, y = dh_sw), alpha=0.65) +
  labs(x = 'Time', y='Proportion of hunting in Sweden') +
  theme_classic() +
  theme(text = element_text(size=font_size))
ggsave(filename='hunting_dist_sw.png', width=8, height=8, path=paste0(output_directory, 'plots/'))

