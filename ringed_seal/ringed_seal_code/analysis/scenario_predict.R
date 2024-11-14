

################### LOAD DATA ################################################################

# extract posterior samples of relevant parameters
samples <- rstan::extract(model, pars=c('N', 'w', 'w_pup', 'alpha', 'd', 'f',
                                        'b_max', 'b_min', 'beta', 'theta_0', 'theta_1',
                                        'phi', 'E_tot_fi', 'E_tot_sw_spring', 'E_tot_sw_fall', 
                                        'sigma_E_fi', 'sigma_E_sw_spring', 'sigma_E_sw_fall',
                                        'psi_fi', 'psi_sw_spring', 'psi_sw_fall', 'psi_bycatch',
                                        'pi_c', 'pi_s', 'kappa', 'r'))

# years to include in future simulation
years.sim <- 2023:2038

# Scenarios for Finnish hunting quotas
Q_fi.scenarios <- matrix(NA, 4, length(years.sim))
Q_fi.scenarios[,1] <- 375 # known Finnish quota for 2023
Q_fi.scenarios[1,] <- Q_fi.scenarios[1,1]
Q_fi.scenarios[2,] <- Q_fi.scenarios[2,1] + c(0:(length(years.sim)-1))*15
Q_fi.scenarios[3,] <- Q_fi.scenarios[3,1] + c(0:(length(years.sim)-1))*35
Q_fi.scenarios[4,2:ncol(Q_fi.scenarios)] <- Q_fi.scenarios[4,1] - 270

# Scenarios for Swedish hunting quotas
Q_sw.scenarios <- matrix(NA, 4, length(years.sim))
Q_sw.scenarios[,1] <- 350 # known Swedish quota for 2023
Q_sw.scenarios[1,] <- Q_sw.scenarios[1,1]
Q_sw.scenarios[2,] <- Q_sw.scenarios[2,1] + c(0:(length(years.sim)-1))*15
Q_sw.scenarios[3,] <- Q_sw.scenarios[3,1] + c(0:(length(years.sim)-1))*35
Q_sw.scenarios[4,2:ncol(Q_sw.scenarios)] <- Q_sw.scenarios[4,1] - 270

idx <- 1:nrow(samples$phi) # indices for posterior samples utilized in the simulation

N.sim <- array(NA, dim=c(nrow(Q_fi.scenarios), length(idx), length(years.sim))) # initialize array for storing population sizes

for(i in 1:length(idx)) {
  # extract posterior samples for parameters
  for(j in 1:nrow(Q_fi.scenarios)) {
    
    # extract posterior sample for all parameters
    params.i <- lapply(samples, function(x) {
      if(ncol(as.matrix(x)) > 1) {x[idx[i],]} else {x[idx[i]]}
    })
    
    # sample future ice extent from historical values
    ice.sim <- c(data$ice[length(data$ice)], sample(data$ice, length(years.sim)-1, replace = T))
    
    # run simulation
    results <- simulate(n0 = round(samples$N[idx[i],,length(data$years)+1]), years.sim = years.sim, 
                                 Q_fi.sim = Q_fi.scenarios[j,], Q_sw.sim = Q_sw.scenarios[j,], 
                                 ice.sim = ice.sim, params = params.i, data = data, indices = indices)
    
    N.sim[j,i,] <- colSums(results) 
  }
}

# est. population sizes for last 6 years
N.past <- apply(samples$N[idx,,(length(data$years)-5):(length(data$years)+1)], c(1,3), sum)

# calculate simulated growth rates
growth.rates <- array(NA, dim=c(nrow(Q_fi.scenarios), length(years.sim), 3))
for(i in 1:nrow(Q_fi.scenarios)) {
  N <- cbind(N.past[,ncol(N.past)-1], N.sim[i,,])
  growth.rates[i,,] <- 100*t(apply(N[,2:ncol(N)] / N[,1:(ncol(N)-1)] - 1, 2, quantile, quantiles))
}
growth.rates.df <- data.frame(Scenario = rep(as.factor(c('No change',
                                                         '+30 /yr',
                                                         '+70 /yr',
                                                         '-540 in 2024')), each=length(years.sim)), 
                              year = rep(years.sim, nrow(Q_fi.scenarios)),
                              lb = c(t(growth.rates[,,1])), 
                              md = c(t(growth.rates[,,2])), 
                              ub = c(t(growth.rates[,,3])))

# calculate historical estimated growth rates
growth.rates.past <- 100*t(apply(N.past[,2:ncol(N.past)]/N.past[,1:(ncol(N.past)-1)]-1, 
                           2, quantile, quantiles))
growth.rates.past.df <- data.frame(Scenario = rep(as.factor('Historical'), nrow(growth.rates.past)), 
                                   year = (2024-nrow(growth.rates.past)):2023,
                                   lb = c(growth.rates.past[,1]), 
                                   md = c(growth.rates.past[,2]), 
                                   ub = c(growth.rates.past[,3]))


# plot historical and simulated growth rates
legend.title <- 'Quota scenario'
scenario.colors <- c('Historical'='black', 
                     'No change'='#009E73', 
                     '+30 /yr'='#0072B2', 
                     '+70 /yr'='#D55E00', 
                     '-540 in 2024'='purple')
scenario.alphas <- c('Historical'=0.2, 
                     'No change'=0.1, 
                     '+30 /yr'=0.1, 
                     '+70 /yr'=0.1, 
                     '-540 in 2024'=0.1)
breaks <- c('Historical', '-540 in 2024', 'No change', '+30 /yr', '+70 /yr')
ggplot(rbind(growth.rates.past.df, growth.rates.df), aes(x=year, group=Scenario, fill=Scenario, alpha=Scenario)) +
  geom_hline(aes(yintercept=0), linetype='dashed') +
  geom_ribbon(aes(ymin=lb, ymax=ub)) +
  geom_line(aes(y=md, color=Scenario), linewidth=1, alpha=1) +
  scale_fill_manual(legend.title, values=scenario.colors, breaks=breaks) +
  scale_color_manual(legend.title, values=scenario.colors, breaks=breaks) +
  scale_alpha_manual(legend.title, values=scenario.alphas, breaks=breaks) +
  scale_x_continuous(breaks=seq(2018,2038,3)) +
  scale_y_continuous(limits = c(-1,8), breaks=seq(-1,8,1)) +
  labs(x='Year', y='% growth') +
  guides(fill = guide_legend(title.position = "top"))+
  theme_classic() +
  theme(legend.position = c(0.25,0.175), text = element_text(size=font_size))

# save plot
ggsave(filename='future_preds.png', width=8, height=8, path=paste(output_directory, 'plots/', sep=''))


