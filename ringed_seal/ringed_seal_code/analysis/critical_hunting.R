

#################### Functions ###################################################

# compute expected normalized composition of the hunting bag
hunting_bag_prop <- function(N0, mu, E, tau_mid, tau, w, ice) {
  n_eff = w*N0/ice
  p_H = E*n_eff*exp(-mu*tau_mid)
  p_H = p_H / sum(p_H)
  return(p_H)
}

# compute left hand side of Euler-Lotka equation when lambda=1 (i.e. R0), including hunting mortality
R0 <- function(H, N, params) {
  with(params, {
    b <- b.N(b_max, theta_0, theta_1, N) # compute density dependent birth rate
    u <- stable.structure.eq(phi) # compute stable stage distribution
    N0 <- N*u # compute population vector
    
    # probabilities of surviving to beginning of hunting season
    phi_tau_0_fi <- phi^tau_0_fi # FI
    phi_tau_0_sw_spring <- phi^tau_0_sw_spring # SW spring
    phi_tau_0_sw_fall <- phi^tau_0_sw_fall # SW fall
    
    # compute expected hunting bag composition in Finland, in Sweden during spring and during fall
    H_fi <- (p_fi*H) * hunting_bag_prop(N0=phi_tau_0_fi*N0, mu=-log(phi), E=E_fi, 
                                        tau_mid=tau_fi-tau_0_fi, w=w_vec, ice=ice)
    H_sw_spring <- ((1-p_fi)*p_sw_spring*H) * hunting_bag_prop(N0=phi_tau_0_sw_spring*N0, mu=-log(phi), E=E_sw_spring, 
                                                               tau_mid=tau_sw_spring-tau_0_sw_spring, w=1, ice=1)
    H_sw_fall <- ((1-p_fi)*(1-p_sw_spring)*H) * hunting_bag_prop(N0=phi_tau_0_fi*(N0-H_fi-H_sw_spring), mu=-log(phi), E=E_sw_fall, 
                                                                 tau_mid=tau_sw_fall-tau_0_sw_fall, w=1, ice=1)
    # compute total expected hunting bag
    H_tot <- H_fi + H_sw_spring + H_sw_fall
    
    # compute probabilities of surviving natural mortality
    rho_S <- (1-H_tot/N0)*phi

    # compute left hand side of Euler-Lotka equation when lambda=1
    EL <- b/2*prod(rho_S[1:5])/(1-rho_S[6])
    
    return(EL)
  })
}

# objective function to solve for H such that R0 = 1
objective <- function(H, N, params) {
  with(params, {
    return((R0(H, N, params)-1)^2)
  })
}


#################### Analysis ###################################################

# posterior samples for haul-out vector during mean ice conditions (See section on haul-out model)
ice_mean <- mean(data$ice)
zeta_0_unscaled <- samples$d + (1-samples$d)*logistic(samples$alpha[,1]) # alpha is on mean-centered space, so alpha[2]*mean(ice) = 0
zeta_0 <- (1-samples$w)/samples$w*zeta_0_unscaled
zeta_0_pup <- (1-samples$w_pup)/samples$w_pup*zeta_0_unscaled
zeta_1 <- ice_mean^2 / (samples$f^2 + ice_mean^2)
w <- zeta_1 / (zeta_0 + zeta_1)
w_pup <- zeta_1 / (zeta_0_pup + zeta_1)
w_vec <- matrix(rep(w, 2*indices$a), nrow=nrow(samples$n0), ncol=2*indices$a)
w_vec[,1] <- w_vec[,indices$a+1] <- w_pup

# historical average for % of hunting that takes place in Finland
p_fi <- data$y_hb_fi / (data$y_hb_fi + data$y_hb_sw)
p_fi <- mean(p_fi[data$y_hb_fi > 0])

# historical average for % of Swedish hunting that takes place in spring
p_sw_spring <- data$y_hb_sw_spring / (data$y_hb_sw_spring + data$y_hb_sw_fall)
p_sw_spring <- mean(p_sw_spring[data$y_hb_sw > 0])

# average historical mean occurrence times of hunting
tau_fi <- mean(data$tau_fi[data$tau_fi>0])
tau_sw_spring <- mean(data$tau_sw_spring[data$tau_sw_spring>0])
tau_sw_fall <- mean(data$tau_sw_fall[data$tau_sw_fall>0])

# current estimated population size
N.current <- rowSums(samples$N[,,dim(samples$N)[3]]) # current estimated total population size
N.grid <- seq(1, 61e3, 6e3) # setup grid of population sizes

H.crit <- matrix(NA, nrow=nrow(samples$n0), ncol=length(N.grid)) # critical hunting level for a range of pop. sizes
H.crit.current <- rep(NA, nrow(samples$n0)) # critical hunting level for current estimated population size
H.target <- rep(NA, nrow(samples$n0)) # approx. hunting level needed for 7% growth based on current pop. size 

# iterate over each posterior sample
for(i in 1:nrow(samples$n0)) {
  
  # parameters needed for computing the critical hunting level
  params <- list(phi=samples$phi[i,], w_vec=w_vec[i,],
                 E_fi = samples$E_fi[i,], E_sw_spring = samples$E_sw_spring[i,], E_sw_fall = samples$E_sw_fall[i,],
                 b_max=samples$b_max[i], theta_0=samples$theta_0[i], theta_1=samples$theta_1[i], 
                 tau_0_fi=indices$tau_0_fi, tau_0_sw_spring=indices$tau_0_sw_spring, tau_0_sw_fall=indices$tau_0_sw_fall,
                 tau_fi=tau_fi, tau_sw_spring=tau_sw_spring, tau_sw_fall=tau_sw_fall,
                 ice=ice_mean, p_fi=p_fi, p_sw_spring=p_sw_spring)
  
  H.crit.current[i] <- optimize(objective, c(0, N.current[i]), tol=1e-5, N=N.current[i], params=params)$minimum
  
  # solve critical hunting level for each population size
  for(j in 1:length(N.grid)) {
    H.crit[i,j] <- optimize(objective, c(0, N.grid[j]), tol=1e-5, N=N.grid[j], params=params)$minimum
  }
  
  params$phi <- params$phi / 1.07 # convert R0 to lhs of Euler-Lotka equation with lambda=1.07
  H.target[i] <- optimize(objective, c(0, N.current[i]), tol=1e-5, N=N.current[i], params=params)$minimum
}

H.crit.quantiles <- t(apply(H.crit, 2, quantile, quantiles))

Q <- data$Q_fi[length(data$Q_fi)] + data$Q_sw[length(data$Q_sw)] # aggregate 2022 quotas
N_2022 <- quantile(rowSums(samples$N[,,length(data$years)]), quantiles) # total population size in 2022


#################### Output ###################################################

estimates <- read.csv(file = paste(output_directory, 'estimates.csv', sep=''), head=T)
estimates <- rbind(estimates, 
                 c('Critical hunting 2023', quantile(H.crit.current, quantiles)), 
                 c('Target hunting 2023', quantile(H.target, quantiles)))
write.csv(estimates, file = paste(output_directory, 'estimates.csv', sep=''))

label_k = function(x) {
  labels = paste(x, "e3", sep='')
  labels[x==0] <- 0
  return(labels)
}
ggplot() +
  geom_ribbon(aes(x=N.grid, ymin=H.crit.quantiles[,1], ymax=H.crit.quantiles[,3]), alpha=0.25) +
  geom_line(aes(x=N.grid, y=H.crit.quantiles[,2]), linewidth=1) +
  geom_point(aes(x=N_2022[2], y=Q), color='red', size=3) +
  geom_errorbarh(aes(xmin=N_2022[1], xmax=N_2022[3], y=Q), color='red', height=150, linewidth=1) +
  scale_x_continuous(limits=c(0,62500), expand=c(0,NA), breaks=seq(0,60e3,10e3), label = function(x){label_k(1e-3*x)}) +
  scale_y_continuous(limits=c(0,5100), breaks=seq(0,5.1e3,1e3), expand=c(0,0), label = function(x){label_k(1e-3*x)}) +
  labs(x = 'Number of ringed seals', y = 'Critical harvest level') +
  theme_classic() +
  theme(legend.position='none', text=element_text(size = font_size))

# save plot
ggsave(filename='critical_hunting_level.png', width=8, height=8, path=paste(output_directory, 'plots/', sep=''))


