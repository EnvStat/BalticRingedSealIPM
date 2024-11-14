
# density dependent birth rate
b.N <- function(b_max, theta_0, theta_1, N) {
    b_N = exp(log(b_max) - theta_0*(expm1(theta_1*sum(N))))
    return(b_N)
}

# compute expected composition of the hunting bag
hunting_bag <- function(N0, mu, E, Q, tau_mid, tau, w, ice) {
  n_eff = w*N0/ice
  H_tot = as.numeric(Q*(1-exp(-E%*%(n_eff*(1-exp(-mu*tau))/mu))))
  p_H = E*n_eff*exp(-mu*tau_mid)
  H = p_H/sum(p_H)*H_tot
  return(H)
}


# generate simulated data with true discrete distribution
simulate <- function(n0, years.sim, Q_fi.sim, Q_sw.sim, ice.sim, 
                              params, data, indices) {
  with(params, {
    
    #################### SETUP DATA STRUCTURES ################################################################
    
    # length of simulation period
    t <- length(years.sim)
    
    # assume future mean occurrence time of hunting is equal to historical average
    tau_fi <- mean(data$tau_fi[data$tau_fi>0])
    tau_sw_spring <- mean(data$tau_sw_spring[data$tau_sw_spring>0])
    tau_sw_fall <- mean(data$tau_sw_fall[data$tau_sw_fall>0])
    
    # standardize simulation years
    years_std <- (years.sim - mean(data$years)) / sd(data$years)
    
    #initialize matrices for storing simulation data
    N <- matrix(0, 12, t) # population state
    U <- matrix(NA, 12*5, t) # intermediate population state

    rho <- matrix(NA, 12, 5) # state transition probabilities (alternative seal fates)

    ######################### GLOBAL PARAMETETRS ###################################################################
    
    # aging matrix
    A = matrix(0, 12, 12)
    A[2:6, 1:5] <- A[8:12, 7:11] <- diag(1, 5)
    A[6, 6] <- A[12, 12] <- 1
    
    # birth rates
    X_b <- cbind(rep(1, length(years_std)), years_std) # covariate matrix for birth rate
    b0 <- (b_min + (b_max-b_min)*logistic(X_b%*%beta)) # birth rate at zero-density
    b <- rep(NA, length(b0)) # density dependent birth rate
    
    # haul-out probability
    X_w <- (ice.sim - mean(data$ice))/sd(data$ice)  # covariate matrix for haul-out
    X_w <- cbind(rep(1, length(years_std)), X_w)
    zeta_0_unscaled <- d + (1-d)*logistic(X_w%*%alpha)
    zeta_0 <- (1-w)/w*zeta_0_unscaled # rate of moving into water (adult)
    zeta_0_pup <- (1-w_pup)/w_pup*zeta_0_unscaled # rate of moving into water (pup)
    zeta_1 <- ice.sim^2 / (f^2 + ice.sim^2) # rate of moving onto ice

    w_t <- zeta_1 / (zeta_0 + zeta_1) # time-varying haul-out probability for 1+ year olds
    w_t_pup <- zeta_1 / (zeta_0_pup + zeta_1) # time-varying haul-out probability for pups
    
  
    ######################## FORWARD PROJECTION #######################################################################

    N[,1] <- n0 # initial population state

    for (i in 1:t) {
      
      # generate births
      if(i > 1) {
        # density dependent birth rate
        b[i] <- b.N(b0[i], theta_0, theta_1,  N[,i-1])

        N[,i] <- A%*%U[1:12, i-1] # aging

        newborns <- rbinom(1, N[6,i], b[i]) # total number of newborn pups
        newborn.females <- rbinom(1, newborns, 0.5) # sex allocation of newborn pups
        N[1,i] = N[1,i] + newborn.females
        N[7,i] = N[7,i] + newborns - newborn.females
      }
      if(sum(N[,i]) <= 0) {break} # stop simulation if population has gone extinct
      
      # probabilities of surviving to beginning of hunting season
      phi_tau_0_fi <- phi^indices$tau_0_fi # FI
      phi_tau_0_sw_spring <- phi^indices$tau_0_sw_spring # SW spring
      phi_tau_0_sw_fall <- phi^indices$tau_0_sw_fall # SW fall
      
      # vector of class-dependent haul-out probabilities
      w_vec <- rep(w_t[i], 12)
      w_vec[c(1,7)] <- w_t_pup[i]

      # stochastic hunting effort
      E_fi <- (E_tot_fi * psi_fi)*exp(sigma_E_fi*rnorm(1))
      E_sw_spring <- (E_tot_sw_spring * psi_sw_spring)*exp(sigma_E_sw_spring*rnorm(1))
      E_sw_fall <- (E_tot_sw_fall * psi_sw_fall)*exp(sigma_E_sw_fall*rnorm(1))
      
      # numerical solutions to hunting
      n_hunted_fi <- n_hunted_sw_spring <- n_hunted_sw_fall <- rep(0, length(N[,i]))
      if(Q_fi.sim[i] > 0) {
        n_hunted_fi <- hunting_bag(N0 = phi_tau_0_fi*N[,i], 
                                          mu = -log(phi), E = E_fi, Q = Q_fi.sim[i], 
                                          tau_mid = tau_fi - indices$tau_0_fi, tau = indices$tau_basking-indices$tau_0_fi,
                                          w=w_vec, ice=ice.sim[i])
      }
      if(Q_sw.sim[i] > 0) {
        n_hunted_sw_spring <- hunting_bag(N0 = phi_tau_0_sw_spring*N[,i], 
                                        mu = -log(phi), E = E_sw_spring, Q = Q_sw.sim[i], 
                                        tau_mid = tau_sw_spring - indices$tau_0_sw_spring, tau = 2/12,
                                        w=1, ice=1)
        n_hunted_sw_fall <- hunting_bag(N0 = phi_tau_0_sw_fall*(N[,i]-n_hunted_sw_spring-n_hunted_fi), 
                                        mu = -log(phi), E = E_sw_fall, Q = Q_sw.sim[i]-sum(n_hunted_sw_spring), 
                                        tau_mid = tau_sw_fall - indices$tau_0_sw_fall, tau = 2/12,
                                        w=1, ice=1)
      }
      
      # log probabilities of alternative seal fates
      rho[,3] <- n_hunted_sw_spring / N[,i] # SW spring hunting
      rho[,4] <- n_hunted_sw_fall / N[,i] # SW fall hunting
      rho[,5] <- n_hunted_fi / N[,i] # FI hunting
      rho[,1] <- (1 - rowSums(rho[,3:5])) * phi
      rho[,2] <- 1 - rowSums(rho[,c(1,3:5)])
      
      if(min(rho) < 0) {rho[rho < 0]=0} # in case of numerical instability
      
      # multi-state survival process
      U_stacked <- matrix(0, nrow = 12, ncol = 5) # vector U arranged in matrix form
       for (j in 1:nrow(U_stacked)) {
         if(N[j,i] > 0) {
           U_stacked[j,] <- t(rmultinom(1, N[j,i], rho[j,]))
         }
       }
      U[,i] <- c(U_stacked)

    }

    return(N)
  })
}


