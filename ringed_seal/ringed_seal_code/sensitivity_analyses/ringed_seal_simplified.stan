
functions {
  // function for computing mortality rates for subadults and males (see Section "Survival and mortality")
  vector mortality_rates(real mu_f0, real mu_f5, real c, vector nu) {
    vector[12] mu;
    mu[1] = mu_f0;
    mu[2:6] = exp( log(mu_f0) - [1^c,2^c,3^c,4^c,5^c]'/5^c*(log(mu_f0)-log(mu_f5)) );
    mu[7] = mu[1]*exp(nu[1]);
    mu[8:12] = mu[2:6]*exp(nu[2]); 
    return mu;
  }
  
  // function for computing the expected composition of the hunting bag (see Appendix S3)
  vector n_hunted(vector N0, vector mu, vector E, real Q, real tau_mean, real tau) {
    int c = num_elements(N0);
    vector[c] H = rep_vector(0, c);
    if(Q > 0) {
      vector[c] N_eff = N0 .* (1-exp(-mu*tau)) ./ mu;
      real log_H_tot = log(Q) + log1m_exp(-E' * N_eff);
      H = exp(log(E) + log(N0 .* exp(-mu*tau_mean)) - log(E' * (N0 .* exp(-mu*tau_mean))) + log_H_tot);
    }
    return H;
  }
  
  // function for computing density dependent birth rate (see Appendix S4)
  real DD(real b0, real N, real theta_0, real theta_1) {
    real b_t = exp( log(b0) - exp( log(theta_0) + log_diff_exp(theta_1*N, 0) ) );
    return b_t;
  }
  
  // function for computing birth rate at carrying capacity (from Euler-Lotka Equation)
  real b_equilibrium(vector phi) {
    real l = exp(sum(log(phi[1:5])));
    real b_eq = 2*(1-phi[6])/l;
    return b_eq;
  }
  
  // function for computing equilibrium demographic structure (analytical calculation of dominant eigenvector, see Caswell (2001))
  vector equilibrium_structure(vector phi_eq) {
    vector[6] l_f = append_row(1, exp(cumulative_sum(log(phi_eq[1:5]))));
    vector[6] l_m = append_row(1, exp(cumulative_sum(log(phi_eq[7:11]))));
    vector[12] l = append_row(l_f, l_m);
    vector[6] c_f = append_row(rep_vector(1, 5), 1/(1-phi_eq[6]));
    vector[6] c_m = append_row(rep_vector(1, 5), 1/(1-phi_eq[12]));
    vector[12] c = append_row(c_f, c_m);
    return (c .* l) / (c'*l);
  }
}

data {
  
  ///////////////////// Miscellaneous ///////////////////////////////////////
  
  // indices
  int<lower=1> a; // number of age classes
  int<lower=1> k; // number of mortality sources
  int<lower=1> t; // length of study period
  int<lower=1> t_h_sw; int<lower=1> t_h_fi; // number of years with legal hunting in Sweden and Finland
  int years[t]; // years included in the study

  int<lower=0, upper=1> exclude_survey[t+1]; // whether or not to exclude survey data in each year (for sensitivity)
  int<lower=0, upper=1> exclude_rep[t]; // whether or not to exclude reproduction data in each year (for sensitivity)
  int<lower=0, upper=1> exclude_hb[t]; // whether or not to exclude hunting bag data in each year (for sensitivity)
  int<lower=0, upper=1> exclude_samples[t]; // whether or not to exclude hunting/bycatch sample data in each year (for sensitivity)
  
  // prior info
  vector[2*a] N0; // initial guess for demographic structure (used as initialization for eigenvector computation)
  matrix[2*a-1, 2*a-1] L_bias; // cholesky decompositions of prior covariance matrix for hunting/bycatch bias
  vector<lower=0, upper=1>[2] w_bounds; // lower and upper bounds on w prior (for sensitivity analysis)
  
  // mean occurrence times of hunting
  vector<lower=0>[t] tau_fi; 
  vector<lower=0>[t] tau_sw_fall; 
  vector<lower=0>[t] tau_sw_spring;

  // covariates
  vector[t] years_std; // years standardized to have mean=0 and sd=1
  vector[t+1] ice; // April sea ice extent in Bothnian Bay (in 000's of km2)
  vector[t] Q_fi; vector[t] Q_sw; // hunting quotas
  
  // time from parturition to the onset of hunting season (in years)
  real<lower=0> tau_0_fi; real<lower=0> tau_0_sw_spring; real<lower=0> tau_0_sw_fall;
  // avg. time between mating and sampling of pregnant females (in years)
  real<lower=0> tau_p;
  // length of ecological seasons (in years)
  real<lower=0> tau_basking; real<lower=0> tau_foraging; real<lower=0> tau_subnivean;

  /////////////////////// Observations ////////////////////////////////////////////
  
  // aerial surveys
  int<lower=0> y_survey[t+1];
  
  // hunting bags
  real<lower=0> y_hb_fi[t]; // size of Finnish hunting bag
  real<lower=0> y_hb_sw_spring[t]; // size of Swedish hunting bag in the spring
  real<lower=0> y_hb_sw_fall[t]; // size of Swedish hunting bag in the fall
  real<lower=0> y_hb_sw[t]; // total size of Swedish hunting bag
  
  // hunting samples
  int<lower=0> y_hs_fi[2*a, t];
  int<lower=0> y_hs_sw_spring[2*a, t];
  int<lower=0> y_hs_sw_fall[2*a, t];

  // hunting samples with missing age or sex info.
  int<lower=0> y_hs_fi_sexNA[a, t];
  int<lower=0> y_hs_fi_ageNA[2, t];
  int<lower=0> y_hs_sw_ageNA_spring[2, t];
  int<lower=0> y_hs_sw_sexNA_spring[a, t];
  int<lower=0> y_hs_sw_ageNA_fall[2, t];
  int<lower=0> y_hs_sw_sexNA_fall[a, t];
  
  // sampling bias data for spring Swedish hunting
  int<lower=0> o_16_19[2, 2]; // number sampled by size group (2016-2019)
  int<lower=0> o_20_21[2, 2]; // number sampled by size group (2020-2021)
  int<lower=0> l[2*a, 2]; // > 100cm by class
  
  // bycaught samples
  int<lower=0> y_bc[2*a, t];
  
  // bycaught samples with missing age or sex info.
  int<lower=0> y_bc_sexNA[a, t]; 
  int<lower=0> y_bc_ageNA[2, t];

  // pregnancy data  
  int<lower=0> pregnant_ss[t]; // sample sizes for embryos
  int<lower=0> y_pregnant[t]; // observed pregnancies
  
  // outcomes of reproductive assessments (joint outcome of CA and placental scars)
  int<lower=0> z[6, t];
    
  // corpus albicans data
  int<lower=0> CA_ss[t]; // sample sizes for CA only
  int<lower=0> y_CA[t]; // observed CA (placental scar not evaluated)
  
  // placental scar data
  int<lower=0> scar_ss[t]; // sample sizes for placental scars only
  int<lower=0> y_scar[t]; // observed placental scars (CA not evaluated)
}

transformed data {

  // model dimension changes with the re-introduction of Swedish and Finnish hunting in 2015 and 2016
  int<lower=1> t0 = t - t_h_sw; // last year without hunting (2014)
  int<lower=1> t1 = t - t_h_fi; // last year with only Swedish hunting (2015)
  
  // aging matrix (moves all seals up one age class)
  matrix[2*a, (1+k)*(2*a)] A = rep_matrix(0, 2*a, (1+k)*(2*a));

  vector[t+1] ice_std; // April sea ice extent in Bothnian Bay transformed to mean=0 and sd=1
  matrix[t+1, 2] X_b; // design matrix for birth rate linear model
  matrix[t+1, 2] X_w; // design matrix for haul-out rate linear model

  // standardized ice cover
  ice_std = (ice - mean(ice)) / sd(ice);

  // design matrix for birth rate linear model
  X_b[,1] = rep_vector(1, t+1); // intercept term
  X_b[,2] = append_row(years_std, years_std[t]+(years_std[2]-years_std[1])); // slope term (2023 is appended)
  
  //design matrix for haul-out linear model
  X_w[,1] = rep_vector(1, t+1); // intercept term
  X_w[,2] = ice_std; // slope term
  
  // construct aging matrix
  A[2:a, 1:(a-1)] = diag_matrix(rep_vector(1, a-1)); A[a, a] = 1; // female aging
  A[(a+2):(2*a), (a+1):(2*a-1)] = diag_matrix(rep_vector(1, a-1)); A[2*a, 2*a] = 1; // male aging

}

parameters {
  
  real<lower=0> n0_unscaled; // initial population size (in 000's)
  
  /////////////////// Reproduction //////////////////////////////
  
  real<lower=0, upper=1> b_max; // asymptotic max. birth rate
  real<lower=0, upper=1> b_scale; // ratio of minimum to maximum birth rate
  vector[2] beta; // slope and intercept parameters for birth rate logistic model
  
  /////////////////// Natural mortality //////////////////////////////
  
  real<lower=0, upper=1> phi_f5_unscaled; // untransformed survival probability for female adults
  real<lower=0, upper=1> phi_scale; // ratio of female pup survival to female adult survival probability
  vector[2] nu_unscaled; // deviations of male survival probabilities from females in logit space (x10 scale)
  real<lower=0, upper=1> c; // shape parameter for survival probability of sub-adults
  
  ////////////////// Anthropogenic mortality //////////////////////////
  
  // reparametrized median hunting efforts (in terms of expected % of quota filled)
  real<lower=0, upper=1> pq_fi; 
  real<lower=0, upper=1> pq_sw_spring; 
  real<lower=0, upper=1> pq_sw_fall;
  
  // noise terms for hunting effort in each year
  vector[t_h_fi] eps_E_fi; 
  vector[t_h_sw] eps_E_sw_spring; 
  vector[t_h_sw] eps_E_sw_fall;
  
  // standard deviations for hunting effort in log-space
  real<lower=0> sigma_E_fi; 
  real<lower=0> sigma_E_sw_spring; 
  real<lower=0> sigma_E_sw_fall;
  
  // untransformed hunting/bycatch bias vectors in logit-space
  vector[2*a-1] g_fi;
  vector[2*a-1] g_sw_spring;
  vector[2*a-1] g_sw_fall;
  vector[2*a-1] g_bycatch;

  ///////////////////// Demographic stochasticity /////////////////////////////
  
  vector[t] u_pup; // stochastic noise for pup births
  vector[t] u_sex; // stochastic noise for sex allocation of pups
  // model dimension changes twice - years divided into pre-hunting, only Swedish hunting and Swedish+Finnish hunting periods
  matrix[(k-3)*(2*a), t0] u0; // stochastic noise for natural mortality only (pre-hunting period)
  matrix[(k-1)*(2*a), t_h_sw-t_h_fi] u1; // stochastic noise for multinomial state transition (period with only Swedish hunting)
  matrix[k*(2*a), t_h_fi] u2; // stochastic noise for multinomial state transition (period with Swedish & Finnish hunting)

  /////////////////////// Haul-out ////////////////////////////////////////////
  
  real<lower=0, upper=1> w_unscaled; // untransformed baseline haul-out probability for sub-adults and adults
  real<lower=0, upper=1> delta; // ratio of pup haul-out % to sub-adult and adult %
  real<lower=-pi()/2, upper=pi()/2> alpha_1_raw; // untransformed slope term for logistic haulout model
  real<lower=-pi()/2, upper=pi()/2> alpha_0_scale_raw; // untransformed ratio of intercept to slope in logistic haulout model
  //real<lower=0> f; // half-saturation constant for movement rate from water to ice
  //real<lower=0, upper=1> d; // ratio of lwr to upr asymptote for movement rate from ice to water
  
  ////////////////////////// Observation process ///////////////////////////////
  
  real<lower=0> sqrt_r_inv; // inverse sqrt of overdispersion parameter for aerial survey estimates
  
  vector<lower=0, upper=1>[2] x_16_19; // sampling probabilities for <100cm and >100cm (2016-2019)
  vector<lower=0, upper=1>[2] x_20_21; // sampling probabilities for <100cm and >100cm (2020-2021)
  vector<lower=0, upper=1>[2*a] q; // probability of length > 100cm by class
  
  real<lower=0, upper=1> pi_s; // detection prob. for placental scars
  real<lower=0, upper=1> pi_c; // detection prob. for CA
  real<lower=0, upper=1> kappa; // probability that seals have CA without having been pregnant
  
  real<lower=0, upper=1> scar_na_0; // prob. scar=NA given CA=0
  real<lower=0, upper=1> scar_na_1; // prob. scar=NA given CA=1

  //////////////////////////// Density dependence /////////////////////////////////
  
  real<lower=0> K_unscaled; // carrying capacity (in scale of 000's)
  real<lower=0, upper=1> theta_0_raw; // untransformed "background" failure rate for pregnancies

}

transformed parameters {
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// VARIABLE DECLARATIONS ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////// Scaled & transformed parameters ///////////////////////////
  
  // initial population size
  real<lower=0> n0;
  
  // natural mortality
  real<lower=0, upper=1> phi_f5; // survival probability of female adults
  real<lower=0, upper=1> phi_f0; // survival probability of female pups
  real<lower=0> mu_f0; // mortality rate of female pups
  real<lower=0> mu_f5; // mortality rate of female adults
  vector[2] nu; // deviation of male mortality rates from females in logit space
  
  // reproduction
  real<lower=0, upper=1> b_min; // minimum birth rate

  // anthropogenic mortality
  real<lower=0> E_tot_fi; // magnitude of Finnish hunting effort
  real<lower=0> E_tot_sw_spring; // magnitude of Swedish hunting effort in spring
  real<lower=0> E_tot_sw_fall; // magnitude of Swedish hunting effort in fall

  // haul-out probability
  real<lower=0, upper=1> w; // baseline haul-out probability of sub-adults and adults
  real<lower=0, upper=1> w_pup; // baseline haul-out probability of pups
  vector[2] alpha; // slope and intercept for haul-out model
  
  // density dependence
  real<lower=0> K; // carrying capacity
  real<lower=0> theta_0; // background failure rate for pregnancies
  
  // observations
  real<lower=0> r; // overdispersion parameter for aerial surveys

  /////////////// State variables /////////////////////////////////////////
  
  matrix[2*a, t+1] N; // population matrix (number of classes)x(number of years)
  real<lower=0> n_pup; // number of pups born
  real<lower=0> n_pup_f; real<lower=0> n_pup_m; // number of female and male pups born
  
  matrix[(1+k)*(2*a), t] U; // intermediate population matrix ([number of states]x[number of classes]) x (number of years)
  matrix[2*a, k+1] U_stacked; // rearranged form of U[,t]. Rows = classes & columns = transition states
  matrix[2*a, k-1] u1_stacked; matrix[2*a, k] u2_stacked; // similar rearranged forms of u1[,t] and u2[,t]

  vector[t+1] N_w; // total hauled-out population size
  matrix[2*a, k] N_dead[t]; // dead seals (number of classes)x(number of mortality sources)
  
  //////////////// Mortality ///////////////////////////////////////////////
  
  // natural mortality
  vector<lower=0>[2*a] mu; // mortality rates
  vector<lower=0, upper=1>[2*a] phi; // survival probabilities

  // bias vectors for hunting & bycatch
  simplex[2*a] psi_fi; 
  simplex[2*a] psi_sw_spring; 
  simplex[2*a] psi_sw_fall;
  simplex[2*a] psi_bycatch;

  // expected composition of the hunting bags
  vector<lower=0>[2*a] n_hunted_fi; 
  vector<lower=0>[2*a] n_hunted_sw_spring; 
  vector<lower=0>[2*a] n_hunted_sw_fall;

  // matrices for state transitions
  matrix[2*a, 2*a] S; // survival probabilities
  matrix[2*a, 2*a] D; // natural mortality probabilities
  matrix[2*a, 2*a] H_sw_spring; // probabilities of being hunted in Sweden in spring
  matrix[2*a, 2*a] H_sw_fall; // probabilities of being hunted in Sweden in fall
  matrix[2*a, 2*a] H_fi; // probabilities of being hunted in Finland
  matrix[(1+k)*(2*a), (2*a)] M; // composite state transition matrix [S', D', H_sw_spring', H_sw_fall', H_fi']'
  matrix[2*a, 1+k] rho[t]; // state transition probabilities (number of classes)x(number of states)

  // median age-and-sex dependent hunting efforts
  vector<lower=0>[2*a] E_fi; 
  vector<lower=0>[2*a] E_sw_spring; 
  vector<lower=0>[2*a] E_sw_fall;
  
  // stochastic age-and-sex dependent hunting efforts
  vector<lower=0>[2*a] E_fi_t; 
  vector<lower=0>[2*a] E_sw_t_spring; 
  vector<lower=0>[2*a] E_sw_t_fall;
  
  // probabilities of surviving until the beginning of the hunting season
  vector<lower=0, upper=1>[2*a] phi_tau_0_fi;
  vector<lower=0, upper=1>[2*a] phi_tau_0_sw_spring; 
  vector<lower=0, upper=1>[2*a] phi_tau_0_sw_fall;

  //---------- logit-normal approximation to multinomial state transitions ---------------------
  
  row_vector<lower=0, upper=1>[1+k] allocation; // expected outcome of state transitions
  matrix<lower=0>[2*a, 1+k] eta_p; row_vector<lower=0>[1+k] eta_p_adj;
  
  // period with no hunting
  real mean_0; // expectation of g
  real<lower=0> sigma2_0; // variance of g

  // period with only Swedish hunting
  row_vector[k-1] mean_1;
  matrix[k-1, k-1] Sigma_1; // covariance matrix for logit-normal distribution 
  matrix[k-1, k-1] L_1; // cholesky decomposition Sigma_1
  
  // period with Swedish & Finnish hunting
  row_vector[k] mean_2;
  matrix[k, k] Sigma_2;
  matrix[k, k] L_2; // cholesky decomposition of Sigma_2
  
  ///////////////// Reproduction ///////////////////////////////////////////
  
  // pregnancy rates
  vector<lower=0, upper=1>[t] p; // fall pregnancy rate over time
  real<lower=0, upper=1> p_max;
  real<lower=0, upper=1> p_min;
  
  // birth rate
  vector<lower=0, upper=1>[t+1] b0; // birth rate at zero population density
  vector<lower=0, upper=1>[t+1] b; // density dependent birth rate
  matrix[2*a, 2*a] B = diag_matrix(rep_vector(1, 2*a)); // matrix projection for births
  
  // probabilities of reproductive assessment outcomes
  matrix<lower=0, upper=1>[6, t] gamma;
  
  /////////////////// Haul-out ///////////////////////////////////////////////
  
  // movement rates
  vector<lower=0>[t+1] zeta_0; // rate of moving into water (sub-adults & adults)
  vector<lower=0>[t+1] zeta_0_pup; // rate of moving into water (pups)
  vector<lower=0, upper=1>[t+1] zeta_1; // rate of moving onto ice
  real<lower=0> f; // half-saturation constant for movement rate from water to ice
  real<lower=0, upper=1> d; // ratio of lwr to upr asymptote for movement rate from ice to water
  
  // haul-out probabilities
  vector<lower=0, upper=1>[t+1] w_t; // yearly expected haul-out probabilities (sub-adults & adults)
  vector<lower=0, upper=1>[t+1] w_t_pup; // yearly expected haul-out probabilities (pups)
  vector[2*a] w_vec; // vector of age-dependent haul-out probabilities
  
  ///////////////////// Density Dependence ////////////////////////////////////

  real<lower=0> theta_1; // strength of density dependence
  
  ///////////////////// Observation processes /////////////////////////////////
  
  vector<lower=0, upper=1>[2*a] v_16_19; // sampling probability by class (2016-2019)
  vector<lower=0, upper=1>[a] v_20_21; // sampling probability by age for females (2020-2021)
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// COMPUTE GLOABAL PARAMETERS /////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // initial population size
  n0 = 1e3*n0_unscaled;
  
  // natural mortality
  phi_f5 = 0.80 + (0.97-0.80)*phi_f5_unscaled;
  phi_f0 = phi_scale*phi_f5;
  mu_f0 = -log(phi_f0);
  mu_f5 = -log(phi_f5);
  nu = 0.1*nu_unscaled;
  mu = mortality_rates(mu_f0, mu_f5, c, nu); // compute mortality rates for all age and sex classes
  phi = exp(-mu);
  phi_tau_0_fi = exp(-mu*tau_0_fi); 
  phi_tau_0_sw_spring = exp(-mu*tau_0_sw_spring); 
  phi_tau_0_sw_fall = exp(-mu*tau_0_sw_fall); 
  
  // anthropogenic mortality
  psi_sw_spring = softmax(append_row(L_bias*g_sw_spring, 0)); 
  psi_sw_fall = softmax(append_row(L_bias*g_sw_fall, 0)); 
  psi_fi = softmax(append_row(L_bias*g_fi, 0)); 
  psi_bycatch = softmax(append_row(L_bias*g_bycatch, 0)); 
  E_tot_fi = -log(1-pq_fi)/5;
  E_tot_sw_spring = -log(1-pq_sw_spring)/100;
  E_tot_sw_fall = -log(1-pq_sw_fall)/100;             
  E_fi = E_tot_fi * psi_fi; 
  E_sw_spring = E_tot_sw_spring * psi_sw_spring; 
  E_sw_fall = E_tot_sw_fall * psi_sw_fall;
  
  // density dependence
  K = 1e5*K_unscaled;
  theta_0 = -log(b_max + (1-b_max)*theta_0_raw);
  theta_1 = 1/K*log(1-(log(b_equilibrium(phi))-log(b_max))/theta_0); // see Appendix S2
  
  // birth
  b_min = b_scale*b_max;
  b0 = (b_min + (b_max-b_min)*inv_logit(X_b*beta));
  p_max = b_max*exp((1-tau_p*(2-tau_p))*theta_0);
  
  // haul-out probability
  d = 0; f = 0;
  w = w_bounds[1]+(w_bounds[2]-w_bounds[1])*w_unscaled;
  w_pup = delta*w;
  alpha = 1.25*[5*tan(alpha_0_scale_raw)*tan(alpha_1_raw), tan(alpha_1_raw)]'; // inverse CDF method to transform uniform prior to Cauchy
  zeta_0 = (1-w)/w*(d + (1-d)*inv_logit(X_w*alpha)); // ice to water (adult)
  zeta_0_pup = (1-w_pup)/w_pup*(d + (1-d)*inv_logit(X_w*alpha));  // ice to water (pup)
  zeta_1 = rep_vector(1, t+1); //(ice .* ice) ./ (f^2 + (ice .* ice)); // water to ice
  w_t = zeta_1 ./ (zeta_0 + zeta_1); // time-varying haul-out probability for sub-adults and adults
  w_t_pup = zeta_1 ./ (zeta_0_pup + zeta_1); // time-varying haul-out probability for pups

  // observation parameters
  r = inv(sqrt_r_inv)^2;
  v_16_19 = (1-q)*x_16_19[1] + q*x_16_19[2];
  v_20_21 = (1-q[1:6])*x_20_21[1] + q[1:6]*x_20_21[2];
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// PROCESS MODEL //////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(i in 1:t) {

    //////////////// Initialization //////////////////////////////////
    if(i == 1) {
      b[i] = DD(b0[i], n0, theta_0, theta_1); // compute density dependent birth rate
      B[1, a] = b[i] / 2; B[a+1, a] = b[i] / 2; // set birth rate of female and male pups
      S = diag_matrix(phi); // construct diagonal survival matrix
      N[,i] = N0; // guess for initial age structure
      
      // iterate matrix projection until convergence to stable stage structure (power method for computing dominant eigenvector)
      for(j in 1:20) { N[,i] = B*A[1:(2*a), 1:(2*a)]*S[1:(2*a), 1:(2*a)]*N[,i]; }
      N[,i] = (N[,i] / sum(N[,i])) * n0; // scale to initial population size
    }

    //////////////// Aging and birth //////////////////////////////////
    if(i > 1) {
      b[i] = DD(b0[i], sum(N[,i-1]), theta_0, theta_1); // density dependent birth rate
      p[i-1] = exp(log(b0[i])+theta_0*(1-tau_p*(2-tau_p)*exp(theta_1*sum(N[,i-1])))); // pregnancy rate during the previous fall

      N[,i] = A*U[,i-1]; // move seals to the next age class

      n_pup = b[i] * N[a,i] + sqrt(N[a,i] * b[i] * (1-b[i])) * u_pup[i-1]; // total number of pups born
      n_pup_f = 0.5*n_pup + sqrt(n_pup * 0.5*(1-0.5)) * u_sex[i-1]; // number of female pups
      n_pup_m = n_pup - n_pup_f; // number of male pups
      N[1,i] = N[1,i] + n_pup_f;
      N[a+1,i] = N[a+1,i] + n_pup_m;
    }


    /////////////////// Hunting (see Appendix S3) ///////////////////////////////////////
    
    // vector of haul-out probabilities
    w_vec = rep_vector(w_t[i], 2*a);
    w_vec[1] = w_t_pup[i]; 
    w_vec[a+1] = w_t_pup[i];

    // initialize empty hunting bags
    n_hunted_sw_spring = rep_vector(0, 2*a);
    n_hunted_sw_fall = rep_vector(0, 2*a);
    n_hunted_fi = rep_vector(0, 2*a);

    // Finnish hunting
    if(i > t1) {
      // stochastic hunting effort
      E_fi_t = E_fi * exp(sigma_E_fi*eps_E_fi[i-t1]);
      
      // expected composition of the hunting bag
      n_hunted_fi = n_hunted((phi_tau_0_fi .* w_vec .* N[,i])/ice[i], mu, E_fi_t, Q_fi[i], 
                              tau_fi[i]-tau_0_fi, tau_basking-tau_0_fi);
    }
    // Swedish hunting
    if(i > t0) {
      // stochastic hunting effort
      E_sw_t_spring = E_sw_spring * exp(sigma_E_sw_spring*eps_E_sw_spring[i-t0]);
      E_sw_t_fall = E_sw_fall * exp(sigma_E_sw_fall*eps_E_sw_fall[i-t0]);
      
      // expected composition of the hunting bags
      n_hunted_sw_spring = n_hunted(phi_tau_0_sw_spring .* N[,i], 
                                    mu, E_sw_t_spring, Q_sw[i], tau_sw_spring[i]-tau_0_sw_spring, 2.0/12);
      n_hunted_sw_fall = n_hunted(phi_tau_0_sw_fall .* (N[,i]-n_hunted_sw_spring-n_hunted_fi), 
                                  mu, E_sw_t_fall, Q_sw[i]-sum(n_hunted_sw_spring), tau_sw_fall[i]-tau_0_sw_fall, 2.0/12);
    }

    //////////////////// Mortality ////////////////////////////////////////
    
    // state transition probabilities
    rho[i,,3] = n_hunted_sw_spring ./ N[,i];
    rho[i,,4] = n_hunted_sw_fall ./ N[,i];
    rho[i,,5] = n_hunted_fi ./ N[,i];
    rho[i,,1] = (1-rho[i,,3]-rho[i,,4]-rho[i,,5]) .* phi;
    rho[i,,2] = 1 - rho[i,,1] - rho[i,,3] - rho[i,,4] - rho[i,,5];

    // create composite state transition matrix
    S = diag_matrix(rho[i,,1]);
    D = diag_matrix(rho[i,,2]);
    H_sw_spring = diag_matrix(rho[i,,3]);
    H_sw_fall = diag_matrix(rho[i,,4]);
    H_fi = diag_matrix(rho[i,,5]);
    M = append_row(S, append_row(D, append_row(H_sw_spring, append_row(H_sw_fall, H_fi))));

    eta_p = to_matrix(M*N[,i], 2*a, 1+k); // expected outcome of state transition


    // logit-Gaussian approximation to multinomial state transitions (see Appendix S7)

    // period without legal hunting
    if(i <= t0) {
      for(j in 1:(2*a)) {
        mean_0 = digamma(eta_p[j,2]) - digamma(eta_p[j,1]);
        sigma2_0 = trigamma(eta_p[j,2]) + trigamma(eta_p[j,1]);
        
        // survivors defined as the reference level
        allocation = append_col(softmax([0, mean_0 + sqrt(sigma2_0)*u0[j,i]]')', [0, 0, 0]);
        U_stacked[j,] = allocation*N[j,i]; // stochastic outcome of state transition
      }
    }
    // period with only Swedish hunting
    if(i > t0 && i <= t1) {
        u1_stacked = to_matrix(u1[,i-t0], 2*a, k-1);
        for(j in 1:(2*a)) {
          eta_p_adj = eta_p[j,]*(1+1/min(eta_p[j,1:4]));
          mean_1 = digamma(eta_p_adj[2:4]) - digamma(eta_p_adj[1]);
          Sigma_1 = rep_matrix(trigamma(eta_p_adj[1]), k-1, k-1) + diag_matrix(trigamma(eta_p_adj[2:4])');
          L_1 = cholesky_decompose(Sigma_1);
          
          // survivors defined as the reference level
          allocation = append_col(softmax(append_col(0, mean_1 + u1_stacked[j,]*L_1')')', 0);
          U_stacked[j,] = allocation*N[j,i]; // stochastic outcome of state transition
        }
    }
    // period with Finnish & Swedish hunting
    if(i > t1) {
        u2_stacked = to_matrix(u2[,i-t1], 2*a, k);
        for(j in 1:(2*a)) {
          eta_p_adj = eta_p[j,]*(1+1/min(eta_p[j,1:5]));
          mean_2 = digamma(eta_p_adj[2:5]) - digamma(eta_p_adj[1]);
          Sigma_2 = rep_matrix(trigamma(eta_p_adj[1]), k, k) + diag_matrix(trigamma(eta_p_adj[2:5])');
          L_2 = cholesky_decompose(Sigma_2);
          
          // survivors defined as the reference level
          allocation = softmax(append_col(0, mean_2 + u2_stacked[j,]*L_2')')';
          U_stacked[j,] = allocation*N[j,i]; // stochastic outcome of state transition
        }
    }

    N_dead[i] = U_stacked[,2:(1+k)];
    U[,i] = to_vector(U_stacked);

    // compute probabilities of each reproductive assessment outcome
    gamma[2,i] = (b[i]*(1-pi_s)*pi_c + (1-b[i])*kappa*pi_c)*(1-scar_na_1); // CA=1 & scar=0
    gamma[3,i] = b[i]*pi_s*(1-pi_c)*(1-scar_na_0); // CA=0 & scar=1
    gamma[4,i] = b[i]*pi_s*pi_c*(1-scar_na_1); // CA=1 & scar=1
    gamma[5,i] = (b[i]*pi_c + (1-b[i])*kappa*pi_c)*scar_na_1; // CA=1 & scar=NA
    gamma[6,i] = (b[i]*(1-pi_c) + (1-b[i])*(kappa*(1-pi_c) + 1-kappa))*scar_na_0; // CA=0 and scar=NA
    gamma[1,i] = 1 - sum(gamma[2:6,i]); // CA=0 & scar=0

    // total hauled-out population size
    // aerial survey assumed to coincide with onset of FI hunting season ~April 16
    N_w[i] = (phi_tau_0_fi .* w_vec)' * N[,i];
  }

  // final population state at the beginning of 2023
  b[t+1] = DD(b0[t+1], sum(N[,t]), theta_0, theta_1);
  N[,t+1] = A*U[,t];
  n_pup = b[t+1] * N[a,t+1] + sqrt(N[a,t+1] * b[t+1] * (1-b[t+1])) * u_pup[t];
  n_pup_f = 0.5*n_pup + sqrt(n_pup * 0.25) * u_sex[t];
  n_pup_m = n_pup - n_pup_f;
  N[1,t+1] = N[1,t+1] + n_pup_f;
  N[a+1,t+1] = N[a+1,t+1] + n_pup_m;
  
  // final hauled-out population size
  w_vec = rep_vector(w_t[t+1], 2*a);
  w_vec[1] = w_t_pup[t+1]; 
  w_vec[a+1] = w_t_pup[t+1];
  N_w[t+1] = (phi_tau_0_fi .* w_vec)' * N[,t+1];
  
  // latest estimate of fall pregnancy rate (2022)
  p[t] = exp( log(b0[t+1]) + theta_0*(1 - tau_p*(2-tau_p)*exp(theta_1*sum(N[,t]))) );
  
  // minimum historical pregnancy rate assuming population size was the same as 1988
  p_min = exp( log(b_min) + theta_0*(1 - tau_p*(2-tau_p)*exp(theta_1*sum(N[,1]))) );

}

model {
  // Swedish spring hunting weighted by sampling bias
  vector[2*a] wgt_samples_sw_16_19;
  vector[a] wgt_samples_sw_20_21;

  // natural mortality weighted by bycatch deviation
  vector[2*a] wgt_samples_bycatch;

  ///////////////////// PRIORS (see Appendix S2) /////////////////////////////////////////////////////

  // initial population size
  n0_unscaled ~ lognormal(log(2.5/0.65), 1);

  // Reproduction
  b_max ~ uniform(0, 1);
  b_scale ~ uniform(0, 1);
  beta[1] ~ cauchy(0, 10);
  beta[2] ~ student_t(5, 0, 1.8);

  // Natural mortality
  phi_f5_unscaled ~ uniform(0, 1);
  phi_scale ~ uniform(0, 1);
  nu_unscaled ~ normal(0, 1);
  c ~ uniform(0,1);

  // median hunting effort (reparametrized)
  pq_fi ~ uniform(0,1);
  pq_sw_spring ~ uniform(0, 1);
  pq_sw_fall ~ uniform(0, 1);

  // standard deviation of hunting effort
  sigma_E_fi ~ cauchy(0, 0.1);
  sigma_E_sw_spring ~ cauchy(0, 0.1);
  sigma_E_sw_fall ~ cauchy(0, 0.1);

  // bias vectors for hunting & bycatch (untransformed)
  g_fi ~ normal(0, 1);
  g_sw_spring ~ normal(0, 1);
  g_sw_fall ~ normal(0, 1);
  g_bycatch ~ normal(0, 1);

  // stochastic noise for hunting
  eps_E_fi ~ normal(0,1);
  eps_E_sw_spring ~ normal(0,1);
  eps_E_sw_fall ~ normal(0,1);

  // haul-out
  w_unscaled ~ beta(4, 4);
  delta ~ uniform(0, 1);
  alpha_0_scale_raw ~ uniform(-pi()/2, pi()/2); // uniform re-parametrization of Cauchy prior
  alpha_1_raw ~ uniform(-pi()/2, pi()/2); // uniform re-parametrization of Cauchy prior
  //d ~ uniform(0, 1);
  //f ~ cauchy(0, 37);

  // observation process
  sqrt_r_inv ~ cauchy(0, 1);
  pi_s ~ uniform(0, 1);
  pi_c ~ uniform(0, 1);
  kappa ~ uniform(0, 1);
  scar_na_0 ~ uniform(0, 1);
  scar_na_1 ~ uniform(0, 1);

  // density dependence
  K_unscaled ~ gamma(18, 8);
  theta_0_raw ~ uniform(0, 1);

  // demographic stochasticity
  u_pup ~ normal(0, 1);
  u_sex ~ normal(0, 1);
  
  // Observation models for sampling bias data (see Appendix S5)
  o_16_19[1,1] ~ binomial(o_16_19[1,2], x_16_19[1]);
  o_16_19[2,1] ~ binomial(o_16_19[2,2], x_16_19[2]);
  o_20_21[1,1] ~ binomial(o_20_21[1,2], x_20_21[1]);
  o_20_21[2,1] ~ binomial(o_20_21[2,2], x_20_21[2]);
  for(i in 1:(2*a)) { l[i, 1] ~ binomial(l[i, 2], q[i]); }

  for(i in 1:t) {
    wgt_samples_sw_16_19 = v_16_19 .* N_dead[i][,2];
    wgt_samples_sw_20_21 = v_20_21 .* N_dead[i][1:a,2];
    wgt_samples_bycatch = psi_bycatch .* N_dead[i][,1];

    // demographic stochasticity (cont.)
    if(i <= t0) { u0[,i] ~ normal(0, 1); }
    if(i > t0 && i <= t1) { u1[,i-t0] ~ normal(0, 1); }
    if(i > t1) { u2[,i-t1] ~ normal(0, 1); }

    ///////////////////////////////////// OBSERVATION MODELS ///////////////////////////////////////////

    // aerial surveys
    if(y_survey[i] > 0 && exclude_survey[i] == 0) { y_survey[i] ~ neg_binomial_2(N_w[i], r); }

    /////////// Reproductive assessments (see Appendix S6 for handling of incomplete CA/scar data) ////////////////////

    // pregnancy observations
    if(pregnant_ss[i] > 0 && exclude_rep[i] == 0) { y_pregnant[i] ~ binomial(pregnant_ss[i], p[i]); }

    // outcomes of reproductive sign assessments
    if(sum(z[,i]) > 0 && exclude_rep[i] == 0) { z[,i] ~ multinomial(gamma[,i]); }

    // observations of CA (placental scar not applicable)
    if(CA_ss[i] > 0 && exclude_rep[i] == 0) { y_CA[i] ~ binomial(CA_ss[i], gamma[2,i]+gamma[4,i]+gamma[5,i]); }

    // observations of placental scars (CA evaluation missing)
    if(scar_ss[i] > 0 && exclude_rep[i] == 0) { y_scar[i] ~ binomial(scar_ss[i], (gamma[3,i]+gamma[4,i])/sum(gamma[1:4,i])); }

    ///////////// Hunting bags //////////////////////////////

    if(y_hb_sw[i] > 0 && exclude_hb[i] == 0) {
      if(years[i] < 2017) { // only total hunting bag is known for 2015-2016
        y_hb_sw[i] ~ normal(sum(N_dead[i][,2]+N_dead[i][,3]), 0.05*sum(N_dead[i][,2]+N_dead[i][,3]));
      } else {
        y_hb_sw_spring[i] ~ normal(sum(N_dead[i][,2]), 0.05*sum(N_dead[i][,2]));
        y_hb_sw_fall[i] ~ normal(sum(N_dead[i][,3]), 0.05*sum(N_dead[i][,3]));
      }
    }
    if(y_hb_fi[i] > 0 && exclude_hb[i] == 0) { 
      y_hb_fi[i] ~ normal(sum(N_dead[i][,4]), 0.05*sum(N_dead[i][,4])); 
    }

    ////////////// Hunted and bycaught samples (see Appendix S6 for handling of incomplete age/sex records) ///////////////////////////

    if(exclude_samples[i] == 0) {

      // Swedish hunting samples from the spring
      if(years[i] == 2015) { // random sampling in 2015
        if(sum(y_hs_sw_spring[,i]) > 0) { 
          y_hs_sw_spring[,i] ~ multinomial(N_dead[i][,2]/sum(N_dead[i][,2]));
        }
        if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
          y_hs_sw_sexNA_spring[,i] ~ multinomial((N_dead[i][1:6,2] + N_dead[i][7:12,2])/sum(N_dead[i][,2]));
        }
        if(sum(y_hs_sw_ageNA_spring[,i]) > 0) {
          y_hs_sw_ageNA_spring[1,i] ~ binomial(sum(y_hs_sw_ageNA_spring[,i]), sum(N_dead[i][1:6,2])/sum(N_dead[i][,2]));
        }
      }
      if(years[i] > 2015 && years[i] <= 2019) { // biased sampling between 2015-2019
        if(sum(y_hs_sw_spring[,i]) > 0) {
          y_hs_sw_spring[,i] ~ multinomial(wgt_samples_sw_16_19/sum(wgt_samples_sw_16_19));
        }
        if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
          y_hs_sw_sexNA_spring[,i] ~ multinomial((wgt_samples_sw_16_19[1:a]+wgt_samples_sw_16_19[(a+1):(2*a)])/sum(wgt_samples_sw_16_19));
        }
        if(sum(y_hs_sw_ageNA_spring[,i]) > 0) {
          y_hs_sw_ageNA_spring[1,i] ~ binomial(sum(y_hs_sw_ageNA_spring[,i]), sum(wgt_samples_sw_16_19[1:6])/sum(wgt_samples_sw_16_19));
        }
      }
      if(years[i] > 2019) { // female only sampling in 2020-2021
        if(sum(y_hs_sw_spring[1:6,i]) > 0) {
          y_hs_sw_spring[1:6,i] ~ multinomial(wgt_samples_sw_20_21/sum(wgt_samples_sw_20_21)); 
        }
        if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
          y_hs_sw_sexNA_spring[,i] ~ multinomial(wgt_samples_sw_20_21/sum(wgt_samples_sw_20_21));
        }
      }


      // Swedish hunting samples from the fall
      if(sum(y_hs_sw_fall[,i]) > 0) {
        y_hs_sw_fall[,i] ~ multinomial(N_dead[i][,3]/sum(N_dead[i][,3]));
      }
      if(sum(y_hs_sw_sexNA_fall[,i]) > 0) {
        y_hs_sw_sexNA_fall[,i] ~ multinomial((N_dead[i][1:6,3] + N_dead[i][7:12,3])/sum(N_dead[i][,3]));
      }
      if(sum(y_hs_sw_ageNA_fall[,i]) > 0) {
        y_hs_sw_ageNA_fall[1,i] ~ binomial(sum(y_hs_sw_ageNA_fall[,i]), sum(N_dead[i][1:6,3])/sum(N_dead[i][,3]));
      }

      // Finnish hunting samples
      if(sum(y_hs_fi[,i]) > 0) {
        y_hs_fi[,i] ~ multinomial(N_dead[i][,4]/sum(N_dead[i][,4]));
      }
      if(sum(y_hs_fi_sexNA[,i]) > 0) {
        y_hs_fi_sexNA[,i] ~ multinomial((N_dead[i][1:6,4]+N_dead[i][7:12,4])/sum(N_dead[i][,4]));
      }
      if(sum(y_hs_fi_ageNA[,i]) > 0) {
        y_hs_fi_ageNA[1,i] ~ binomial(sum(y_hs_fi_ageNA[,i]), sum(N_dead[i][1:6,4])/sum(N_dead[i][,4]));
      }

      // Bycaught samples
      if(sum(y_bc[,i]) > 0) {
        y_bc[,i] ~ multinomial(wgt_samples_bycatch / sum(wgt_samples_bycatch));
      }
      if(sum(y_bc_sexNA[,i]) > 0) {
        y_bc_sexNA[,i] ~ multinomial((wgt_samples_bycatch[1:a]+wgt_samples_bycatch[(a+1):(2*a)])/sum(wgt_samples_bycatch));
      }
      if(sum(y_bc_ageNA[,i]) > 0) {
        y_bc_ageNA[,i] ~ binomial(sum(y_bc_ageNA[,i]), sum(wgt_samples_bycatch[1:6])/sum(wgt_samples_bycatch));
      }
    }

  }
  // final year aerial survey estimate
  if(y_survey[t+1] > 0 && exclude_survey[t+1] == 0) { y_survey[t+1] ~ neg_binomial_2(N_w[t+1], r); }

}

generated quantities {

  ///////////////////// Likelihood calculations & Posterior predictive checks //////////////////////////////

  // Swedish spring hunting weighted by sampling bias
  vector[2*a] wgt_samples_sw_16_19;
  vector[a] wgt_samples_sw_20_21;
  
  // natural mortality weighted by bycatch deviation
  vector[2*a] wgt_samples_bycatch;

  // aerial survey estimates
  int<lower=0> y_survey_sim[t+1] = rep_array(0, t+1);
  real log_lik_survey[t+1] = rep_array(0.0, t+1);

  // hunting bags
  real y_hb_fi_sim[t] = rep_array(0.0, t);
  real y_hb_sw_sim[t] = rep_array(0.0, t);
  real log_lik_hb_fi[t] = rep_array(0.0, t);
  real log_lik_hb_sw[t] = rep_array(0.0, t);

  // hunted & bycaught samples
  int<lower=0> y_hs_fi_sim[2*a, t] = rep_array(0, 2*a, t);
  int<lower=0> y_hs_sw_fall_sim[2*a, t] = rep_array(0, 2*a, t);
  int<lower=0> y_hs_sw_spring_sim[2*a, t] = rep_array(0, 2*a, t);
  int<lower=0> y_bc_sim[2*a, t] = rep_array(0, 2*a, t);
  real log_lik_hs_fi[t] = rep_array(0.0, t);
  real log_lik_hs_sw_spring[t] = rep_array(0.0, t);
  real log_lik_hs_sw_fall[t] = rep_array(0.0, t);
  real log_lik_bc[t] = rep_array(0.0, t);

  // reproductive assessments
  int<lower=0> y_pregnant_sim[t] = rep_array(0, t);
  int<lower=0> y_CA_sim[t] = rep_array(0, t);
  int<lower=0> y_scar_sim[t] = rep_array(0, t);
  int<lower=0> z_sim[6, t] = rep_array(0, 6, t);
  int<lower=0> z_non_na_sim[4, t] = rep_array(0, 4, t);
  real log_lik_reproduction[t] = rep_array(0.0, t);
  

  for(i in 1:t) {
    wgt_samples_sw_16_19 = v_16_19 .* N_dead[i][,2];
    wgt_samples_sw_20_21 = v_20_21 .* N_dead[i][1:a,2];
    wgt_samples_bycatch = psi_bycatch .* N_dead[i][,1];
    
    // aerial survey estimate
    y_survey_sim[i] = neg_binomial_2_rng(N_w[i], r);
    if(y_survey[i] > 0) {log_lik_survey[i] = neg_binomial_2_lpmf(y_survey[i] | N_w[i], r);}

    // hunting bags
    if(i > t0) {
      y_hb_sw_sim[i] = normal_rng(sum(N_dead[i][,2]+N_dead[i][,3]), 0.05*sum(N_dead[i][,2]+N_dead[i][,3]));
      log_lik_hb_sw[i] += normal_lpdf(y_hb_sw[i] | sum(N_dead[i][,2]+N_dead[i][,3]), 0.05*sum(N_dead[i][,2]+N_dead[i][,3]));
    }
    if(i > t1) { 
      y_hb_fi_sim[i] = normal_rng(sum(N_dead[i][,4]), 0.05*sum(N_dead[i][,4]));
      log_lik_hb_fi[i] += normal_lpdf(y_hb_fi[i] | sum(N_dead[i][,4]), 0.05*sum(N_dead[i][,4]));
    }

    // Swedish hunting samples from the spring
    if(years[i] == 2015) {
      if(sum(y_hs_sw_spring[,i]) > 0) {
        y_hs_sw_spring_sim[,i] = multinomial_rng(N_dead[i][,2]/sum(N_dead[i][,2]), sum(y_hs_sw_spring[,i]));
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_spring[,i] | N_dead[i][,2]/sum(N_dead[i][,2]));
      }
      if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_sexNA_spring[,i] | (N_dead[i][1:6,2] + N_dead[i][7:12,2])/sum(N_dead[i][,2]));
      }
      if(sum(y_hs_sw_ageNA_spring[,i]) > 0) {
        log_lik_hs_sw_spring[i] += binomial_lpmf(y_hs_sw_ageNA_spring[1,i] | sum(y_hs_sw_ageNA_spring[,i]), sum(N_dead[i][1:6,2])/sum(N_dead[i][,2]));
      }
    }
    if(years[i] > 2015 && years[i] <= 2019) {
      if(sum(y_hs_sw_spring[,i]) > 0) {
        y_hs_sw_spring_sim[,i] = multinomial_rng(wgt_samples_sw_16_19/sum(wgt_samples_sw_16_19), sum(y_hs_sw_spring[,i]));
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_spring[,i] | wgt_samples_sw_16_19/sum(wgt_samples_sw_16_19));
      }
      if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_sexNA_spring[,i] | (wgt_samples_sw_16_19[1:a]+wgt_samples_sw_16_19[(a+1):(2*a)])/sum(wgt_samples_sw_16_19));
      }
      if(sum(y_hs_sw_ageNA_spring[,i]) > 0) {
        log_lik_hs_sw_spring[i] += binomial_lpmf(y_hs_sw_ageNA_spring[1,i] | sum(y_hs_sw_ageNA_spring[,i]), sum(wgt_samples_sw_16_19[1:6])/sum(wgt_samples_sw_16_19));
      }
    }
    if(years[i] > 2019) {
      if(sum(y_hs_sw_spring[1:a,i]) > 0) {
        y_hs_sw_spring_sim[1:a,i] = multinomial_rng(wgt_samples_sw_20_21/sum(wgt_samples_sw_20_21), sum(y_hs_sw_spring[1:a,i]));
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_spring[1:a,i] | wgt_samples_sw_20_21/sum(wgt_samples_sw_20_21));
      }
      if(sum(y_hs_sw_sexNA_spring[,i]) > 0) {
        log_lik_hs_sw_spring[i] += multinomial_lpmf(y_hs_sw_sexNA_spring[,i] | wgt_samples_sw_20_21/sum(wgt_samples_sw_20_21));
      }
    }


    // Swedish hunting samples from the fall
    if(sum(y_hs_sw_fall[,i]) > 0) {
      y_hs_sw_fall_sim[,i] = multinomial_rng(N_dead[i][,3]/sum(N_dead[i][,3]), sum(y_hs_sw_fall[,i]));
      log_lik_hs_sw_fall[i] += multinomial_lpmf(y_hs_sw_fall[,i] | N_dead[i][,3]/sum(N_dead[i][,3]));
    }
    if(sum(y_hs_sw_sexNA_fall[,i]) > 0) {
      log_lik_hs_sw_fall[i] += multinomial_lpmf(y_hs_sw_sexNA_fall[,i] | (N_dead[i][1:6,3] + N_dead[i][7:12,3])/sum(N_dead[i][,3]));
    }
    if(sum(y_hs_sw_ageNA_fall[,i]) > 0) {
      log_lik_hs_sw_fall[i] += binomial_lpmf(y_hs_sw_ageNA_fall[1,i] | sum(y_hs_sw_ageNA_fall[,i]), sum(N_dead[i][1:6,3])/sum(N_dead[i][,3]));
    }

    // Finnish hunting samples
    if(sum(y_hs_fi[,i]) > 0) { 
      y_hs_fi_sim[,i] = multinomial_rng(N_dead[i][,4]/sum(N_dead[i][,4]), sum(y_hs_fi[,i]));
      log_lik_hs_fi[i] += multinomial_lpmf(y_hs_fi[,i] | N_dead[i][,4]/sum(N_dead[i][,4]));
    }
    if(sum(y_hs_fi_sexNA[,i]) > 0) {
      log_lik_hs_fi[i] += multinomial_lpmf(y_hs_fi_sexNA[,i] | (N_dead[i][1:6,4]+N_dead[i][7:12,4])/sum(N_dead[i][,4]));
    }
    if(sum(y_hs_fi_ageNA[,i]) > 0) {
      log_lik_hs_fi[i] += binomial_lpmf(y_hs_fi_ageNA[1,i] | sum(y_hs_fi_ageNA[,i]), sum(N_dead[i][1:6,4])/sum(N_dead[i][,4]));
    }

    // bycaught samples
    if(sum(y_bc[,i]) > 0) {
      y_bc_sim[,i] = multinomial_rng(wgt_samples_bycatch/sum(wgt_samples_bycatch), sum(y_bc[,i]));
      log_lik_bc[i] += multinomial_lpmf(y_bc[,i] | wgt_samples_bycatch/sum(wgt_samples_bycatch));
    }
    if(sum(y_bc_sexNA[,i]) > 0) {
      log_lik_bc[i] += multinomial_lpmf(y_bc_sexNA[,i] | (wgt_samples_bycatch[1:a]+wgt_samples_bycatch[(a+1):(2*a)])/sum(wgt_samples_bycatch));
    }
    if(sum(y_bc_ageNA[,i]) > 0) {
      log_lik_bc[i] += binomial_lpmf(y_bc_ageNA[,i] | sum(y_bc_ageNA[,i]), sum(wgt_samples_bycatch[1:6])/sum(wgt_samples_bycatch));
    }

    // reproductive assessments
    if(pregnant_ss[i] > 0) {
      y_pregnant_sim[i] = binomial_rng(pregnant_ss[i],  p[i]);
      log_lik_reproduction[i] += binomial_lpmf(y_pregnant[i] | pregnant_ss[i],  p[i]);
    }
    if(CA_ss[i] > 0) {
      y_CA_sim[i] = binomial_rng(CA_ss[i],  gamma[2,i]+gamma[4,i]+gamma[5,i]);
      log_lik_reproduction[i] += binomial_lpmf(y_CA[i] | CA_ss[i],  gamma[2,i]+gamma[4,i]+gamma[5,i]);
    }
    if(scar_ss[i] > 0) {
      y_scar_sim[i] = binomial_rng(scar_ss[i],  (gamma[3,i]+gamma[4,i])/sum(gamma[1:4,i]));
      log_lik_reproduction[i] += binomial_lpmf(y_scar[i] | scar_ss[i],  (gamma[3,i]+gamma[4,i])/sum(gamma[1:4,i]));
    }
    if(sum(z[,i]) > 0) {
      z_sim[,i] = multinomial_rng(gamma[,i], sum(z[,i]));
      log_lik_reproduction[i] += multinomial_lpmf(z[,i] | gamma[,i]);
    }
    if(sum(z[1:4,i]) > 0) {z_non_na_sim[,i] = multinomial_rng(gamma[1:4,i]/sum(gamma[1:4,i]), sum(z[1:4,i]));}

  }
  // aerial survey estimate (final year)
  y_survey_sim[t+1] = neg_binomial_2_rng(N_w[t+1], r);
  log_lik_survey[t+1] = neg_binomial_2_lpmf(y_survey[t+1] | N_w[t+1], r);

}
