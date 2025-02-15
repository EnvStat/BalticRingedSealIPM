
# set font and figure sizes
font_size <- 25
size <- 8

# initialize list to store posterior samples from all model fits
samples <- list()

# full model posterior samples
model <- readRDS('ringed_seal_results/ringed_seal_post.rds')
samples[[1]] <- rstan::extract(model)

# reduced dataset samples
for(i in 1:length(sensitivity_data)) {
  output_directory <- paste0('ringed_seal_results/sensitivities/',i,'/')
  model <- readRDS(paste0(output_directory, 'ringed_seal_post_sens',i,'.rds'))
  samples[[i+1]] <- rstan::extract(model)
}

# prior and model sensitivity samples
for(i in 1:4) {
  output_directory <- paste0('ringed_seal_results/sensitivities/m',i-1,'/')
  model <- readRDS(paste0(output_directory, 'ringed_seal_post_sens_m',i-1,'.rds'))
  samples[[1+length(sensitivity_data)+i]] <- rstan::extract(model)
}

# function for plotting prior-posterior comparison and compute overlaps
prior_posterior <- function(x, prior_dens, posterior_samples, xlim,
                            xlab = "x", label_x = 3/4, label_y = 3.5/4) {
  
  idx <- 1:8 # indices for models to include in plots

  lb <- min(x); ub <- max(x); n <- length(x);
  model_labels <- c('Full model', 'Survey (thinned)', 'Survey (< 2012)', 'Survey (>= 2012)',
                    'Reproduction (thinned)', 'Reproduction (> 2006)', 'Hunting bag', 'Dead samples')
  
  # data frame for storing summary statistics
  summaries <- data.frame(sensitivity = factor(c(1:ncol(posterior_samples))), 
                         prior_overlap = rep(NA, ncol(posterior_samples)),
                         posterior_overlap = rep(NA, ncol(posterior_samples)),
                         median = rep(NA, ncol(posterior_samples)),
                         sd = rep(NA, ncol(posterior_samples)))
  
  # data frame for storing densities (initialized with prior density)
  df <- data.frame(sensitivity = factor(rep(0, n)), x = x, y = prior_dens)
  
  # compute posterior density for full model
  posterior_main <- density(posterior_samples[,1], from = lb, to = ub, n = n, na.rm = T)$y
  
  # compute posterior densities and summary statistics for other models
  for(i in 1:ncol(posterior_samples)) {
    posterior_dens <- density(posterior_samples[,i], from = lb, to = ub, n = n, na.rm = T)$y
    summaries$prior_overlap[i] <- round(1-sum(pmax(posterior_dens - prior_dens, 0)*(ub-lb)/n), 2)*100
    summaries$posterior_overlap[i] <- round(1-sum(pmax(posterior_dens - posterior_main, 0)*(ub-lb)/n), 2)*100
    summaries$median[i] <- round(median(posterior_samples[,i], na.rm = T), 2)
    summaries$sd[i] <- round(sd(posterior_samples[,i], na.rm = T), 2)
    df.post <- data.frame(sensitivity = factor(rep(i, n)), x = x, y = posterior_dens)
    df <- rbind(df, df.post)
  }
  
  # y-coordinates for labels
  y_coords <- df$y[df$sensitivity %in% idx] # posterior densities for data sens. to adjust y-axis scale
  label_ys <- label_y*max(c(prior_dens, y_coords)) - seq(0, 0.4*max(c(prior_dens, y_coords)), length.out = length(idx))

  # filter data for plotting
  df.plot <- df[which(df$sensitivity %in% c(0,idx)),]
  summaries.plot <- summaries[which(summaries$sensitivity %in% idx),]
  
  # print plot
  print(
    ggplot(df.plot, aes(x = x, y = y, color = sensitivity)) +
      geom_line(aes(size = sensitivity, linetype = sensitivity)) +
      geom_text(data = summaries.plot, show.legend = F, size = 6, aes(label = paste0(model_labels[idx], ' = ', prior_overlap, '%'), 
                                     x = xlim[1] + label_x*(xlim[2]-xlim[1]), 
                                     y = label_ys,
                                     color = sensitivity)) +
      scale_linetype_manual(values=c("dashed", rep("solid", length(idx)))) +
      scale_size_manual(values=c(1, 1, rep(0.5, length(idx)-1))) +
      scale_color_manual(values=c("black", "black", "#0072B2","#56B4E9","#009E73",
                                  "#D55E00", "#E69F00", "#8E2071", "#f36ccf")) +
      xlim(xlim[1], xlim[2]) +
      labs(x = xlab, y = 'Density') +
      theme_classic() +
      theme(text = element_text(size=font_size), legend.position = 'none')
    )
  
  return(summaries)
}

# directory for saving outputs
output_directory <- 'ringed_seal_results/sensitivities/prior_posterior/'

n_samps <- nrow(samples[[1]]$n0) # number of posterior draws

# posterior samples to exclude
excl <- matrix(FALSE, nrow = n_samps, ncol = length(samples))
excl[1:(n_samps/4),6:7] <- TRUE # remove chain 1 due to poor mixing
excl[(1/2*n_samps):(3/4*n_samps),8] <- TRUE # remove chain 3 due to poor mixing

# initial population size
lb <- 0; ub <- 30000; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- dlnorm(x, 8.25, 1)

posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$n0
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(0, 12500), expression(n[1988]))
ggsave(filename='n0.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'n0.csv'))

# survival parameters
lb <- 0.8; ub <- 0.97; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0.8, 0.97)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$phi_f5
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(phi['f,5+']), label_x = 1/4)
ggsave(filename='phi_5.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'phi_5.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$phi_scale
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(phi['f,0']/phi['f,5+']), label_x = 0.8/4, label_y = 4/4)
ggsave(filename='phi_0_5.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'phi_0_5.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$c
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), "c", label_x = 2.5/4, label_y = 1.35/2)
ggsave(filename='c.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'c.csv'))

lb <- -0.4; ub <- 0.6; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dnorm(x, 0, 0.1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$nu[,1]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(nu['0']), label_x = 3.25/4)
ggsave(filename='nu_0.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'nu_0.csv'))

lb <- -0.4; ub <- 0.6; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dnorm(x, 0, 0.1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$nu[,2]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(nu['1+']), label_x = 3.25/4)
ggsave(filename='nu_1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'nu_1.csv'))


# hunting parameters
lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$pq_fi
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(p[Q]^'fi'), label_x = 1/4)
ggsave(filename='p_Q_fi.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'p_Q_fi.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$pq_sw_spring
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(p[Q]^'sw(1)'))
ggsave(filename='p_Q_sw1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'p_Q_sw1.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$pq_sw_fall
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(p[Q]^'sw(2)'))
ggsave(filename='p_Q_sw2.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'p_Q_sw2.csv'))

lb <- 0; ub <- 20; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dcauchy(x, scale = 0.1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$sigma_E_fi
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(0, 2), expression(sigma[E^'fi']))
ggsave(filename='sigma_E_fi.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'sigma_E_fi.csv'))

lb <- 0; ub <- 2.5; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dcauchy(x, scale = 0.1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$sigma_E_sw_spring
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(sigma[E^'sw(1)']))
ggsave(filename='sigma_E_sw1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'sigma_E_sw1.csv'))

lb <- 0; ub <- 20; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dcauchy(x, scale = 0.1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$sigma_E_sw_fall
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, 3), expression(sigma[E^'sw(2)']))
ggsave(filename='sigma_E_sw2.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'sigma_E_sw2.csv'))

n_draws <- 100000
draws <- matrix(rnorm(n_draws*ncol(data$L_bias), 0, sd=1), ncol(data$L_bias), n_draws)
prior_samples_raw <- data$L_bias%*%draws
prior_samples <- cbind(t(prior_samples_raw), rep(0, n_draws))
prior_samples <- apply(prior_samples, 1, softmax)

lb <- 0; ub <- 0.4; n <- 500;
x <- seq(lb, ub, length.out = n)

classes <- paste0(rep(c('f,', 'm,'), each=6), c(0:4,"5+"))
for(j in 1:12) {
  prior <- density(prior_samples[j,], from = lb, to = ub, n = n)$y
  posterior_samples <- matrix(NA, n_samps, length(samples))
  for(i in 1:length(samples)) {
    posterior_samples[,i] <- samples[[i]]$psi_fi[,j]
    posterior_samples[excl[,i],i] <- NA
  }
  output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), substitute(psi[j]^{'fi'}, list(j=classes[j])))
  ggsave(filename=paste0('psi_fi_',j,'.png'), width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
  write.csv(output, file = paste0(output_directory, 'psi_fi_',j,'.csv'))
  
  posterior_samples <- matrix(NA, n_samps, length(samples))
  for(i in 1:length(samples)) {
    posterior_samples[,i] <- samples[[i]]$psi_sw_spring[,j]
    posterior_samples[excl[,i],i] <- NA
  }
  output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), substitute(psi[j]^{'sw(1)'}, list(j=classes[j])))
  ggsave(filename=paste0('psi_sw1_',j,'.png'), width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
  write.csv(output, file = paste0(output_directory, 'psi_sw1_',j,'.csv'))
  
  posterior_samples <- matrix(NA, n_samps, length(samples))
  for(i in 1:length(samples)) {
    posterior_samples[,i] <- samples[[i]]$psi_sw_fall[,j]
    posterior_samples[excl[,i],i] <- NA
  }  
  output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), substitute(psi[j]^{'sw(2)'}, list(j=classes[j])))
  ggsave(filename=paste0('psi_sw2_',j,'.png'), width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
  write.csv(output, file = paste0(output_directory, 'psi_sw2_',j,'.csv'))
}

lb <- 0; ub <- 0.6; n <- 500;
x <- seq(lb, ub, length.out = n)
for(j in 1:12) {
  prior <- density(prior_samples[j,], from = lb, to = ub, n = n)$y
  posterior_samples <- matrix(NA, n_samps, length(samples))
  for(i in 1:length(samples)) {
    posterior_samples[,i] <- samples[[i]]$psi_bycatch[,j]
    posterior_samples[excl[,i],i] <- NA
  }  
  output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), substitute(psi[j]^{'bc'}, list(j=classes[j])))
  ggsave(filename=paste0('psi_bc_',j,'.png'), width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
  write.csv(output, file = paste0(output_directory, 'psi_bc_',j,'.csv'))
}


# reproduction parameters
lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$b_max
  posterior_samples[excl[,i],i] <- NA
} 
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(b['max']), label_x = 1/4)
ggsave(filename='b_max.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'b_max.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$b_scale
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(b['min']/b['max']), label_x = 1/4, label_y = 1/4)
ggsave(filename='b_min_max.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'b_min_max.csv'))

lb <- -50; ub <- 50; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- dcauchy(x, scale = 10)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$beta[,1]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(-15, 15), expression(beta[0]), label_x = 1/4)
ggsave(filename='beta_0.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'beta_0.csv'))

lb <- -20; ub <- 20; n <- 1000;
x <- seq(lb, ub, length.out = n)
prior <- 1/1.8*dt(x/1.8, 5)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$beta[,2]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(-12, 10), expression(beta[1]), label_x = 1/4)
ggsave(filename='beta_1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'beta_1.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$pi_s
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(pi[s]), label_x = 1.75/4)
ggsave(filename='pi_s.png', width=0.65*size, height=size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'pi_s.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$pi_c
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(pi[c]), label_x = 1.75/4)
ggsave(filename='pi_c.png', width=0.65*size, height=size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'pi_c.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$kappa
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(kappa), label_x = 2.5/4)
ggsave(filename='kappa.png', width=0.65*size, height=size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'kappa.csv'))


# haul-out parameters
lb <- 0.4; ub <- 0.9; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dbeta((x-0.4)/0.5, 4, 4)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$w
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, 1), expression(hat(omega)), label_x = 3.25/4, label_y = 4/4)
ggsave(filename='w.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'w.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$delta
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(delta), label_x = 1/4, label_y = 4/4)
ggsave(filename='delta.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'delta.csv'))

lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$d
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), "d")
ggsave(filename='d.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'd.csv'))

lb <- 0; ub <- 100; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dcauchy(x, scale=37)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$f
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), "f")
ggsave(filename='f.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'f.csv'))

lb <- -70; ub <- 70; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- dcauchy(x, scale = 1.25)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$alpha[,2]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(-20, 40), expression(alpha[1]))
ggsave(filename='alpha_1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'alpha_1.csv'))

lb <- -20; ub <- 20; n <- 1500;
x <- seq(lb, ub, length.out = n)
prior <- dcauchy(x, scale = 5)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$alpha[,1]/samples[[i]]$alpha[,2]
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(-1, 3), expression(alpha[0]/alpha[1]), label_x = 1/4)
ggsave(filename='alpha_0_1.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'alpha_0_1.csv'))

# observation parameters
lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- 2*dcauchy(x, scale = 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$sqrt_r_inv
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), expression(1/sqrt(r)))
ggsave(filename='r.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'r.csv'))


# density dependence parameters
lb <- 0; ub <- 1; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dunif(x, 0, 1)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$theta_0_raw
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, ub), paste(expression(theta[0]), '(untransformed)'), 
                label_x = 1/4, label_y = 1/2.5)
ggsave(filename='theta_0.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'theta_0.csv'))

lb <- 0; ub <- 5e5; n <- 500;
x <- seq(lb, ub, length.out = n)
prior <- dgamma(x, 18, 8e-5)
posterior_samples <- matrix(NA, n_samps, length(samples))
for(i in 1:length(samples)) {
  posterior_samples[,i] <- samples[[i]]$K
  posterior_samples[excl[,i],i] <- NA
}
output <- prior_posterior(x, prior, posterior_samples, c(lb, 1.05*ub), expression(K), label_x = 3.15/4, label_y = 4/4)
ggsave(filename='K.png', width=size, height=0.65*size, path=paste0(output_directory, 'plots/'))
write.csv(output, file = paste0(output_directory, 'K.csv'))


