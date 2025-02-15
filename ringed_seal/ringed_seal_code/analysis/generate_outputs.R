

###################################################################################################### 
###################################### GENERATE PLOTS ###############################################
###################################################################################################### 

######################## Plot total population size ############################

label_k = function(x, s) {
  labels = paste(x, s, sep='')
  labels[x==0] <- 0
  return(labels)
}
N.tot.samples <- apply(samples$N, c(1,3), sum)
N.quantiles <- t(apply(N.tot.samples, 2, quantile, quantiles))

K.scale <- 0.12 # scaling of second y-axis for the carrying capacity
ggplot() +
  geom_ribbon(aes(x=c(data$years, 2023), ymin=N.quantiles[,1], ymax=N.quantiles[,3], fill='gray16'), alpha=0.25) +
  geom_line(aes(x=c(data$years, 2023), y=N.quantiles[,2]), linewidth=1) +
  geom_ribbon(aes(x=c(2025,2026), ymin=0, ymax=0, fill='#369ACC')) +
  geom_ribbon(aes(x=c(2025,2026), 
                  ymin=K.scale*quantile(samples$K, quantiles[1]), ymax=K.scale*quantile(samples$K, quantiles[3]), 
                  fill='#75CE87',), alpha=0.75) +
  geom_line(aes(x=c(2025,2026), y=K.scale*median(samples$K)), linewidth=1, color='#005612') +
  scale_y_continuous(limits=c(0,4.2e4), breaks=seq(0,4e4,5e3), label = function(x){label_k(1e-3*x, 'e3')}, 
    sec.axis = sec_axis(~ . * 1/K.scale, breaks=seq(0,4e5,5e4), label = function(x){label_k(1e-3*x, 'e3')}, name = "")) +
  scale_x_continuous(limits=c(1987.5,2026), breaks=c(seq(1990,2020,5),2025), expand=c(0.02,0), labels=c(seq(1990,2020,5), expression(infinity))) +
  labs(x='Year', y='Number of ringed seals') +
  scale_fill_identity(name = 'Asymptotic projections', guide = 'legend',
                      labels = c('without density dependence', 'with density dependence', 'historical development')) +
  guides(fill = guide_legend(reverse = T)) +
  theme_classic() +
  theme(legend.position = c(0.3,0.8), text = element_text(size=font_size), 
        axis.text.y.right = element_text(colour = "#005612"))
ggsave(filename='pop_size.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))


##################### Plot growth rate ################################

# historical growth rate
growth.rates <- (N.tot.samples[,2:ncol(N.tot.samples)] / N.tot.samples[,1:(ncol(N.tot.samples)-1)] - 1)*100
growth.rate.quantiles <- t(apply(growth.rates, 2, quantile, quantiles))

# asymptotic growth rate
g_asymptotic <- sapply(c(1:nrow(samples$n0)), function(i) {
  return(lambda(samples$phi[i,], samples$b_max[i]))
})
g_asymptotic_quantiles <- quantile(100*(g_asymptotic-1), quantiles)

ggplot() +
  geom_ribbon(aes(x=c(data$years[-1], 2023)-0.5, ymin=growth.rate.quantiles[,1], ymax=growth.rate.quantiles[,3], fill='gray16'), alpha=0.25) +
  geom_line(aes(x=c(data$years[-1], 2023)-0.5, y=growth.rate.quantiles[,2]), linewidth=1) +
  geom_ribbon(aes(x=c(2024,2025), ymin=g_asymptotic_quantiles[1], ymax=g_asymptotic_quantiles[3], fill='#369ACC'), alpha=0.75) +
  geom_line(aes(x=c(2024,2025), y=g_asymptotic_quantiles[2]), linewidth=1, color='#006598') +
  geom_ribbon(aes(x=c(2025,2026), ymin=-0.1, ymax=0.1, fill='#75CE87'), alpha=0.75) +
  geom_line(aes(x=c(2025,2026), y=0), linewidth=1, color='#005612') +
  geom_hline(aes(yintercept=0), linetype='dashed') +
  scale_y_continuous(limits=c(-2,10)) + 
  scale_x_continuous(limits=c(1987.5,2031), breaks=c(seq(1990,2020,5),2025), labels=c(seq(1990,2020,5), expression(infinity))) +
  labs(x='Year', y='% growth') +
  scale_fill_identity() +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='growth_rate.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))


###################### Plot survival probabilities #########################

# trim 0.5% off each tail for better aesthetics
phi <- apply(samples$phi, 2, function(phi.i) {
  phi.i[which(phi.i < quantile(phi.i, 0.005) | phi.i > quantile(phi.i, 0.995))] <- NA
  return(phi.i)
})
df.phi <- data.frame(sex=factor(rep(c('Female', 'Male'), each=6*nrow(samples$phi))), 
                   age = rep(rep(factor(c(0:4, '5+')), each=nrow(samples$phi)), 2), 
                   phi = c(phi))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.phi, aes(x=age, y=phi, fill=sex, split=sex)) +
  geom_half_violin(position='identity', scale='width', adjust=1.75, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Survival probability') +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='survival_rates.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))

samps.phi.thinned <- samples$phi[seq(1, nrow(samples$phi), length.out = 3000),1:6]
df.phi.cor <- data.frame(age = rep(0:5, each = nrow(samps.phi.thinned)), 
                         phi = c(samps.phi.thinned),
                         sample = rep(1:nrow(samps.phi.thinned), 6),
                         col = samps.phi.thinned[,1])
ggplot(df.phi.cor, aes(x = age, y = phi, group = sample, col=col)) +
  geom_line(linewidth = 0.025) +
  scale_x_continuous(breaks = 0:5, labels = c(0:4, '5+')) +
  labs(x = 'Age', y = 'Survival probability') +
  scale_color_viridis_c(option = "rocket") +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='survival_corr.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))

########################## Plot birth rate ####################################

b.quantiles <- t(apply(samples$b, 2, quantile, quantiles))

# birth rate at carrying capacity (from Euler-Lotka eq.)
b_eq <- 2*(1-samples$phi[,6])/exp(rowSums(log(samples$phi[,1:5])))
b_eq.quantiles <- quantile(b_eq, quantiles)

ggplot() +
  geom_ribbon(aes(x=c(data$years,2023), ymin=b.quantiles[,1], ymax=b.quantiles[,3]), alpha=0.25) +
  geom_line(aes(x=c(data$years,2023), y=b.quantiles[,2]), linewidth=1) +
  geom_ribbon(aes(x=c(2024,2025), 
                  ymin=quantile(samples$b_max, quantiles[1]), ymax=quantile(samples$b_max, quantiles[3])), 
              fill='#369ACC', alpha=0.75) +
  geom_line(aes(x=c(2024,2025), y=median(samples$b_max)), linewidth=1, color='#006598') +
  geom_ribbon(aes(x=c(2025,2026), 
                  ymin=b_eq.quantiles[1], ymax=b_eq.quantiles[3]), 
              fill='#75CE87', alpha=0.75) +
  geom_line(aes(x=c(2025,2026), y=b_eq.quantiles[2]), linewidth=1, color='#005612') + 
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
  scale_x_continuous(limits=c(1987.5,2026), breaks=c(seq(1990,2020,5),2025), labels=c(seq(1990,2020,5), expression(infinity))) +
  labs(x='Year', y='Birth rate') +
  theme_classic() +
  theme(text = element_text(size=font_size))
ggsave(filename='birth_rate.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))


################### Plot Finnish hunting bias #########################

# calculate relative probability of becoming hunted in Finland (see Appendix S3)
tau_fi <- mean(data$tau_fi[data$tau_fi>0]) # average timing of hunting
psi.fi <- samples$psi_fi * samples$phi^tau_fi * samples$w_vec
psi.fi <- psi.fi / rowSums(psi.fi)

# trim 0.5% off the tails for improved aesthetics
psi.fi.trim <- apply(psi.fi, 2, function(psi.i) {
  psi.i[which(psi.i < quantile(psi.i, 0.005) | psi.i > quantile(psi.i, 0.995))] <- NA
  return(psi.i)
})
df.bias.fi <- data.frame(sex=as.factor(rep(c('Female', 'Male'), each=6*nrow(samples$psi_fi))), 
                         age_group = factor(rep(rep(c(0:4, '5+'), each=nrow(samples$psi_fi)), 2), levels=c(0:4, '5+')), 
                         psi = c(psi.fi.trim))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.bias.fi, aes(x=age_group, y=psi, fill=sex, split=sex)) +
  geom_hline(aes(yintercept=1/12), linetype='dashed') +
  geom_half_violin(position='identity', scale='width', adjust=1.5, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Relative probability') +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='h_bias_fi.png', width=6.5, height=6, path=paste(output_directory, 'plots/', sep=''))


################### Plot Swedish hunting bias in spring #########################

# calculate relative probability of becoming hunted in Sweden in spring (see Appendix S3)
tau_sw_spring <- mean(data$tau_sw_spring[data$tau_sw_spring>0]) # average timing of hunting
psi.sw.spring <- samples$psi_sw_spring * samples$phi^tau_sw_spring
psi.sw.spring <- psi.sw.spring / rowSums(psi.sw.spring)

# trim 0.5% off the tails for improved aesthetics
psi.sw.spring.trim <- apply(psi.sw.spring, 2, function(psi.i) {
  psi.i[which(psi.i < quantile(psi.i, 0.005) | psi.i > quantile(psi.i, 0.995))] <- NA
  return(psi.i)
})
df.bias.sw.spring <- data.frame(sex=as.factor(rep(c('Female', 'Male'), each=6*nrow(samples$psi_sw_spring))), 
                         age_group = factor(rep(rep(c(0:4, '5+'), each=nrow(samples$psi_sw_spring)), 2), levels=c(0:4, '5+')), 
                         psi = c(psi.sw.spring.trim))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.bias.sw.spring, aes(x=age_group, y=psi, fill=sex, split=sex)) +
  geom_hline(aes(yintercept=1/12), linetype='dashed') +
  geom_half_violin(position='identity', scale='width', adjust=1.25, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Relative probability') +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='h_bias_sw_spring.png', width=6.5, height=6, path=paste(output_directory, 'plots/', sep=''))


################### Plot Swedish sampling bias in the spring (see Appendix S5) #########################

# compute relative sampling probabilities
v.sw <- samples$v_16_19 / rowSums(samples$v_16_19)

# trim 0.5% off the tails for improved aesthetics
v.sw <- apply(v.sw, 2, function(v.i) {
  v.i[which(v.i < quantile(v.i, 0.0075) | v.i > quantile(v.i, 0.9925))] <- NA
  return(v.i)
})
df.v.sw <- data.frame(Sex=as.factor(rep(c('Female', 'Male'), each=6*nrow(samples$v_16_19))), 
                              age_group = factor(rep(rep(c(0:4, '5+'), each=nrow(samples$v_16_19)), 2), levels=c(0:4, '5+')), 
                              v = c(v.sw))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.v.sw, aes(x=age_group, y=v, fill=Sex, split=Sex)) +
  geom_hline(aes(yintercept=1/12), linetype='dashed') +
  geom_half_violin(position='identity', scale='width', adjust=1.25, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Relative sampling probability') +
  scale_y_continuous(limits = c(0, 0.13), breaks=seq(0,0.13,0.03)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(legend.position = c(0.8,0.15),text = element_text(size=font_size))
ggsave(filename='sampling_bias_sw.png', width=6.5, height=6, path=paste(output_directory, 'plots/', sep=''))


################### Plot Swedish hunting bias in fall #########################

# calculate relative probability of becoming hunted in Sweden in the fall (see Appendix S3)
tau_sw_fall <- mean(data$tau_sw_fall[data$tau_sw_fall>0]) # average timing of hunting
psi.sw.fall <- samples$psi_sw_fall * samples$phi^tau_sw_fall
psi.sw.fall <- psi.sw.fall / rowSums(psi.sw.fall)

# trim 0.5% off the tails for improved aesthetics
psi.sw.fall.trim <- apply(psi.sw.fall, 2, function(psi.i) {
  psi.i[which(psi.i < quantile(psi.i, 0.0075) | psi.i > quantile(psi.i, 0.9925))] <- NA
  return(psi.i)
})
df.bias.sw.fall <- data.frame(sex=as.factor(rep(c('Female', 'Male'), each=6*nrow(samples$psi_sw_fall))), 
                                age_group = factor(rep(rep(c(0:4, '5+'), each=nrow(samples$psi_sw_fall)), 2), levels=c(0:4, '5+')), 
                                psi = c(psi.sw.fall.trim))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.bias.sw.fall, aes(x=age_group, y=psi, fill=sex, split=sex)) +
  geom_hline(aes(yintercept=1/12), linetype='dashed') +
  geom_half_violin(position='identity', scale='width', adjust=1.25, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Relative probability') +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(text = element_text(size=font_size), legend.position = 'none')
ggsave(filename='h_bias_sw_fall.png', width=6.5, height=6, path=paste(output_directory, 'plots/', sep=''))


##################### Plot bycatch bias ###############################

# calculate relative probability of becoming bycaught (see Appendix S3)
psi.bycatch <- samples$psi_bycatch*(1-samples$phi)
psi.bycatch <- psi.bycatch / rowSums(psi.bycatch)

# trim 0.5% off the tails for improved aesthetics
psi.bycatch.trimmed <- apply(psi.bycatch, 2, function(psi.i) {
  psi.i[which(psi.i < quantile(psi.i, 0.005) | psi.i > quantile(psi.i, 0.995))] <- NA
  return(psi.i)
})
df.bias.bycatch <- data.frame(Sex=as.factor(rep(c('Female', 'Male'), each=6*nrow(samples$psi_bycatch))), 
                              age_group = factor(rep(rep(c(0:4, '5+'), each=nrow(samples$psi_bycatch)), 2), levels=c(0:4, '5+')), 
                              psi = c(psi.bycatch.trimmed))

group.colors <- c(Female='#E43B57', Male='#198DCE')
ggplot(df.bias.bycatch, aes(x=age_group, y=psi, fill=Sex, split=Sex)) +
  geom_hline(aes(yintercept=1/12), linetype='dashed') +
  geom_half_violin(position='identity', scale='width', adjust=1.5, trim=F) +
  geom_boxplot(position=position_dodge(0.25), outlier.shape = NA, coef=0, width=0.11, fill='ivory2') +
  labs(x = 'Age', y = 'Relative probability') +
  scale_y_continuous(limits = c(0, 0.55), breaks=seq(0,0.55,0.1)) +
  scale_fill_manual(values=group.colors) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8),text = element_text(size=font_size))
ggsave(filename='bycatch_bias.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))

######################## Plot age structure ##############################

# mean historical age structure
demographic.str <- aperm(samples$N, c(1,3,2)) / c(N.tot.samples)
age.str <- demographic.str[,,1:6] + demographic.str[,,7:12]
age.str.mean <- t(apply(age.str, c(2,3), mean))

# stable age structure during asymptotic exponential growth
stable.str.r.samps <- t(sapply(1:nrow(samples$n0), function(i) {
  stable.str <- matrix(stable.structure.eq(samples$phi[i,]/g_asymptotic[i]), ncol=2)
  return(rowSums(stable.str))}))
stable.str.r <- apply(stable.str.r.samps, 2, mean)

# stable age structure at carrying capacity
stable.str.K.samps <- t(sapply(1:nrow(samples$n0), function(i) {
  stable.str <- matrix(stable.structure.eq(samples$phi[i,]), ncol=2)
  return(rowSums(stable.str))}))
stable.str.K <- apply(stable.str.K.samps, 2, mean)

df.age.str <- data.frame(group=as.factor(rep(1, 6*(length(data$years)+1))), year=rep(c(data$years,2023), 6), 
                         age=rep(c(0:4,'5+'), each=1+length(data$years)),
                         proportion=c(t(age.str.mean)))
df.age.str.r <- data.frame(group=as.factor(rep(2, 6*(length(data$years)+1))), year=rep(c(2024,2025), each=6),
                           age=rep(c(0:4,'5+'), 2),
                           proportion=rep(stable.str.r, 2))
df.age.str.K <- data.frame(group=as.factor(rep(3, 6*(length(data$years)+1))), year=rep(c(2025,2026), each=6),
                           age=rep(c(0:4,'5+'), 2),
                           proportion=rep(stable.str.K, 2))
ggplot(df.age.str, aes(x=year, y=proportion, group=rev(age), alpha=rev(age))) +
  geom_area(col='black') +
  geom_area(data=df.age.str.r, aes(x=year, y=proportion, group=rev(age), alpha=rev(age)), fill='#0572A9', col='black') +
  geom_area(data=df.age.str.K, aes(x=year, y=proportion, group=rev(age), alpha=rev(age)), fill='#005612', col='black') +
  geom_vline(aes(xintercept=2025), color='white', linewidth=0.5) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  scale_x_continuous(limits=c(1987.5, 2026), breaks=c(seq(1990,2020,5),2025), labels=c(seq(1990,2020,5), expression(infinity))) +
  labs(x = 'Year', y='Relative frequency') +
  scale_alpha_discrete(name = 'Age', guide = 'legend') +
  guides(alpha = guide_legend(
    override.aes = list(fill = "gray16", alpha=rev(c(0.05, 0.1, 0.2, 0.35, 0.5, 1))), 
    nrow=1, title.position="top", label.position = "bottom")) +
  theme_classic() +
  theme(legend.position = c(0.23,0.75), legend.direction = "horizontal", text = element_text(size=font_size), 
        legend.background = element_rect(fill = NA), legend.key = element_rect(color = 'black'))
ggsave(filename='age_str.png', width=10, height=5, path=paste(output_directory, 'plots/', sep=''))


################## Plot expected haul-out probability over time ##############################

w_t <- samples$N_w / c(N.tot.samples)#[,1:35])
w_t_quantiles <- t(apply(w_t, 2, quantile, quantiles))
error_bars(w_t_quantiles, shape=8, cex=2, alpha = 1,
           x_ticks=c(data$years, data$years[length(data$years)]+1), y_data=w_t_quantiles[,2], 
           xlab = 'Year', ylab='Expected \n hauled-out proportion') +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  scale_x_continuous(limits=c(1987.5,2023.5), breaks=seq(1990,2024,5)) +
  theme(text = element_text(size=font_size))
ggsave(filename='annual_haulout_prob.png', width=13, height=4.5, path=paste(output_directory, 'plots/', sep=''))

################## Plot expected haul-out probability by ice cover ##############################

C <- seq(1e-9,38,0.01) # setup grid of ice extents
C.std <- (C - mean(data$ice)) / sd(data$ice) # standardized ice extent
X_w <- cbind(rep(1, length(C)), C.std) # design matrix of covariates for haul-out

w_ice <- matrix(NA, nrow(samples$n0), length(C))
for(i in 1:nrow(w_ice)) {
  zeta_0 <- (1-samples$w[i])/samples$w[i]*(samples$d[i]+(1-samples$d[i])*logistic(X_w%*%samples$alpha[i,]))
  zeta_0_pup <- (1-samples$w_pup[i])/samples$w_pup[i]*(samples$d[i]+(1-samples$d[i])*logistic(X_w%*%samples$alpha[i,]))
  zeta_1 <- C^2 / (samples$f[i]^2 + C^2)
  w_ice[i,] <- (1-age.str.mean[1,length(data$years)])*zeta_1/(zeta_1 + zeta_0) + age.str.mean[1,length(data$years)]*zeta_1/(zeta_1 + zeta_0_pup)
}
wt_obs <- apply(t(data$y_survey/t(N.tot.samples[,1:ncol(N.tot.samples)])), 2, median)
w_ice_quantiles <- t(apply(w_ice, 2, quantile, quantiles))

ggplot() +
  geom_ribbon(aes(x=1e3*C, ymin=w_ice_quantiles[,1], ymax=w_ice_quantiles[,3]), alpha=0.15, fill='#CC6677') +
  geom_line(aes(x=1e3*C, y=w_ice_quantiles[,2]), linewidth=0.3, color='#CC6677') +
  geom_point(aes(x=1e3*data$ice, y=wt_obs), shape=20, size=2, alpha=1, color='#882255') +
  geom_path(aes(x=rev(1e3*data$ice), y=seq(max(wt_obs, na.rm=T),0,length.out=(length(data$years)+1))), color='#0072B2', linewidth=2) +
  scale_x_continuous(breaks=seq(0,36e3,6e3), label = function(x){label_k(1e-3*x, 'e3')},
                     name=expression('Ice extent (km'^2*')')) +
  scale_y_continuous(limits = c(0,max(wt_obs, na.rm=T)), breaks=seq(0,1,0.1), name='Expected hauled-out proportion', position='right',
                      sec.axis = sec_axis(~ 1987.5+(. / max(wt_obs, na.rm=T))*((max(data$years)+1)-1987.5), 
                                         breaks=seq(1990,2023,5), name = "Year")) +
  coord_flip() + theme_classic() +
  theme(text = element_text(size=font_size),
        axis.title.x.top = element_text(colour = "#CC6677"),
        axis.text.x.top = element_text(colour = "#CC6677"),
        axis.title.x.bottom = element_text(colour = "#0072B2"),
        axis.text.x.bottom = element_text(colour = "#0072B2"))
ggsave(filename='haulout_ice.png', width=13, height=4.5, path=paste(output_directory, 'plots/', sep=''))


###################################################################################################### 
################################## SAVE PARAMETER ESTIMATES ##########################################
###################################################################################################### 

output.param.names <- c('N_1988', 'N_2023', 
                   'growth_rate_1989', 'growth_rate_2021', 'growth_rate_2023', 'growth_rate_max', 'growth_rate_asymptotic', 'R0_asymptotic',
                   paste('1988_freq', 0:5, sep=''), paste('2023_freq', 0:5, sep=''), paste('K_freq_', c(0,'sub','5+'), sep=''), 'f_freq',
                   paste('phi_', 1:12, sep=''), 'phi_scale', 'nu0', 'nu1', 'c',
                   'b_1988', 'b0_2023', 'b_2023', 'b_min', 'b_max', 'b_scale', 'b_eq', 'beta_0', 'beta_1',
                   'p_1988', 'p_2022', 'p_max', 'p_min', 'abortion_rate_2022', 
                   'K', 'theta_0', 'theta_0_raw',
                   'w', 'w_pup', 'delta', 'w_baseline', 'wt_max', 'w_ice_max', 'w_ice_min', 
                   'd', 'f', 'alpha_0', 'alpha_1',
                   'E_fi', 'E_sw_spring', 'E_sw_fall', 'sigma_fi', 'sigma_sw_spring', 'sigma_sw_fall',
                   'pi_s', 'pi_c', 'kappa', 'r',
                   'rel_bc_0_1', 'rel_bc0_m_f',
                   'max_Rhat', 'min_ESS'
                   )

estimates <- matrix(NA, length(output.param.names), 3)
rownames(estimates) <- output.param.names
colnames(estimates) <- c('2.5', '50', '97.5')

# Population size and growth
estimates['N_1988',] <- N.quantiles[1,]
estimates['N_2023',] <- N.quantiles[nrow(N.quantiles),]
estimates['growth_rate_1989',] <- growth.rate.quantiles[1,]
estimates['growth_rate_2021',] <- growth.rate.quantiles[nrow(growth.rate.quantiles)-2,]
estimates['growth_rate_2023',] <- growth.rate.quantiles[nrow(growth.rate.quantiles),]
estimates['growth_rate_max',] <- apply(growth.rate.quantiles, 2, max)
estimates['growth_rate_asymptotic',] <- g_asymptotic_quantiles
estimates['R0_asymptotic',] <- quantile(0.5*samples$b_max*exp(rowSums(log(samples$phi[,1:5])))/(1-samples$phi[,6]), quantiles)

# Age structure
demographic.str.1988 <- samples$N[,,1] / rowSums(samples$N[,,1])
age.str.1988 <- demographic.str.1988[,1:6] + demographic.str.1988[,7:12]
estimates[9:14,] <- t(apply(age.str.1988, 2, quantile, quantiles))
demographic.str.2023 <- samples$N[,,length(data$years)+1] / rowSums(samples$N[,,length(data$years)+1])
age.str.2023 <- demographic.str.2023[,1:6] + demographic.str.2023[,7:12]
estimates[15:20,] <- t(apply(age.str.2023, 2, quantile, quantiles))
stable.str.K.simplified <- matrix(NA, nrow(stable.str.K.samps), 3)
stable.str.K.simplified[,c(1,3)] <- stable.str.K.samps[,c(1,6)]
stable.str.K.simplified[,2] <- rowSums(stable.str.K.samps[,2:5])
estimates[21:23,] <- t(apply(stable.str.K.simplified, 2, quantile, quantiles))
estimates['f_freq',] <- quantile(rowSums(demographic.str.2023[,1:6]), quantiles)

# Survival rates
estimates[25:36,] <- t(apply(samples$phi, 2, quantile, quantiles))
estimates['phi_scale',] <- quantile(samples$phi_sc, quantiles)
estimates['nu0',] <- quantile(samples$nu[,1], quantiles)
estimates['nu1',] <- quantile(samples$nu[,2], quantiles)
estimates['c',] <- quantile(samples$c, quantiles)

# Haul-out
estimates['w',] <- quantile(samples$w, quantiles)
estimates['w_pup',] <- quantile(samples$w*samples$delta, quantiles)
estimates['delta',] <- quantile(samples$delta, quantiles)
estimates['w_baseline',] <- w_ice_quantiles[which.min(abs(C - median(data$ice))),]
estimates['wt_max',] <- quantile(apply(w_t, 1, max), quantiles)
estimates['w_ice_max',] <- quantile(apply(w_ice, 1, max), quantiles)
estimates['w_ice_min',] <- quantile(apply(w_ice[,(0.2*ncol(w_ice)):ncol(w_ice)], 1, min), quantiles)
estimates['f',] <- quantile(samples$f, quantiles)
estimates['alpha_0',] <- quantile(samples$alpha[,1], quantiles)
estimates['alpha_1',] <- quantile(samples$alpha[,2], quantiles)
estimates['d',] <- quantile(samples$d, quantiles)

# Reproduction
estimates['p_1988',] <- quantile(samples$p[,1], quantiles)
estimates['p_2022',] <- quantile(samples$p[,35], quantiles)
estimates['p_max',] <- quantile(samples$p_max, quantiles)
estimates['p_min',] <- quantile(samples$p_min, quantiles)
estimates['b_1988',] <- b.quantiles[1,]
estimates['b0_2023',] <- quantile(samples$b0[,1+length(data$years)], quantiles)
estimates['b_2023',] <- b.quantiles[1+length(data$years),]
estimates['b_min',] <- quantile(samples$b_max*samples$b_scale, quantiles)
estimates['b_max',] <- quantile(samples$b_max, quantiles)
estimates['b_scale',] <- quantile(samples$b_scale, quantiles)
estimates['b_eq',] <- quantile(b_eq, quantiles)
estimates['beta_0',] <- quantile(samples$beta[,1], quantiles)
estimates['beta_1',] <- quantile(samples$beta[,2], quantiles)
estimates['abortion_rate_2022',] <- quantile(1-samples$b[,36]/samples$p[,35], quantiles)

# Density dependence
estimates['K',] <- quantile(samples$K, quantiles)
estimates['theta_0',] <- quantile(samples$theta_0, quantiles)
estimates['theta_0_raw',] <- quantile(samples$theta_0_raw, quantiles)

# Hunting
estimates['E_fi',] <- quantile(samples$E_tot_fi, quantiles)
estimates['E_sw_spring',] <- quantile(samples$E_tot_sw_spring, quantiles)
estimates['E_sw_fall',] <- quantile(samples$E_tot_sw_fall, quantiles)
estimates['sigma_fi',] <- quantile(samples$sigma_E_fi, quantiles)
estimates['sigma_sw_spring',] <- quantile(samples$sigma_E_sw_spring, quantiles)
estimates['sigma_sw_fall',] <- quantile(samples$sigma_E_sw_fall, quantiles)

# Observations
estimates['pi_s',] <- quantile(samples$pi_s, quantiles)
estimates['pi_c',] <- quantile(samples$pi_c, quantiles)
estimates['kappa',] <- quantile(samples$kappa, quantiles)
estimates['r',] <- quantile(samples$r, quantiles)

# ratio of bycatch probability of pups to 1 yo seals
estimates['rel_bc_0_1',] <- quantile(rowSums(psi.bycatch[,c(1,7)])/rowSums(psi.bycatch[,c(2,8)]), quantiles)
# ratio of bycatch probability of male pups to female pups
estimates['rel_bc0_m_f',] <- quantile(psi.bycatch[,7]/psi.bycatch[,1], quantiles)

# Diagnostics
estimates['max_Rhat',] <- max(summary(model)$summary[,'Rhat'], na.rm = T)
estimates['min_ESS',] <- min(summary(model)$summary[,'n_eff'], na.rm = T)


#### Save output ####

write.csv(estimates, file = paste(output_directory, 'estimates.csv', sep=''))

