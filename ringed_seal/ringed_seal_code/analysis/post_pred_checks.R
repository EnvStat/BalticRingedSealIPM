

############### Aerial survey estimates ###################

label_k = function(x, s) {
  labels = paste(x, s, sep='')
  labels[x==0] <- 0
  return(labels)
}
data$y_survey[data$y_survey==0] <- NA
pp_raster(samples$y_survey_sim, c(data$years, data$years[length(data$years)]+1), 200, 
          c(data$years, data$years[length(data$years)]+1), data$y_survey, as.factor(data$exclude_survey), cex=2,
          x_lab='Year', y_lab='Aerial survey estimate', interpolate = F) + 
  scale_y_continuous(limits = c(0, 28e3), breaks=seq(0,25e3,5e3), label = function(x){label_k(1e-3*x, 'e3')}) +
  scale_x_continuous(breaks=seq(1990,2020,5)) +
  guides(fill=guide_colorbar(title.position='right', direction='horizontal'),
         colour=guide_legend(override.aes = list(size=3))) +
  theme(legend.position = c(0.3,0.7), legend.spacing.y = unit(0.05, "cm"), 
        legend.background = element_rect(fill = NA), panel.grid.minor = element_blank(), text = element_text(size=font_size))

# save plot
ggsave(filename='survey_pp.png', width=13, height=4.5, path=paste(output_directory,'plots/', sep=''))


################ Pregnancy data #################

data$y_pregnant[data$pregnant_ss==0] <- NA

# posterior quantiles for true pregnancy rate
p_preg.quantiles <- apply(t(t(samples$y_pregnant_sim) / data$pregnant_ss), 2, quantile, quantiles, na.rm=T)
p.quantiles <- apply(samples$p[,1:length(data$years)], 2, quantile, quantiles)

ggplot() +
  geom_errorbar(aes(x=data$years, ymin=p_preg.quantiles[1,], ymax=p_preg.quantiles[3,]), alpha=0.35) +
  geom_linerange(aes(x=data$years, ymin=p_preg.quantiles[1,], ymax=p_preg.quantiles[3,], color = "gray16"), alpha=0) +
  geom_ribbon(aes(x=data$years, ymin=p.quantiles[1,], ymax=p.quantiles[3,], fill='gray16'), alpha=0.15) +
  geom_line(aes(x=data$years, y=p.quantiles[2,]), col='gray16', alpha=0.5) +
  geom_point(aes(x=data$years, y=data$y_pregnant / data$pregnant_ss, color='black', size=data$pregnant_ss, group=as.factor(data$exclude_rep), shape=as.factor(data$exclude_rep))) +
  scale_color_identity(name='', guide='legend', labels=c('Observed', '95% Posterior \n predictive interval')) +
  scale_fill_identity(name='', guide='legend', labels=c('Posterior estimate \n of true value')) +
  scale_size_continuous(name='Sample size', guide='legend', 
                      breaks=round(seq(min(data$pregnant_ss[data$pregnant_ss>0]), max(data$pregnant_ss), length.out=4))) +
  scale_shape_manual(values=c(16, 1), guide='none') +
  guides(color=guide_legend(override.aes = list(alpha=c(1, 0.5), linetype=c('blank', 'solid'), shape=c(19, NA))), 
         fill=guide_legend(override.aes = list(alpha=0.25))) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x='Year', y='Pregnancy rate') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=font_size), legend.position = 'bottom') +
  guides(size=guide_legend(title.position="top", label.position = 'bottom'),
         fill='none', color='none')

# save plot
ggsave(filename='preg_rate_pp.png', width=5, height=10, path=paste(output_directory, 'plots/', sep=''))

################ Corpus albicans data #################

CA_ss <- data$CA_ss + colSums(data$z) # sample size for CA
n_CA_obs <- data$y_CA + colSums(data$z[c(2,4,5),]) # number of observed CA
n_CA_sim <- samples$y_CA_sim + apply(samples$z_sim[,c(2,4,5),], c(1,3), sum) # simulated number of CA
CA_prop_sim <- t(t(n_CA_sim) / CA_ss) # simulated % of samples with CA

# posterior quantiles for expected % of samples with CA
CA_rate <- apply(samples$gamma[,c(2,4,5),], c(1,3), sum)
CA_rate.quantiles <- apply(CA_rate, 2, quantile, quantiles, na.rm=T)
CA_prop_sim.quantiles <- apply(CA_prop_sim, 2, quantile, quantiles, na.rm=T)

# plot posterior predictive check
ggplot() +
  geom_errorbar(aes(x=data$years, ymin=CA_prop_sim.quantiles[1,], ymax=CA_prop_sim.quantiles[3,]), alpha=0.35) +
  geom_ribbon(aes(x=data$years, ymin=CA_rate.quantiles[1,], ymax=CA_rate.quantiles[3,]), alpha=0.15) +
  geom_line(aes(x=data$years, y=CA_rate.quantiles[2,]), alpha=0.5) +
  geom_point(aes(x=data$years, y=n_CA_obs / CA_ss, size=CA_ss, group=as.factor(data$exclude_rep), shape=as.factor(data$exclude_rep)), color='black') +
  scale_size_continuous(name='Sample size', guide='legend', 
                        breaks=round(seq(min(CA_ss[CA_ss>0]), max(CA_ss), length.out=4))) +
  scale_shape_manual(values=c(16, 1), guide='none') +
  scale_y_continuous(limits=c(0,1)) +
  labs(x='Year', y='Proportion with CA') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=font_size), legend.position = 'bottom') +
  guides(size=guide_legend(title.position="top", label.position = 'bottom'))

# save plot
ggsave(filename='CA_prop_pp.png', width=5, height=10, path=paste(output_directory,'plots/', sep=''))

################ Placental scar data #################

scar_ss <- data$scar_ss + colSums(data$z[1:4,]) # sample size for placental scars
n_scars_obs <- data$y_scar + colSums(data$z[c(3,4),]) # observed number of scars
n_scars_sim <- samples$y_scar_sim + apply(samples$z_non_na_sim[,c(3,4),], c(1,3), sum) # simulated number of scars
scar_prop_sim <- t(t(n_scars_sim) / scar_ss) # simulated scar %

# posterior quantiles for expected % of samples with placental scars
scar_rate <- apply(samples$gamma[,c(3,4),], c(1,3), sum) / apply(samples$gamma[,1:4,], c(1,3), sum)
scar_rate.quantiles <- apply(scar_rate, 2, quantile, quantiles, na.rm=T)
scar_prop_sim.quantiles <- apply(scar_prop_sim, 2, quantile, quantiles, na.rm=T)

# plot posterior predictive check
ggplot() +
  geom_errorbar(aes(x=data$years, ymin=scar_prop_sim.quantiles[1,], ymax=scar_prop_sim.quantiles[3,]), alpha=0.35) +
  geom_ribbon(aes(x=data$years, ymin=scar_rate.quantiles[1,], ymax=scar_rate.quantiles[3,]), alpha=0.15) +
  geom_line(aes(x=data$years, y=scar_rate.quantiles[2,]), alpha=0.5) +
  geom_point(aes(x=data$years, y=n_scars_obs / scar_ss, size=scar_ss, group=as.factor(data$exclude_rep), shape=as.factor(data$exclude_rep)), color='black') +
  scale_size_continuous(name='Sample size', guide='legend', 
                      breaks=round(seq(min(scar_ss[scar_ss>0]), max(scar_ss), length.out=4))) +
  scale_shape_manual(values=c(16, 1), guide='none') +
  scale_y_continuous(limits=c(0,1)) +
  labs(x='Year', y='Proportion with placental scar') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=font_size), legend.position = 'bottom') +
  guides(size=guide_legend(title.position="top", label.position = 'bottom'))

# save plot
ggsave(filename='scar_prop_pp.png', width=5, height=10, path=paste(output_directory, 'plots/', sep=''))


############### Size of the Finnish hunting bag ######################

pp_raster(samples$y_hb_sw_sim, data$years, 3, data$years, data$y_hb_sw, as.factor(data$exclude_hb), cex=2,
          x_lab='Year', y_lab='Total harvest in Sweden', interpolate = F) + 
  scale_x_continuous(limits = c(2014.5, 2022.5), breaks=2015:2022) +
  guides(fill=guide_colorbar(title.position='top', direction='horizontal'),
         colour=guide_legend(override.aes = list(size=3))) +
  theme(legend.position = c(0.35,0.85), legend.spacing.y = unit(0.05, "cm"), 
        legend.background = element_rect(fill = NA), panel.grid.minor = element_blank(), text = element_text(size=font_size))

# save plot
ggsave(filename='hbag_sw_pp.png', width=8, height=8, path=paste(output_directory, 'plots/', sep=''))


############### Size of the Finnish hunting bag ######################

pp_raster(samples$y_hb_fi_sim, data$years, 3, data$years, data$y_hb_fi, as.factor(data$exclude_hb), cex=2,
          x_lab='Year', y_lab='Total harvest in Finland', interpolate = F) + 
  scale_x_continuous(limits = c(2015.5, 2022.5), breaks=2016:2022) +
  theme(panel.grid.minor = element_blank(), legend.position = "none", text = element_text(size=font_size))

# save plot
ggsave(filename='hbag_fi_pp.png', width=8, height=8, path=paste(output_directory, 'plots/', sep=''))


################ Finnish hunting samples ###############

# indexes for years with available hunting sample data
if(sum(data$exclude_samples) == 0) {
  idx <- which(colSums(data$y_hs_fi) > 0)
} else { # out-of-sample predictions for data sensitivities
  idx <- which(colSums(data$y_hs_fi) > 0 & as.logical(data$exclude_samples))
}
  
# simulated samples from the Finnish hunting bag
samples_fi <- samples$y_hs_fi_sim[floor(seq(1,nrow(samples$n0),length.out=100)),,idx]

# plot posterior predictive check
df.samples_fi <- data.frame(n.samps = c(samples_fi),
                            draw = rep(1:dim(samples_fi)[1], dim(samples_fi)[2]*dim(samples_fi)[3]))
ggplot() +
  stat_density(data=df.samples_fi, aes(n.samps, group=draw, color='#0D65B3'), 
               geom="line", position="identity", linewidth=0.035) +
  stat_density(aes(c(data$y_hs_fi[,idx]), color='black'), 
               geom="line", position="identity", linewidth=1) +
  scale_color_identity(name='', guide='legend', labels=c('Posterior samples', 'True distribution')) +
  guides(color=guide_legend(override.aes = list(linewidth=3))) +
  labs(x='Number of samples', y='Density') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=font_size), legend.position = 'none')

# save plot
ggsave(filename='h_samples_fi_pp.png', width=7, height=10, path=paste(output_directory, 'plots/', sep=''))


################ Swedish hunting samples ###############

# indexes for years with available hunting sample data
if(sum(data$exclude_samples) == 0) {
  idx <- which(colSums(data$y_hs_sw_spring+data$y_hs_sw_fall) > 0)
} else { # out-of-sample predictions for data sensitivities
  idx <- which(colSums(data$y_hs_sw_spring+data$y_hs_sw_fall) > 0 & as.logical(data$exclude_samples))
}

samples_sw_spring <- samples$y_hs_sw_spring_sim[floor(seq(1,nrow(samples$n0),length.out=100)),,idx]
samples_sw_fall <- samples$y_hs_sw_fall_sim[floor(seq(1,nrow(samples$n0),length.out=100)),,idx]
samples_sw <- samples_sw_fall + samples_sw_spring

# plot posterior predictive check
df.samples_sw <- data.frame(n.samps = c(samples_sw),
                            draw = rep(1:nrow(samples_sw), ncol(samples_sw)))
ggplot() +
  geom_density(data=df.samples_sw, aes(n.samps, group=draw), linewidth=0.035, col='#0D65B3') +
  geom_density(aes(c(data$y_hs_sw_spring[,idx] +
                     data$y_hs_sw_fall[,idx])), 
               linewidth=1) +
  labs(x='Number of samples', y='Density') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), legend.position = "none", text = element_text(size=font_size))

# save plot
ggsave(filename='h_samples_sw_pp.png', width=7, height=10, path=paste(output_directory, 'plots/', sep=''))


################ Bycatch samples ###############

# indexes for years with available bycatch sample data
if(sum(data$exclude_samples) == 0) {
  idx <- which(colSums(data$y_bc) > 0)
} else { # out-of-sample predictions for data sensitivities
  idx <- which(colSums(data$y_bc) > 0 & as.logical(data$exclude_samples))
}

# simulated bycatch samples
samples_bycatch <- samples$y_bc_sim[floor(seq(1,nrow(samples$n0),length.out=100)),,idx]

# plot posterior predictive check
df.samples_bycatch <- data.frame(n.samps = c(samples_bycatch),
                            draw = rep(1:dim(samples_bycatch)[1], dim(samples_bycatch)[2]*dim(samples_bycatch)[3]))
ggplot() +
  stat_density(data=df.samples_bycatch, aes(n.samps, group=draw, col='#0D65B3'), 
               geom="line",position="identity", bw=0.75, linewidth=0.035) +
  stat_density(aes(c(data$y_bc[,idx]), col='black'), geom="line",position="identity", bw=0.75, linewidth=1) +
  scale_y_continuous(trans = 'sqrt') +
  scale_color_identity(name='', guide='legend', labels=c('Simulated', 'Observed')) +
  labs(x='Number of samples', y='Density') +
  theme_classic() +
  guides(color=guide_legend(override.aes = list(linewidth=2))) +
  theme(legend.position = c(0.7,0.8), legend.key.width = unit(3, "line"),
        panel.grid.minor = element_blank(), text = element_text(size=font_size))

# save plot
ggsave(filename='bycatch_samples_pp.png', width=7, height=10, path=paste(output_directory, 'plots/', sep=''))

