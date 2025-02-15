
# load packages
library(loo)
library(ggplot2)

font_size <- 25 # set font size for plots

# function for computing loo
compute_k <- function(log_lik) {
  ll <- log_lik[,which(colSums(log_lik) != 0)]
  loocv <- loo(ll)
  pareto_k <- loocv$diagnostic$pareto_k
  looic <- loocv$looic
  waic <- waic(ll)$waic
  return(list(pareto_k=pareto_k, looic=looic, waic=waic))
}

# compute loo-cv for each dataset
survey_loo <- compute_k(samples$log_lik_survey[,!as.logical(data$exclude_survey)])
reproduction_loo <- compute_k(samples$log_lik_reproduction[,!as.logical(data$exclude_rep)])
hb_fi_loo <- compute_k(samples$log_lik_hb_fi[,!as.logical(data$exclude_hb)])
hb_sw_loo <- compute_k(samples$log_lik_hb_sw[,!as.logical(data$exclude_hb)])
hs_fi_loo <- compute_k(samples$log_lik_hs_fi[,!as.logical(data$exclude_samples)])
hs_sw_loo <- compute_k(samples$log_lik_hs_sw_spring[,!as.logical(data$exclude_samples)] +
                       samples$log_lik_hs_sw_fall[,!as.logical(data$exclude_samples)])
bc_loo <- compute_k(samples$log_lik_bc[,!as.logical(data$exclude_samples)])

# aggregate outputs
labs <- c("Aerial survey", "Reproduction", "Hunting bag (SW)", "Hunting bag (FI)",
          "Hunting samples (SW)", "Hunting samples (FI)", "Bycatch samples")
k_mean <- c(mean(survey_loo$pareto_k), mean(reproduction_loo$pareto_k), mean(hb_sw_loo$pareto_k), 
            mean(hb_fi_loo$pareto_k), mean(hs_sw_loo$pareto_k), mean(hs_fi_loo$pareto_k), 
            mean(bc_loo$pareto_k))
k_sum <- c(sum(survey_loo$pareto_k), sum(reproduction_loo$pareto_k), sum(hb_sw_loo$pareto_k), 
           sum(hb_fi_loo$pareto_k), sum(hs_sw_loo$pareto_k), sum(hs_fi_loo$pareto_k), 
           sum(bc_loo$pareto_k))
looic <- c(sum(survey_loo$looic), sum(reproduction_loo$looic), sum(hb_sw_loo$looic), 
           sum(hb_fi_loo$looic), sum(hs_sw_loo$looic), sum(hs_fi_loo$looic), 
           sum(bc_loo$looic))
waic <- c(sum(survey_loo$waic), sum(reproduction_loo$waic), sum(hb_sw_loo$waic), 
           sum(hb_fi_loo$waic), sum(hs_sw_loo$waic), sum(hs_fi_loo$waic), 
           sum(bc_loo$waic))

# create plots
mean.df <- data.frame(Data = factor(labs, levels=labs), Effect = k_mean)
ggplot(data=mean.df, aes(x=Data, y=Effect)) +
  geom_bar(stat="identity") +
  labs(x = 'Data type', y = 'Mean Pareto-k') + 
  theme_classic() +
  theme(text = element_text(size=font_size), 
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename='pareto_k_mean.png', width=6, height=10, path=paste0(output_directory, 'plots/'))

sum.df <- data.frame(Data = factor(labs, levels=labs), Effect = k_sum)
ggplot(data=sum.df, aes(x=Data, y=Effect)) +
  geom_bar(stat="identity") +
  labs(x = 'Data type', y = 'Total Pareto-k') +
  theme_classic() +
  theme(text = element_text(size=font_size), 
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename='pareto_k_sum.png', width=6, height=10, path=paste0(output_directory, 'plots/'))


# output loo-cv results
if(compute_loo) {
  loo.df <- data.frame(data = labs, looic = looic, waic = waic)
  write.csv(loo.df, file = paste0(output_directory, 'loo.csv')) 
}

