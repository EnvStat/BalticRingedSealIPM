

thin <- function(x) {
  x[(c(1:length(x)) %% 2) == 1] <- 1
  return(x)
}

data_1 <- data_2 <- data_3 <- data_4 <- data_5 <- data_6  <- data_7 <- data

# thinned aerial survey data
data_1$exclude_survey <- thin(data_1$exclude_survey)

# only < 2012 survey data
data_2$exclude_survey[which(c(1988:2023) >= 2012)] <- 1

# only >= 2012 survey data
data_3$exclude_survey[which(c(1988:2023) < 2012)] <- 1

# thinned reproduction data
data_4$exclude_rep <- thin(data_4$exclude_rep)

# only >= 2007 reproduction data
data_5$exclude_rep[which(data$years < 2007)] <- 1

# thinned hunting bag data
data_6$exclude_hb <- thin(data_6$exclude_hb)

# thinned hunting & bycatch sample data
data_7$exclude_samples <- thin(data_7$exclude_samples)

# gather all datasets in a list
sensitivity_data <- list(data_1, data_2, data_3, data_4, data_5, data_6, data_7)

