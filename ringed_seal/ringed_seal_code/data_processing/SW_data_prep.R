
raw_data_directory <- paste(data_directory, 'source_data/raw_data/Sweden/', sep='')
output_data_directory <- paste(data_directory, 'source_data/', sep='')

# load latest age determination data
sample_ages <- read.csv(paste(raw_data_directory, 'SW_aging.csv', sep=''), header = T)
colnames(sample_ages) <- c('accnr', 'age_det')
sample_ages$age_det <- as.numeric(sample_ages$age_det)

# load seal samples from sweden
samples <- read.csv(paste(raw_data_directory, 'SW_samples_raw.csv', sep=''), header = T)

# merge sample data with age determination data using sample ID
samples.complete <- merge(samples, sample_ages, by = "accnr", all.x = TRUE)

# fill in any missing age in the sample data using age determination data
samples.complete$age[!is.na(samples.complete$age_det)] <- samples.complete$age_det[!is.na(samples.complete$age_det)]

# prepare clean Swedish sample data
samples.sw <- samples.complete[,c('year', 'month', 'date', 'source', 'sex', 'age', 'body_length', 'gravid', 'CA', 'placental_scar')]
samples.sw$year <- as.numeric(samples.sw$year)
samples.sw$month <- as.numeric(samples.sw$month)
samples.sw$age <- as.numeric(samples.sw$age)
samples.sw$body_length <- as.numeric(samples.sw$body_length)
samples.sw$gravid <- as.numeric(samples.sw$gravid)
samples.sw$CA <- as.numeric(samples.sw$CA)
samples.sw$placental_scar <- as.numeric(samples.sw$placental_scar)
colnames(samples.sw) <- c('year', 'month', 'date', 'source', 'sex', 'age', 'body_length', 'pregnant', 'CA', 'placental_scar')

# prepare clean Swedish hunting sample data
hunted.samples.sw <- samples.sw[which(samples.sw$source == 'hunted' & !is.na(samples.sw$year)), 
                                which(!(colnames(samples.sw) %in% c('source', 'pregnant', 'CA', 'placental_scar')))]
hunted.samples.sw$spring <- (hunted.samples.sw$month >= 4 & hunted.samples.sw$month < 8)

# prepare clean Swedish bycatch sample data
bycaught.samples.sw <- samples.sw[which(samples.sw$source == 'bycaught' & !is.na(samples.sw$year)), 
                                  which(!(colnames(samples.sw) %in% c('source', 'body_length', 'pregnant', 'CA', 'placental_scar')))]

# prepare clean Swedish reproduction data
reproduction.data.sw <- samples.sw[which(!is.na(samples.sw$year) & 
                                           !is.na(samples.sw$month) & 
                                           samples.sw$sex=='f' &
                                           (!is.na(samples.sw$pregnant) | !is.na(samples.sw$CA) | !is.na(samples.sw$placental_scar)) &
                                           samples.sw$age >= 4), 
                                   c('year', 'month', 'age', 'CA', 'placental_scar', 'pregnant')]

# save files
write.csv(samples.sw, file = paste(output_data_directory, 'samples_SW.csv', sep=''), row.names = F)
write.csv(hunted.samples.sw, file = paste(output_data_directory, 'hunting_samples_SW.csv', sep=''), row.names = F)
write.csv(bycaught.samples.sw, file = paste(output_data_directory, 'bycatch_samples_SW.csv', sep=''), row.names = F)
write.csv(reproduction.data.sw, file = paste(output_data_directory, 'reproduction_data_SW.csv', sep=''), row.names = F)


