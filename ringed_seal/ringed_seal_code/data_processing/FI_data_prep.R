

raw_data_directory <- paste(data_directory, 'source_data/raw_data/Finland/', sep='')
output_data_directory <- paste(data_directory, 'source_data/', sep='')

# load seal samples from Finland
samples <- read.csv(paste(raw_data_directory, 'FI_samples_raw_clean.csv', sep=''), header = T)

# prepare clean Finnish sample data
samples.fi <- samples[,c('year', 'month', 'day', 'source', 'sex', 'age', 'body_length', 'CA', 'placental_scar')]
samples.fi$sex[which(samples.fi$sex == 1)] <- 'm'
samples.fi$sex[which(samples.fi$sex == 2)] <- 'f'
samples.fi$sex[which(samples.fi$sex == 0)] <- NA

# prepare clean Finnish hunting sample data (source = 1,2,3,5 are hunting)
hunted.samples.fi <- samples.fi[which(samples.fi$source %in% c(1:3,5)), c('year', 'month', 'day', 'sex', 'age', 'body_length')]

# prepare clean Finnish bycatch sample data (source = 4 is bycatch)
bycaught.samples.fi <- samples.fi[which(samples.fi$source == 4), c('year', 'month', 'day', 'sex', 'age')]

# prepare clean Finnish reproduction data
reproduction.data.fi <- samples.fi[which(!is.na(samples.fi$month) & 
                                           samples.fi$sex=='f' &
                                           (!is.na(samples.fi$CA) | !is.na(samples.fi$placental_scar)) &
                                           samples.fi$age >= 4), 
                                   c('year', 'month', 'age', 'CA', 'placental_scar')]

# save files
write.csv(samples.fi, file = paste(output_data_directory, 'samples_FI.csv', sep=''), row.names = F)
write.csv(hunted.samples.fi, file = paste(output_data_directory, 'hunting_samples_FI.csv', sep=''), row.names = F)
write.csv(bycaught.samples.fi, file = paste(output_data_directory, 'bycatch_samples_FI.csv', sep=''), row.names = F)
write.csv(reproduction.data.fi, file = paste(output_data_directory, 'reproduction_data_FI.csv', sep=''), row.names = F)


