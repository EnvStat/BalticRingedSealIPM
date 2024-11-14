

############################ Global parameters ###############################

source_directory <- paste(data_directory, 'source_data/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

# years included in the model
years <- 1988:2022

############################ Load data ###############################

hunting.records <- read.csv(paste(source_directory, 'raw_data/Sweden/SW_hunting_records_clean.csv', sep=''), header=T)
sample.data.fi <- read.csv(paste(source_directory, 'hunting_samples_FI.csv', sep=''), header=T)
sample.data.sw <- read.csv(paste(source_directory, 'hunting_samples_SW.csv', sep=''), header=T)

# exclude hunting records without sampling or body length information
hunting.records <- hunting.records[which(!is.na(hunting.records$sampled) & !is.na(hunting.records$body_length)), 
                                   c('year', 'sex', 'body_length', 'sampled')]
hunting.records$over_100 <- as.numeric(hunting.records$body_length > 100)

# separate hunting records from 2019 and 2020-2021 due to different sampling protocols
hunting.records.19 <- hunting.records[which(hunting.records$year==2019), c('sampled', 'over_100')]
hunting.records.20_21 <- hunting.records[which(hunting.records$year %in% 2020:2021 &
                                                 hunting.records$sex == 'f'), c('sampled', 'over_100')]

# summarize sampling by body length
sampling_summary_19 <- merge(aggregate(sampled ~ over_100, FUN='sum', data=hunting.records.19), 
                          aggregate(sampled ~ over_100, FUN='length', data=hunting.records.19), by='over_100')
colnames(sampling_summary_19) <- c('over_100', 'sampled', 'total')

sampling_summary_20_21 <- merge(aggregate(sampled ~ over_100, FUN='sum', data=hunting.records.20_21), 
                             aggregate(sampled ~ over_100, FUN='length', data=hunting.records.20_21), by='over_100')
colnames(sampling_summary_20_21) <- c('over_100', 'sampled', 'total')

# filter Swedish samples for spring hunting
sample.data.sw <- sample.data.sw[which(sample.data.sw$spring==T & 
                                         !is.na(sample.data.sw$sex) & 
                                         !is.na(sample.data.sw$age) & 
                                         !is.na(sample.data.sw$body_length)), c('sex', 'age', 'body_length')]

# filter Finnish samples for spring hunting
sample.data.fi <- sample.data.fi[which(sample.data.fi$month > 4 & 
                                         sample.data.fi$month < 8 &
                                         !is.na(sample.data.fi$sex) & 
                                         !is.na(sample.data.fi$age) & 
                                         !is.na(sample.data.fi$body_length)), c('sex', 'age', 'body_length')]

# combine SW and FI samples
sample.data <- sample.data.fi #rbind(sample.data.fi, sample.data.sw)

# determine demographic class (1-12) for each seal
sample.data$age_class <- sample.data$age
sample.data$age_class[which(sample.data$age_class > 5)] <- 5
sample.data$class <- 1 + sample.data$age_class + 6*(sample.data$sex == 'm')

# determine if seal meets sampling criteria
hunting.records$over_100 <- as.numeric(hunting.records$body_length > 100)
sample.data$over_100 <- as.numeric(sample.data$body_length > 100)

# simplify sample data
sample.data <- sample.data[, c('class', 'over_100')]

# summarize body length by demographic class
length_summary <- merge(aggregate(over_100 ~ class, FUN='sum', data=sample.data), 
                          aggregate(over_100 ~ class, FUN='length', data=sample.data), by='class')
colnames(length_summary) <- c('class', 'over_100', 'total')


# Save outputs
write.csv(sampling_summary_19, file = paste(output_directory, 'sampling_by_length_16_19.csv', sep=''), row.names = F)
write.csv(sampling_summary_20_21, file = paste(output_directory, 'sampling_by_length_20_21.csv', sep=''), row.names = F)
write.csv(length_summary, file = paste(output_directory, 'length_by_class.csv', sep=''), row.names = F)


