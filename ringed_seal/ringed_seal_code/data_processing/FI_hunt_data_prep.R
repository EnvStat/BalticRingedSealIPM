

############################ Global parameters ###############################

source_directory <- paste(data_directory, 'source_data/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

# years included in the model
years <- 1988:2022

############################ Load and setup raw data ###############################

samples.raw <- read.csv(file = paste(source_directory, 'hunting_samples_FI.csv', sep=''), header = T)

# assign to demographic classes
samples.raw$age_class <- samples.raw$age
samples.raw$age_class[samples.raw$age_class > 5] <- 5
samples.raw$class <- 1 + 6*(samples.raw$sex == 'm') + samples.raw$age_class

# compute time between parturition (~March 1) and hunting
samples.raw$date <- paste(samples.raw$year,'-',samples.raw$month,'-',samples.raw$day, sep='')
samples.raw$date[which(is.na(samples.raw$month) | is.na(samples.raw$day))] <- NA
samples.raw$date <- as.Date(samples.raw$date)

samples.raw$hunting_time <- NA
samples.raw$hunting_time[which(samples.raw$month < 3)] <- 
  samples.raw$date[which(samples.raw$month < 3)]-as.Date(paste(samples.raw$year[which(samples.raw$month < 3)]-1, '-3-1', sep=''))
samples.raw$hunting_time[which(samples.raw$month >= 3)] <- 
  samples.raw$date[which(samples.raw$month >= 3)]-as.Date(paste(samples.raw$year[which(samples.raw$month >= 3)], '-3-1', sep=''))
samples.raw$hunting_time[which(is.na(samples.raw$date) | samples.raw$hunting_time < 0)] <- NA
samples.raw$hunting_time <- samples.raw$hunting_time / 365

# select years included in Finnish hunting model
samples.raw <- samples.raw[which(samples.raw$year >= 2016 & samples.raw$year < 2023),]


########################## Prepare hunting sample data ###############################

# separate data to groups based on available information (exclude 2022, age determination is incomplete)
samples.fi <- samples.raw[which(!is.na(samples.raw$class) & samples.raw$year < 2022), c('year', 'class')]
samples.fi.missing_age <- samples.raw[which(is.na(samples.raw$age) & !is.na(samples.raw$sex) & samples.raw$year < 2022), c('year', 'sex')]
samples.fi.missing_sex <- samples.raw[which(is.na(samples.raw$sex) & !is.na(samples.raw$age) & samples.raw$year < 2022), c('year', 'age_class')]

# hunting samples with both age & sex available
hunting.samples <- matrix(0, length(years), 12)
rownames(hunting.samples) <- years
colnames(hunting.samples) <- 1:12

hunting.incomplete <- table(samples.fi)
hunting.samples[which(rownames(hunting.samples) %in% rownames(hunting.incomplete)), 
               which(colnames(hunting.samples) %in% colnames(hunting.incomplete))] <- hunting.incomplete

# samples with missing sex info
hunting.samples.missing_sex <- matrix(0, length(years), 6)
rownames(hunting.samples.missing_sex) <- years
colnames(hunting.samples.missing_sex) <- 1:6

hunting.incomplete.missing_sex <- table(samples.fi.missing_sex)
hunting.samples.missing_sex[which(rownames(hunting.samples.missing_sex) %in% rownames(hunting.incomplete.missing_sex)), 
                           which(colnames(hunting.samples.missing_sex) %in% colnames(hunting.incomplete.missing_sex))] <- hunting.incomplete.missing_sex

# samples with missing age info
hunting.samples.missing_age <- matrix(0, length(years), 2)
rownames(hunting.samples.missing_age) <- years
colnames(hunting.samples.missing_age) <- c('f', 'm')

# spring samples with missing age info
hunting.incomplete.missing_age <- table(samples.fi.missing_age)
hunting.samples.missing_age[which(rownames(hunting.samples.missing_age) %in% rownames(hunting.incomplete.missing_age)), 
                           which(colnames(hunting.samples.missing_age) %in% colnames(hunting.incomplete.missing_age))] <- hunting.incomplete.missing_age


############################# Prepare data on timing of hunting ####################################

# mean occurrence time of hunting measured from March 1 (census point)
hunting.times.incomplete <- aggregate(hunting_time ~ year, FUN='mean', data=samples.raw)
hunting.times <- matrix(0, length(years), 1)
rownames(hunting.times) <- years

hunting.times[which(years %in% hunting.times.incomplete$year)] <- hunting.times.incomplete$hunting_time

hunting.times.samples <- samples.raw$hunting_time[!is.na(samples.raw$hunting_time)]

############################# Save data #######################################

write.csv(hunting.times, file = paste0(output_directory, 'hunting_times_fi.csv'))
write.csv(hunting.times.samples, file = paste0(source_directory, 'hunting_times_samples_fi.csv'))

# save data
write.csv(hunting.samples, file = paste0(output_directory, 'hunting_samples_fi.csv'))
write.csv(hunting.samples.missing_age, file = paste0(output_directory, 'hunting_samples_fi_ageNA.csv'))
write.csv(hunting.samples.missing_sex, file = paste0(output_directory, 'hunting_samples_fi_sexNA.csv'))



