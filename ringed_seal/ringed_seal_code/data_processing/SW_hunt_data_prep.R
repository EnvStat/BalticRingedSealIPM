

############################ Global parameters ###############################

source_directory <- paste(data_directory, 'source_data/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

# years included in the model
years <- 1988:2022

############################ Load and setup raw data ###############################

samples.raw <- read.csv(file = paste(source_directory, 'hunting_samples_SW.csv', sep=''), header = T)

# assign to demographic classes
samples.raw$age_class <- samples.raw$age
samples.raw$age_class[samples.raw$age_class > 5] <- 5
samples.raw$class <- 1 + 6*(samples.raw$sex == 'm') + samples.raw$age_class

# compute time between parturition (~March 1) and hunting
samples.raw$hunting_time <- NA
samples.raw$hunting_time[which(samples.raw$month < 3)] <- 
  as.Date(samples.raw$date[which(samples.raw$month < 3)])-
  as.Date(paste(samples.raw$year[which(samples.raw$month < 3)]-1, '-3-1', sep=''))
samples.raw$hunting_time[which(samples.raw$month >= 3)] <- 
  as.Date(samples.raw$date[which(samples.raw$month >= 3)])-
  as.Date(paste(samples.raw$year[which(samples.raw$month >= 3)], '-3-1', sep=''))
samples.raw$hunting_time[which(is.na(samples.raw$date) | samples.raw$hunting_time < 0)] <- NA
samples.raw$hunting_time <- samples.raw$hunting_time / 365

# select years included in Swedish hunting model
samples.raw <- samples.raw[which(samples.raw$year >= 2015),]


########################## Prepare hunting sample data ###############################

# separate data to groups based on available information (exclude 2022, age determination is incomplete)
samples.sw <- samples.raw[which(!is.na(samples.raw$class) & samples.raw$year < 2022), 
                          c('year', 'class', 'spring')]
samples.sw.missing_age <- samples.raw[which(is.na(samples.raw$age) & !is.na(samples.raw$sex) & samples.raw$year < 2022), 
                                      c('year', 'sex', 'spring')]
samples.sw.missing_sex <- samples.raw[which(is.na(samples.raw$sex) & !is.na(samples.raw$age) & samples.raw$year < 2022), 
                                      c('year', 'age_class', 'spring')]

# setup data structure for samples
spring.samples <- fall.samples <- matrix(0, length(years), 12)
rownames(spring.samples) <- rownames(fall.samples) <- years
colnames(spring.samples) <- colnames(fall.samples) <- 1:12

# spring samples
spring.incomplete <- table(samples.sw[which(samples.sw$spring==T), c('year', 'class')])
spring.samples[which(rownames(spring.samples) %in% rownames(spring.incomplete)), 
               which(colnames(spring.samples) %in% colnames(spring.incomplete))] <- spring.incomplete

# fall samples
fall.incomplete <- table(samples.sw[which(samples.sw$spring==F), c('year', 'class')])
fall.samples[which(rownames(fall.samples) %in% rownames(fall.incomplete)), 
             which(colnames(fall.samples) %in% colnames(fall.incomplete))] <- fall.incomplete

# setup data structure for samples with missing sex info
spring.samples.missing_sex <- fall.samples.missing_sex <- matrix(0, length(years), 6)
rownames(spring.samples.missing_sex) <- rownames(fall.samples.missing_sex) <- years
colnames(spring.samples.missing_sex) <- colnames(fall.samples.missing_sex) <- 1:6

# spring samples with missing sex info
spring.incomplete.missing_sex <- table(samples.sw.missing_sex[which(samples.sw.missing_sex$spring==1), c('year', 'age_class')])
spring.samples.missing_sex[which(rownames(spring.samples.missing_sex) %in% rownames(spring.incomplete.missing_sex)), 
                           which(colnames(spring.samples.missing_sex) %in% colnames(spring.incomplete.missing_sex))] <- spring.incomplete.missing_sex

# fall samples with missing sex info
fall.incomplete.missing_sex <- table(samples.sw.missing_sex[which(samples.sw.missing_sex$spring==0),
                                                            c('year', 'age_class')])
fall.samples.missing_sex[which(rownames(fall.samples.missing_sex) %in% rownames(fall.incomplete.missing_sex)), 
                         which(colnames(fall.samples.missing_sex) %in% colnames(fall.incomplete.missing_sex))] <- fall.incomplete.missing_sex

# setup data structure for samples with missing age info
spring.samples.missing_age <- fall.samples.missing_age <- matrix(0, length(years), 2)
rownames(spring.samples.missing_age) <- rownames(fall.samples.missing_age) <- years
colnames(spring.samples.missing_age) <- colnames(fall.samples.missing_age) <- c('f', 'm')

# spring samples with missing age info
spring.incomplete.missing_age <- table(samples.sw.missing_age[which(samples.sw.missing_age$spring==T), c('year', 'sex')])
spring.samples.missing_age[which(rownames(spring.samples.missing_age) %in% rownames(spring.incomplete.missing_age)), 
                           which(colnames(spring.samples.missing_age) %in% colnames(spring.incomplete.missing_age))] <- spring.incomplete.missing_age

# fall samples with missing age info
fall.incomplete.missing_age <- table(samples.sw.missing_age[which(samples.sw.missing_age$spring==F), c('year', 'sex')])
fall.samples.missing_age[which(rownames(fall.samples.missing_age) %in% rownames(fall.incomplete.missing_age)), 
                         which(colnames(fall.samples.missing_age) %in% colnames(fall.incomplete.missing_age))] <- fall.incomplete.missing_age


############################# Prepare data on timing of hunting ####################################

# mean occurrence time of hunting in spring and summer measured from March 1 (census point)
hunting.times.spring.incomplete <- aggregate(hunting_time ~ year, FUN='mean', data=samples.raw[which(samples.raw$spring==1),])
hunting.times.fall.incomplete <- aggregate(hunting_time ~ year, FUN='mean', data=samples.raw[which(samples.raw$spring==0),])

hunting.times.spring <- hunting.times.fall <- matrix(0, length(years), 1)
rownames(hunting.times.spring) <- rownames(hunting.times.fall) <- years

hunting.times.spring[which(years %in% hunting.times.spring.incomplete$year)] <- hunting.times.spring.incomplete$hunting_time
hunting.times.fall[which(years %in% hunting.times.fall.incomplete$year)] <- hunting.times.fall.incomplete$hunting_time

hunting.times <- cbind(hunting.times.spring, hunting.times.fall)
colnames(hunting.times) <- c('spring', 'fall')

hunting.times.samples <- samples.raw$hunting_time[!is.na(samples.raw$hunting_time)]

############################# Save data #######################################

write.csv(hunting.times, file = paste0(output_directory, 'hunting_times_sw.csv'))
write.csv(hunting.times.samples, file = paste0(source_directory, 'hunting_times_samples_sw.csv'))

# save spring data
write.csv(spring.samples, file = paste0(output_directory, 'hunting_samples_sw_spring.csv'))
write.csv(spring.samples.missing_age, file = paste0(output_directory, 'hunting_samples_sw_spring_ageNA.csv'))
write.csv(spring.samples.missing_sex, file = paste0(output_directory, 'hunting_samples_sw_spring_sexNA.csv'))

# save fall data
write.csv(fall.samples, file = paste(output_directory, 'hunting_samples_sw_fall.csv'))
write.csv(fall.samples.missing_age, file = paste0(output_directory, 'hunting_samples_sw_fall_ageNA.csv'))
write.csv(fall.samples.missing_sex, file = paste0(output_directory, 'hunting_samples_sw_fall_sexNA.csv'))

