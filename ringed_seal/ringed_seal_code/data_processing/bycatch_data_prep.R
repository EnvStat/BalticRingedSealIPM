

############################ Global parameters ###############################

source_directory <- paste(data_directory, 'source_data/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

# years included in the model
years <- 1988:2022

############################ Load and setup raw data ###############################

samples.fi.raw <- read.csv(file = paste(source_directory, 'bycatch_samples_FI.csv', sep=''), header = T)
samples.sw.raw <- read.csv(file = paste(source_directory, 'bycatch_samples_SW.csv', sep=''), header = T)
samples.raw <- rbind(samples.fi.raw[,c('year','sex','age')], samples.sw.raw[,c('year','sex','age')])

# assign to demographic classes
samples.raw$age_class <- samples.raw$age
samples.raw$age_class[samples.raw$age_class > 5] <- 5
samples.raw$class <- 1 + 6*(samples.raw$sex == 'm') + samples.raw$age_class

# select years included in the model (exclude 2022, age determination is incomplete)
samples.raw <- samples.raw[which(samples.raw$year >= 1988 & samples.raw$year < 2022),]


########################## Prepare bycatch sample data ###############################

# separate data to groups based on available information
samples.fi <- samples.raw[which(!is.na(samples.raw$class)), c('year', 'class')]
samples.fi.missing_age <- samples.raw[which(is.na(samples.raw$age) & !is.na(samples.raw$sex)), c('year', 'sex')]
samples.fi.missing_sex <- samples.raw[which(is.na(samples.raw$sex) & !is.na(samples.raw$age)), c('year', 'age_class')]

# bycatch samples with both age & sex available
bycatch.samples <- matrix(0, length(years), 12)
rownames(bycatch.samples) <- years
colnames(bycatch.samples) <- 1:12

bycatch.incomplete <- table(samples.fi)
bycatch.samples[which(rownames(bycatch.samples) %in% rownames(bycatch.incomplete)), 
               which(colnames(bycatch.samples) %in% colnames(bycatch.incomplete))] <- bycatch.incomplete

# samples with missing sex info
bycatch.samples.missing_sex <- matrix(0, length(years), 6)
rownames(bycatch.samples.missing_sex) <- years
colnames(bycatch.samples.missing_sex) <- 1:6

bycatch.incomplete.missing_sex <- table(samples.fi.missing_sex)
bycatch.samples.missing_sex[which(rownames(bycatch.samples.missing_sex) %in% rownames(bycatch.incomplete.missing_sex)), 
                           which(colnames(bycatch.samples.missing_sex) %in% colnames(bycatch.incomplete.missing_sex))] <- bycatch.incomplete.missing_sex

# samples with missing age info
bycatch.samples.missing_age <- matrix(0, length(years), 2)
rownames(bycatch.samples.missing_age) <- years
colnames(bycatch.samples.missing_age) <- c('f', 'm')

# spring samples with missing age info
bycatch.incomplete.missing_age <- table(samples.fi.missing_age)
bycatch.samples.missing_age[which(rownames(bycatch.samples.missing_age) %in% rownames(bycatch.incomplete.missing_age)), 
                           which(colnames(bycatch.samples.missing_age) %in% colnames(bycatch.incomplete.missing_age))] <- bycatch.incomplete.missing_age


############################# Save data #######################################

write.csv(bycatch.samples, file = paste(output_directory, 'bycatch_samples.csv', sep=''))
write.csv(bycatch.samples.missing_age, file = paste(output_directory, 'bycatch_samples_ageNA.csv', sep=''))
write.csv(bycatch.samples.missing_sex, file = paste(output_directory, 'bycatch_samples_sexNA.csv', sep=''))



