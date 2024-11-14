

############################ Global parameters ###############################

source_directory <- paste(data_directory, 'source_data/', sep='')
output_directory <- paste(data_directory, 'model_data/', sep='')

# years included in the model
years <- 1988:2022


############################ Load data ###############################

samples.sw.raw <- read.csv(file = paste(source_directory, 'reproduction_data_SW.csv', sep=''), header = T)
samples.fi.raw <- read.csv(file = paste(source_directory, 'reproduction_data_FI.csv', sep=''), header = T)

# sample processing for 2022 is incomplete so exclude those
samples.sw.raw <- samples.sw.raw[which(samples.sw.raw$year < 2022),]
samples.fi.raw <- samples.fi.raw[which(samples.fi.raw$year < 2022),]

###################### Prepare pregnancy data ###################

pregnancy.data <- samples.sw.raw[which(samples.sw.raw$month %in% c(8:12) & !is.na(samples.sw.raw$pregnant)), 
                                 c('year', 'pregnant')]

pregnancies.incomplete <- merge(aggregate(pregnant ~ year, FUN=sum, data=pregnancy.data),
                                aggregate(pregnant ~ year, FUN=length, data=pregnancy.data), by='year')

pregnancies <- matrix(0, nrow=length(years), ncol=3)
colnames(pregnancies) <- colnames(pregnancies.incomplete) <- c('year', 'pregnant', 'sample_size')
pregnancies[which(years %in% pregnancies.incomplete$year),] <- as.matrix(pregnancies.incomplete)
pregnancies[,'year'] <- years


###################### Prepare post-partum signs data ###################

# combine Finnish and Swedish samples for placental scar and CA
postpartum.data <- rbind(samples.fi.raw, samples.sw.raw[,-ncol(samples.sw.raw)])
postpartum.data <- postpartum.data[which(postpartum.data$age >= 5),] # remove 4 year olds

# check for fading in placental scars using data from last 10 years
scars.april <- postpartum.data$placental_scar[which(!is.na(postpartum.data$placental_scar) &
                                                      postpartum.data$year > 2012 & 
                                                      postpartum.data$month == 4)]
scars.may <- postpartum.data$placental_scar[which(!is.na(postpartum.data$placental_scar) & 
                                                    postpartum.data$year > 2012 & 
                                                    postpartum.data$month == 5)]
prop.test(c(sum(scars.april), sum(scars.may)), c(length(scars.april), length(scars.may)))

# check for fading in CA using data from last 10 years
CA.april <- postpartum.data$CA[which(!is.na(postpartum.data$CA) & 
                                       postpartum.data$year > 2012 & 
                                       postpartum.data$month == 4)]
CA.may <- postpartum.data$CA[which(!is.na(postpartum.data$CA) & 
                                     postpartum.data$year > 2012 & 
                                     postpartum.data$month == 5)]
CA.june <- postpartum.data$CA[which(!is.na(postpartum.data$CA) & 
                                     postpartum.data$year > 2012 & 
                                     postpartum.data$month == 6)]
prop.test(c(sum(CA.april), sum(CA.may)), c(length(CA.april), length(CA.may)))
fisher.test(matrix(c(sum(CA.may), sum(CA.june), length(CA.may)-sum(CA.may), length(CA.june)-sum(CA.june)), nrow=2))


################## Prepare joint assessment data for scars & CA #########################

CA.scar.joint <- postpartum.data[which(postpartum.data$month == 4), # use only April data since scars fade
                                 c('year', 'CA', 'placental_scar')]

# compare proportions of CA between samples with missing and non-missing scar
# (use only data after 2007 since scars weren't assessed before that)
CA_scar.NA <- CA.scar.joint$CA[which(CA.scar.joint$year >= 2007 & 
                                   is.na(CA.scar.joint$placental_scar) &
                                   !is.na(CA.scar.joint$CA))]
CA_scar.nonNA <- CA.scar.joint$CA[which(CA.scar.joint$year >= 2007 & 
                                      !is.na(CA.scar.joint$placental_scar) &
                                      !is.na(CA.scar.joint$CA))]
fisher.test(matrix(c(sum(CA_scar.NA), sum(CA_scar.nonNA), length(CA_scar.NA)-sum(CA_scar.NA), 
                     length(CA_scar.nonNA)-sum(CA_scar.nonNA)), nrow=2)) # significant difference

# compare proportions of placental scars between samples with missing and non-missing CA
# (use only data after 2007 since scars weren't assessed before that)
scar_CA.NA <- CA.scar.joint$placental_scar[which(CA.scar.joint$year >= 2007 & 
                                       !is.na(CA.scar.joint$placental_scar) &
                                       is.na(CA.scar.joint$CA))]
scar_CA.nonNA <- CA.scar.joint$placental_scar[which(CA.scar.joint$year >= 2007 & 
                                          !is.na(CA.scar.joint$placental_scar) &
                                          !is.na(CA.scar.joint$CA))]
fisher.test(matrix(c(sum(scar_CA.NA), sum(scar_CA.nonNA), length(scar_CA.NA)-sum(scar_CA.NA), 
                     length(scar_CA.nonNA)-sum(scar_CA.nonNA)), nrow=2)) # too small sample size (no evidence of diff.)


# assign samples to groups based on assessment outcome
CA.scar.joint$group <- NA
CA.scar.joint$group[which(CA.scar.joint$placental_scar == 0 & CA.scar.joint$CA == 0)] <- 1
CA.scar.joint$group[which(CA.scar.joint$placental_scar == 0 & CA.scar.joint$CA == 1)] <- 2
CA.scar.joint$group[which(CA.scar.joint$placental_scar == 1 & CA.scar.joint$CA == 0)] <- 3
CA.scar.joint$group[which(CA.scar.joint$placental_scar == 1 & CA.scar.joint$CA == 1)] <- 4
CA.scar.joint$group[which(is.na(CA.scar.joint$placental_scar) & CA.scar.joint$CA == 1)] <- 5
CA.scar.joint$group[which(is.na(CA.scar.joint$placental_scar) & CA.scar.joint$CA == 0)] <- 6

# create summary table (include only data after 2007 when placental scars started being evaluated)
z.incomplete <- table(CA.scar.joint[which(CA.scar.joint$year >= 2007), c('year', 'group')])
z <- matrix(0, nrow=length(years), ncol=6)
rownames(z) <- years
colnames(z) <- 1:6
z[which(years %in% rownames(z.incomplete)),] <- z.incomplete


################## Prepare placental scar only data ##########################

scar.data <- CA.scar.joint[which(!is.na(CA.scar.joint$placental_scar) & 
                                   is.na(CA.scar.joint$CA)), c('year', 'placental_scar')]

scar.incomplete <- merge(aggregate(placental_scar ~ year, FUN=sum, data=scar.data),
                                aggregate(placental_scar ~ year, FUN=length, data=scar.data), by='year')

scars <- matrix(0, nrow=length(years), ncol=3)
colnames(scars) <- colnames(scar.incomplete) <- c('year', 'scars', 'sample_size')
scars[which(years %in% scar.incomplete$year),] <- as.matrix(scar.incomplete)
scars[,'year'] <- years


################### Prepare CA only data ############################

# CA data from April after 2007 already included in the joint CA-scar data
# CA only data includes either samples from May-June, or samples before 2007
CA.data <- postpartum.data[which(postpartum.data$age >= 5 &
                                   !is.na(postpartum.data$CA) &
                                   (postpartum.data$year < 2007 | postpartum.data$month %in% c(5:6))), c('year', 'CA')]

CA.incomplete <- merge(aggregate(CA ~ year, FUN=sum, data=CA.data),
                         aggregate(CA ~ year, FUN=length, data=CA.data), by='year')

CA <- matrix(0, nrow=length(years), ncol=3)
colnames(CA) <- colnames(CA.incomplete) <- c('year', 'CA', 'sample_size')
CA[which(years %in% CA.incomplete$year),] <- as.matrix(CA.incomplete)
CA[,'year'] <- years


############################# Save data #######################################

write.csv(pregnancies, file = paste(output_directory, 'pregnancies.csv', sep=''), row.names = F)
write.csv(z, file = paste(output_directory, 'reproductive_assessments.csv', sep=''))
write.csv(scars, file = paste(output_directory, 'scars.csv', sep=''), row.names = F)
write.csv(CA, file = paste(output_directory, 'CA.csv', sep=''), row.names = F)


