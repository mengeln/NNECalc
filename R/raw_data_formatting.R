###Read in data###
#data <- read.csv("P:\\PartTimers\\MarkEngeln\\SWAMPLabResults7_24_2012.csv")
load("data/data.RData")

###Load functions###
source("R/NNE_functions.R")

###Create smaller, workable subset####
randomIDs <- sample(data$LabSampleID[data$LabSampleID!="None"], 500)
testlabdata <- data[which(data$LabSampleID %in% randomIDs), ]

###Create fake phab data###
StationCode <- unique(testlabdata$StationCode)
Latitude <- rep(34.6, times=length(StationCode))
WaterTemperature <- rep(21, times=length(StationCode))
WaterDepth <- rep(10, times=length(StationCode))
CanopyClosure <- rep(15, times=length(StationCode))
accural <- rep(80, times=length(StationCode))
testphabdata <- data.frame(StationCode, Latitude, WaterTemperature, WaterDepth,
                           CanopyClosure, accural)

###Error check and format lab data###
testdataf <- NNEformat(testlabdata)

###Merge data###
calcdata <- merge(testdataf, testphabdata, by="StationCode", all=F, all.x=T, all.y=F)

###Calc misc data###

calcdata$SolarRadiation <- SolarRadiation(calcdata)

calcdata$LightFactor <- LightFactor(calcdata)

###Use calculators###

density_qual2k <- MaxAlgaeDensity_standardQual2k(calcdata)
benthic_qual2k <- BenthicChlora(density_qual2k)

density_qual2krevised <- MaxAlgaeDensity_revisedQual2k(calcdata)
benthic_qual2krevised <- BenthicChlora(density_qual2krevised)

density_qual2kaccural <- MaxAlgaeDensity_accrual(calcdata, density_qual2krevised)
benthic_qual2kaccural <- BenthicChlora_accrual(calcdata, benthic_qual2krevised)

###Bind results###
qual2k_results <- data.frame(density_qual2k, benthic_qual2k, density_qual2krevised, benthic_qual2krevised,
           density_qual2kaccural, benthic_qual2kaccural)
qual2k_results <- cbind(
  calcdata[, c("StationCode", "LabSampleID", "SampleDate", "Replicate")],
  qual2k_results)
View(qual2k_results)
