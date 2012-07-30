###Read in data###
load("data/data.RData")


smalldata <- function (data) {
  ###Create smaller, workable subset
  randomIDs <- sample(data$LabSampleID[data$LabSampleID!="None"], 500)
  testlabdata <- data[which(data$LabSampleID %in% randomIDs), ]
  return(testlabdata)
}

fakedata <- function (data) {
  ###Create fake phab data###
  StationCode <- unique(data$StationCode)
  Latitude <- rep(34.6, times=length(StationCode))
  WaterTemperature <- rep(21, times=length(StationCode))
  WaterDepth <- rep(10, times=length(StationCode))
  CanopyClosure <- rep(15, times=length(StationCode))
  accural <- rep(80, times=length(StationCode))
  testphabdata <- data.frame(StationCode, Latitude, WaterTemperature, WaterDepth,
                             CanopyClosure, accural)
}

runall <- function (testlabdata, testphabdata) {
  ###Load functions###
  source("R/NNE_functions.R")
  
  ###Error check and format lab data###
  testdataf <- NNEformat(testlabdata)
  
  ###Merge data###
  calcdata <- merge(testdataf, testphabdata, by="StationCode", all=F, all.x=T, all.y=F)
  
  ###Calc misc data###
  
  calcdata$SolarRadiation <- SolarRadiation(calcdata)
  
  calcdata$LightFactor <- LightFactor(calcdata)
  
  ###Use calculators###
  
  density_qual2k <- MaxAlgaeDensity_standardQual2k(calcdata)
  standardQual2k_BenthicChlora <- BenthicChlora(density_qual2k)
  
  density_qual2krevised <- MaxAlgaeDensity_revisedQual2k(calcdata)
  RevisedQual2k_BenthicChlora <- RevisedQual2k_BenthicChlora(density_qual2krevised)
  
  density_qual2kaccural <- MaxAlgaeDensity_accrual(calcdata, density_qual2krevised)
  benthic_qual2kaccural <- BenthicChlora_accrual(calcdata, RevisedQual2k_BenthicChlora)
  
  ###Bind results###
  qual2k_results <- data.frame(density_qual2k, standardQual2k_BenthicChlora, density_qual2krevised, RevisedQual2k_BenthicChlora,
            density_qual2kaccural, benthic_qual2kaccural)
  qual2k_results <- cbind(
    calcdata[, c("StationCode", "LabSampleID", "SampleDate", "Replicate")],
    qual2k_results)
  return(merge(calcdata, qual2k_results, all=F, all.y=T, all.x=T))
}

testdata <- smalldata(data)
results <- runall(testdata, fakedata(testdata))
#View(results[, c("StationCode", "LabSampleID", "SampleDate", 
                 #"Replicate", "standardQual2k_MaxAlgaeDensity", "standardQual2k_BenthicChlora")])
View(results[results$standardQual2k_MaxAlgaeDensity != 0, ])
