
fakedata <- function (data) {
  ###Create fake phab data for testing purposes###
  StationCode <- unique(data$StationCode)
  Latitude <- rep(34.6, times=length(StationCode))
  WaterTemperature <- rep(21, times=length(StationCode))
  WaterDepth <- rep(10, times=length(StationCode))
  CanopyClosure <- rep(15, times=length(StationCode))
  accrual <- rep(80, times=length(StationCode))
  Turbidity <- rep(.6, time=length(StationCode))
  LightFraction <- rep(.9, times=length(StationCode))
  testphabdata <- data.frame(StationCode, Latitude, WaterTemperature, WaterDepth,
                             CanopyClosure, accrual, Turbidity, LightFraction)
}

###Read in data###
load("data/data.RData")
validate <- read.csv("P:\\PartTimers\\MarkEngeln\\NNECalc\\Data\\validate.csv", header=T)

###Run##
results <- Qual2k(data, validate)
#View(results[, c("StationCode", "LabSampleID", "SampleDate", 
                 #"Replicate", "standardQual2k_MaxAlgaeDensity", "standardQual2k_BenthicChlora")])
View(results[results$standardQual2k_MaxAlgaeDensity != 0, ])
