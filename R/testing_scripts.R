###Create a set of fake phab data for testing purposes###
fakedata <- function (data) {
  ###Create fake phab data for testing purposes###
  StationCode <- unique(data$StationCode)
  Latitude <- rep(34.6, times=length(StationCode))
  WaterTemperature <- rep(21, times=length(StationCode))
  WaterDepth <- rep(10, times=length(StationCode))
  CanopyClosure <- rep(15, times=length(StationCode))
  accrual <- rep(120, times=length(StationCode))
  Turbidity <- rep(.6, time=length(StationCode))
  LightFraction <- rep(.9, times=length(StationCode))
  testphabdata <- data.frame(StationCode, Latitude, WaterTemperature, WaterDepth,
                             CanopyClosure, accrual, Turbidity, LightFraction)
}

###Read in data###
load("data/data.RData")
validate <- read.csv("P:\\PartTimers\\MarkEngeln\\NNECalc\\Data\\validate.csv", header=T)

###Run##
source("R/NNE_functions.R")
results <- NNEcalc(data, validate)
View(results)


###Target calculator testing###
testdata <- Qual2k(SWAMPformat(data, validate))
View(Qual2k_targets(testdata, 100))

library(rrcovNA)
imputed <- impSeqRob(imputedata[,c("Latitude", "CanopyClosure", 	"WaterDepth",
                                "WaterTemperature",	"Turbidity.NTU.")])

data[,c("Latitude", "CanopyClosure",
        "WaterDepth", "WaterTemperature",	"Turbidity.NTU.")] <- imputed