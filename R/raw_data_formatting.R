###Read in data###
data <- read.csv("P:\\PartTimers\\MarkEngeln\\SWAMPLabResults7_24_2012.csv")


###Select Relevant Columns###
workingdata<-data[, c("StationCode", "LabSampleID", "SampleDate", "DWC_Month", "Replicate",
                      "AnalyteName", "Result" )]

###Convert from flat to wide format###
library(reshape)
workingdata <- cast(workingdata, formula = StationCode + LabSampleID + SampleDate + DWC_Month +
  Replicate ~ AnalyteName, value="Result", fun.aggregate=sum)
workingdata <- data.frame(workingdata)

###Sum for total nitrogen###
workingdata$nitrogen <- apply(workingdata[, 6:11], 1, sum)

###Remove observations with missing nitrogen###
workingdata$nitrogen[workingdata$nitrogen==-88]<-NA
workingdata <- workingdata[!is.na(workingdata$nitrogen),]

###Check for P errors and remove###
Pcheck <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0)
Pcheck2 <- which(workingdata$OrthoPhosphate.as.P>workingdata$Phosphorus.as.P)
workingdata <- workingdata[!(1:length(workingdata[[1]]) %in% intersect(Pcheck, Pcheck2)),]

###Calculate Organic P###
Pcheck <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0)
workingdata$OrganicP <- rep(NA, length(workingdata[[1]]))
workingdata$OrganicP[Pcheck] <- workingdata$Phosphorus.as.P[Pcheck] - 
  workingdata$OrthoPhosphate.as.P[Pcheck]
workingdata$OrganicP[is.na(workingdata$OrganicP)] <- 0


###Check for outliers###

outliercheck(workingdata, 12:15)


