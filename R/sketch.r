###Calculate all models from NNE data###
NNEdata <- read.csv("C:\\Documents and Settings\\gisuser\\Desktop\\NNEdata.csv")
workingNNEdata <- NNEdata[which(apply(is.na(NNEdata[, c(17:19, 21)]), 1, sum)==0),]
workingNNEdata$WaterDepth[workingNNEdata$Project=="SMC"] <- 
  workingNNEdata$WaterDepth[workingNNEdata$Project=="SMC"]/100
workingNNEdata <- miscData(workingNNEdata)
NNEdata_predictions <- merge(Qual2k(workingNNEdata), Dodds(workingNNEdata))

write.csv(NNEdata_predictions, 
          file="C:\\Documents and Settings\\gisuser\\Desktop\\NNEcalc_predictions.csv")

###SWAMP PHAB formatting###

SWAMPturbid <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\SWAMP\\turbiditySWAMP.csv")

SWAMPtemp <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\SWAMP\\tempSWAMP.csv")

SWAMPcanopy <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\SWAMP\\SWAMP PHab calculated.csv")

library(reshape)
SWAMPcanopy <- data.frame(cast(data=SWAMPcanopy, StationCode + SampleDate + Latitude~ Phab_variable,
                    value="Result", fun.aggregate=mean))

SWAMPtemp <- data.frame(cast(data=SWAMPtemp, StationCode + SampleDate ~ AnalyteName,
                              value="Result", fun.aggregate=mean))

SWAMPturbid <- data.frame(cast(data=SWAMPturbid, StationCode + SampleDate ~ AnalyteName,
                             value="Result", fun.aggregate=mean))

SWAMPphab <- merge(merge(SWAMPcanopy, SWAMPtemp), SWAMPturbid)

SWAMPchem <- read.csv("C:\\Documents and Settings\\gisuser\\Desktop\\SWAMPchem.csv")

###change SWAMPformat so that it outputs SWAMPchem2 global variable before merging###
SWAMPdata <- SWAMPformat(SWAMPchem, SWAMPphab)

fullSWAMPdata <- data.frame(merge(SWAMPchem2, SWAMPphab, all=T))


###SMC###
SMCchem <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\SMC\\tblExtract_SMC_ChemistryResults.csv")
SMCphab <- read.csv("C:\\Documents and Settings\\gisuser\\Desktop\\SMCphab.csv")
###change SWAMPformat so that it outputs SMCchem2 global variable before merging###
SMCdata <- SWAMPformat(SMCchem, SMCphab)
SMCchem2$StationCode <- as.character(SMCchem2$StationCode)
SMCchem2$SampleDate <- as.character(SMCchem2$SampleDate)
SMCphab$StationCode <- as.character(SMCphab$StationCode)
SMCphab$SampleDate <- as.character(SMCphab$SampleDate)
fullSMCdata <- data.frame(merge(SMCchem2, SMCphab, all=F, all.x=T, all.y=F, by=c(
  "StationCode", "SampleDate")))
fullSMCdata$WaterDepth <- fullSMCdata$WaterDepth/100
SMCdata$WaterDepth <- SMCdata$WaterDepth/100
###Prop 50###
p50 <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\Prop50\\Prop50_formatted.csv")

p50 <- colRename(p50)

###Make master table###
fullSMCdata <- fullSMCdata[,!(colnames(fullSMCdata) %in% c("X"))]
colnames(fullSWAMPdata)[18:19] <- c("CanopyClosure", "WaterDepth")
colnames(fullSWAMPdata)[21] <- "WaterTemperature"
colnames(p50)[15] <- "StreamVelocity.m.s."
colnames(fullSMCdata)[9] <- "OrthoPhosphate.as.P"
colnames(fullSMCdata)[8] <- "Nitrogen..Total"
p50 <- p50[,!(colnames(p50) %in% c("SampleID"))]
p50$StationCode <- as.character(p50$StationCode)
p50$SampleDate <- as.character(p50$SampleDate)

library(plyr)
NNEdata <- rbind.fill(fullSWAMPdata, fullSMCdata, p50)
write.csv(NNEdata, file="C:\\Documents and Settings\\gisuser\\Desktop\\NNEdata.csv")

###Usable data###
workingNNEdata <- NNEdata[which(apply(is.na(NNEdata[, c(17:19, 21)]), 1, sum)==0),]
workingNNEdata <- miscData(workingNNEdata)
NNEdata_predictions <- merge(Qual2k(workingNNEdata), Dodds(workingNNEdata))

write.csv(NNEdata_predictions, 
          file="C:\\Documents and Settings\\gisuser\\Desktop\\NNEcalc_predictions.csv")

validation <- read.delim("clipboard", header=T)

fullvalidation <- merge(NNEdata_predictions, validation, all=F)

plot(data=fullvalidation, spreadsheet_Standard.QUAL2K ~standardQual2k_BenthicChlora)

plot(data=fullvalidation, RevisedQual2k_BenthicChlora ~spreadsheet_Revised.QUAL2K)

plot(data=fullvalidation, MaxAlgalDen_Dodds02 ~spreadsheet_Dodds..02..06..max.Chl.a)

plot(data=fullvalidation, benthic_qual2kaccrual ~spreadsheet_Revised.QUAL2K.with.accrual.adj)