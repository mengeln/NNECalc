###read
validate <- read.csv("P:\\PartTimers\\MarkEngeln\\NNECalc\\Data\\validate.csv", header=T)
results <- read.csv("P:\\PartTimers\\MarkEngeln\\NNECalc\\Data\\results.csv", header=T)



###Merge Data together###
alldata <- merge(results, validate, all=F, all.x=T, all.y=F, by.x=(c("StationCode", "SampleDate")),
                 by.y=(c("StationCode", "SampleDate")))

###Remove missing data###
alldata<- na.omit(alldata)
good <- which(complete.cases(alldata$chla_mg_per_m2) &
  alldata$nitrogen > 0 & alldata$Phosphorus.as.P>0 & alldata$OrthoPhosphate.as.P>0)

alldata$ratio_by_weight <- alldata$nitrogen / alldata$Phosphorus.as.P
alldata$N_in_M <- (alldata$nitrogen/1000)/14.0067
alldata$P_in_M <- (alldata$Phosphorus.as.P/1000)/30.973761
alldata$ratio_in_M <- mean(alldata$N_in_M/alldata$P_in_M)
range(alldata$N_in_M/alldata$P_in_M)

###Create model###
model <- lm(log10(chla_mg_per_m2) ~ log10(nitrogen) +  log10(OrthoPhosphate.as.P) +
  log10(Phosphorus.as.P) + Temperature.C. + Turbidity.NTU. + Depth.m. + CanopyCover... +
  Velocity.m.s. + ratio_by_weight + ratio_in_M + Latitude + Phi_lb + SolarRadiation
  , data=alldata[good,])

summary(model)
library(MASS)
stepAIC(model)

model1 <- lm(log10(chla_mg_per_m2) ~  ratio_by_weight + Velocity.m.s. + log10(OrthoPhosphate.as.P)
             + Latitude
             , data=alldata[good,])
model2 <- lm(log10(chla_mg_per_m2) ~ log10(nitrogen) + log10(Phosphorus.as.P), data=alldata[good,])
anova(model1, model2)


#  plot(log10(Phosphorus.as.P) ~  log10(OrthoPhosphate.as.P)
#       , data=alldata[good,])
 
# plot(log10(chla_mg_per_m2) ~  log10(Phosphorus.as.P)
#      , data=alldata[good,])
# 
# plot(log10(chla_mg_per_m2) ~ log10(nitrogen), data=alldata[good,])