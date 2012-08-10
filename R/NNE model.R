observed <- read.delim("clipboard", header=T, stringsAsFactors=F)
alldata <- merge(NNEdata_predictions, observed, all=F, all.x=T, all.y=F)
alldata <- alldata[which(!is.na(alldata$observedChla.mg.m2.)),]

good <- which(complete.cases(alldata$observedChla.mg.m2.) &
  alldata$nitrogen > 0 & alldata$Phosphorus.as.P>0 & alldata$OrthoPhosphate.as.P>0)

alldata <- alldata[good,]

alldata$ratio_by_weight <- alldata$nitrogen / alldata$Phosphorus.as.P
alldata$N_in_M <- (alldata$nitrogen/1000)/14.0067
alldata$P_in_M <- (alldata$Phosphorus.as.P/1000)/30.973761
alldata$ratio_in_M <- alldata$N_in_M/alldata$P_in_M

#validation <- alldata[1:100,]

Nlimited <- alldata[which(alldata$ratio_in_M<=10),]
normal <- alldata[which(alldata$ratio_in_M>10 &
  alldata$ratio_in_M<20),]
Plimited <- alldata[which(alldata$ratio_in_M>=20),]


###Create model###
model <- glm(log10(observedChla.mg.m2.) ~ log10(nitrogen) +  log10(OrthoPhosphate.as.P) +
   WaterTemperature +  WaterDepth + CanopyClosure + log10(Phosphorus.as.P) +
    ratio_in_M + Latitude + SolarRadiation
  ,family = gaussian, data=Nlimited)

summary(model)
library(MASS)
stepAIC(model)

####DELTA: A vector of length two. The first component is the raw cross-validation 
#estimate of prediction error. The second component is the adjusted cross-validation 
#estimate. The adjustment is designed to compensate for the bias introduced 
#by not using leave-one-out cross-validation.

###see ?cv.glm

library(boot)
###Plimited model###
Plimited_model <- glm(formula = log10(observedChla.mg.m2.) ~ log10(nitrogen) + 
  WaterTemperature + log10(Phosphorus.as.P) + Latitude, family = gaussian, 
                      data = Plimited)
cv.glm(Plimited, Plimited_model, K=3)$delta

###Nlimited model
Nlimited_model <- glm(formula = log10(observedChla.mg.m2.) ~ log10(OrthoPhosphate.as.P) + 
  ratio_in_M, , family = gaussian, data = Nlimited)
cv.glm(Nlimited, Nlimited_model, K=3)$delta

###normal model##
normal_model <-glm(formula = log10(observedChla.mg.m2.) ~ log10(OrthoPhosphate.as.P) + 
  WaterTemperature + CanopyClosure + Latitude, family = gaussian, 
                  data = normal)
cv.glm(normal, normal_model, K=3)$delta


library(neuralnet)
dev.off()
nn <- neuralnet(log10(observedChla.mg.m2.) ~ log10(nitrogen) +  log10(OrthoPhosphate.as.P) +
  +     WaterTemperature +  WaterDepth + CanopyClosure + log10(Phosphorus.as.P) +
  +     ratio_in_M + Latitude + SolarRadiation, hidden = 2,
                linear.output = FALSE, data=normal)

plot.nn(nn)
dev.off()
prediction(net.sum, list.glm=list(main=main, full=full))
