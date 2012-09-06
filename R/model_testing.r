source("r/NNE_functions.r")
source("r/sketch.r")
source("r/NNE model.r")
library(randomForest)
###Dodds variable importance
set.seed(1234)
rf <- randomForest(logChl ~ lognitrogen + logPhosphorus.as.P,
                   data=alldata, ntree=500, importance=T)
importance(rf)

rf_all <- randomForest(observedChla.mg.m2. ~ nitrogen + Phosphorus.as.P,
                   data=alldata, ntree=500, importance=T)
importance(rf_all)

rf_lat <- randomForest(logChl ~ lognitrogen + logPhosphorus.as.P + Latitude,
                       data=alldata, ntree=500, importance=T)
importance(rf_lat)

###Canopy Cover
chla_canopy_rsq <- sapply(1:85, function(i){
  summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
             alldata[alldata$CanopyClosure<(i+15) & alldata$CanopyClosure>i,]))$r.squared
})
canopycover <- sapply(1:85, function(i)mean(c(i+15, i)))
plot(chla_canopy_rsq ~ canopycover)

###Latitude
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$Latitude<34,]))
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$Latitude<36,]))
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$Latitude<38 & alldata$Latitude>36,]))
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$Latitude>38,]))
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$Latitude>40,]))
chla_latitude_rsq <- sapply(32:40, function(i){
  summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
             alldata[alldata$Latitude<(i+2) & alldata$Latitude>i,]))$r.squared
})
latitude <- sapply(32:40, function(i)mean(c(i+2, i)))
plot(chla_latitude_rsq ~ latitude)

###Water Depth###
summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
           alldata[alldata$WaterDepth<.11 & alldata$WaterDepth>.07,]))