source("r/NNE_functions.r")
source("r/sketch.r")
source("r/NNE model.r")
library(randomForest)
###Plots###
library(ggplot2)
ggplot(alldata[alldata$observedChla.mg.m2. < 600,], aes(x=observedChla.mg.m2.,
                                                        y=standardQual2k_BenthicChlora, colour=Latitude)) + 
                                                          geom_point() + stat_smooth()

ggplot(alldata[alldata$ratio_in_M < 50,], aes(x=observedChla.mg.m2., y=standardQual2k_BenthicChlora, colour=ratio_in_M)) + 
  geom_point() + stat_smooth()

ggplot(alldata[alldata$observedChla.mg.m2. < 600,], aes(x=observedChla.mg.m2.,
                                                        y=MeanBenthicChlA_Dodds97, colour=Latitude)) + 
                                                          geom_point() + stat_smooth()

ggplot(alldata[alldata$ratio_in_M < 50,], aes(x=observedChla.mg.m2., y=MeanBenthicChlA_Dodds97, colour=ratio_in_M)) + 
  geom_point() + stat_smooth()



###Qual2K Variable Importance###
set.seed(1234)
alldata$DWC_Month <- as.factor(alldata$DWC_Month)
rf_qual <- randomForest(observedChla.mg.m2. ~ nitrogen + OrthoPhosphate.as.P + DWC_Month +
  Latitude + CanopyClosure + WaterDepth + WaterTemperature + Turbidity, data=alldata,
                        ntree=500, importance=T, proximity=T)
importance(rf_qual)

\varImpPlot(rf_qual)
MDSplot(rf_qual, alldata$observedChla.mg.m2.)

###Dodds variable importance

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
  summary(lm(observedChla.mg.m2. ~ MeanBenthicChlA_Dodds97, 
             alldata[alldata$CanopyClosure<(i+15) & alldata$CanopyClosure>i,]))$r.squared
})
canopycover <- sapply(1:85, function(i)mean(c(i+15, i)))
plot(chla_canopy_rsq ~ canopycover)

chla_canopy_msq <- sapply(1:85, function(i){
  var(alldata$observedChla.mg.m2.[which(alldata$CanopyClosure<(i+15) & alldata$CanopyClosure>i)],
      alldata$MeanBenthicChlA_Dodds97[which(alldata$CanopyClosure<(i+15) & alldata$CanopyClosure>i)])
})
plot(chla_canopy_msq ~ canopycover)
###Latitude
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
#            alldata[alldata$Latitude<34,]))
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
#            alldata[alldata$Latitude<36,]))
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
#            alldata[alldata$Latitude<38 & alldata$Latitude>36,]))
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
#            alldata[alldata$Latitude>38,]))
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
#            alldata[alldata$Latitude>40,]))
chla_latitude_rsq <- sapply(32:40, function(i){
  summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, 
             alldata[alldata$Latitude<(i+2) & alldata$Latitude>i,]))$r.squared
})
latitude <- sapply(32:40, function(i)mean(c(i+2, i)))
plot(chla_latitude_rsq ~ latitude)

###Turbidity###
chla_Turbidity_msq <- sapply(1:20, function(i){
  var(alldata$observedChla.mg.m2.[which(alldata$Turbidity<(i*3/20 + .3) & alldata$Turbidity>(i*3/20))],
      alldata$MeanBenthicChlA_Dodds97[which(alldata$Turbidity<(i*3/20 + .3) & alldata$Turbidity>(i*3/20))])
})

turbid <- sapply(1:20, function(i)mean(c((i*3/20 + .3), (i*3/20))))
plot(chla_Turbidity_msq ~ turbid)


###Polynomial Regression###
set.seed(1234)
validate <- sample(1:596, 596*.7)
holdout <- !(1:596 %in% validate)

poly <- lm(logChl ~ poly(lognitrogen, 2) + logPhosphorus.as.P, data=alldata[validate,])
summary(lm(alldata$logChl[holdout] ~ predict(poly, alldata[holdout,])))

linear <- lm(logChl ~ lognitrogen + logPhosphorus.as.P, data=alldata[validate,])
summary(lm(alldata$logChl[holdout] ~ predict(linear, alldata[holdout,])))


vwdata <- data.frame(alldata$logChl, predict(poly, alldata))
colnames(vwdata) <- c("Observed", "Predicted")
vwReg(Predicted ~ Observed, data=vwdata,
      shade.alpha=0, slices=400, palette=colorRampPalette(c("black", "green",
                                                            "yellow", "red"), bias=5)(20),
      method=lm)

