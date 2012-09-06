observed <- read.csv("P:\\PartTimers\\MarkEngeln\\from Betty\\dataForNNEmodels\\observedChla.csv")
alldata <- merge(NNEdata_predictions, observed, all=F, all.x=T, all.y=F)
alldata <- alldata[which(!is.na(alldata$observedChla.mg.m2.)),]

good <- which(complete.cases(alldata$observedChla.mg.m2.) &
  alldata$nitrogen > 0 & alldata$Phosphorus.as.P>0 & alldata$OrthoPhosphate.as.P>0)

alldata <- alldata[good,]

alldata$ratio_by_weight <- alldata$nitrogen / alldata$Phosphorus.as.P
alldata$N_in_M <- (alldata$nitrogen/1000)/14.0067
alldata$P_in_M <- (alldata$Phosphorus.as.P/1000)/30.973761
alldata$ratio_in_M <- alldata$N_in_M/alldata$P_in_M

###log transform###

alldata$lognitrogen <- log10(alldata$nitrogen)
alldata$logPhosphorus.as.P <- log10(alldata$Phosphorus.as.P)
alldata$logOrthoPhosphate <- log10(alldata$OrthoPhosphate)
alldata$logChl <- log10(alldata$observedChla.mg.m2.)
# 
# # #validation <- alldata[1:100,]
# # 
# # Nlimited <- alldata[which(alldata$ratio_in_M<=10),]
# # normal <- alldata[which(alldata$ratio_in_M>10 &
# #   alldata$ratio_in_M<20),]
# # Plimited <- alldata[which(alldata$ratio_in_M>=20),]
# 
# 
# summary(lm(observedChla.mg.m2. ~ standardQual2k_BenthicChlora, data=alldata))
# plot(observedChla.mg.m2. ~ standardQual2k_BenthicChlora, data=alldata)
# 
# summary(lm(observedChla.mg.m2. ~ revisedQual2k_MaxAlgaeDensity, data=alldata))
# plot(observedChla.mg.m2. ~ revisedQual2k_MaxAlgaeDensity, data=alldata)
# 
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, data=alldata))
# plot(observedChla.mg.m2. ~ MeanAlgalDen_Dodds97, data=alldata)
# 
# summary(lm(observedChla.mg.m2. ~ MeanAlgalDen_Dodds02, data=alldata))
# plot(observedChla.mg.m2. ~ MeanAlgalDen_Dodds02, data=alldata)
# # ###Create model###
# # model <- lm(log10(observedChla.mg.m2.) ~ log10(nitrogen) +  log10(OrthoPhosphate.as.P) +
# #    WaterTemperature +  WaterDepth + CanopyClosure +
# #     ratio_in_M + Latitude
# #   , data=normal)
# # 
# # summary(model)
# # library(MASS)
# # stepAIC(model)
# 
# ####DELTA: A vector of length two. The first component is the raw cross-validation 
# #estimate of prediction error. The second component is the adjusted cross-validation 
# #estimate. The adjustment is designed to compensate for the bias introduced 
# #by not using leave-one-out cross-validation.
# 
# ###see ?cv.glm
# 
# # library(boot)
# # ###Plimited model###
# # Plimited_model <- glm(formula = log10(observedChla.mg.m2.) ~ log10(nitrogen) + 
# #   WaterTemperature + log10(Phosphorus.as.P) + Latitude, family = gaussian, 
# #                       data = Plimited)
# # cv.glm(Plimited, Plimited_model, K=3)$delta
# # 
# # ###Nlimited model
# # Nlimited_model <- glm(formula = log10(observedChla.mg.m2.) ~ log10(nitrogen) + 
# #   WaterTemperature + Latitude, family = gaussian, data = Nlimited)
# # cv.glm(Nlimited, Nlimited_model, K=3)$delta
# # 
# # ###normal model##
# # normal_model <-glm(formula = log10(observedChla.mg.m2.) ~ log10(OrthoPhosphate.as.P) + 
# #   WaterTemperature + CanopyClosure + Latitude, family = gaussian, 
# #                   data = normal)
# # cv.glm(normal, normal_model, K=3)$delta
# # 
# # 
# # ###random forest##
# #alldata2 <- alldata[complete.cases(alldata[,c(17:19, 21:22, 49:53)]),]
# library(caret)
# registerDoSEQ()
# ctrl <- rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "all")
# 
# set.seed(10)
# caretResults_nolat <- rfe(trainData[,c(18:19, 21:22, 49:53)], trainData[,54], sizes=c(2, 4, 6), rfeControl=ctrl)
# caretResults$fit
# rownames(caretResults$fit$importance)
# caretResults$fit$importance
# rfNNE <- caretResults
# # #save(rfNNE, file="Data/rfNNE.RData")
# 
# load("Data/rfNNE.RData")
# plot(log10(holdback$observedChla.mg.m2.) ~ predict(caretResults, holdback))
# summary(lm(log10(holdback$observedChla.mg.m2.) ~ predict(caretResults, holdback)))
# plotdata <- data.frame(predict(caretResults, holdback), log10(holdback$observedChla.mg.m2.))
# plotdata$accuracy <- ifelse(holdback$rf/holdback$observedChla.mg.m2. > 2, "over", ifelse(
#   holdback$rf/holdback$observedChla.mg.m2. < 0.5, "under", "good"))
# pdf(file="randomForest_plot.pdf")
# ggplot(plotdata, aes(x=log10.holdback.observedChla.mg.m2.., y=predict.caretResults..holdback.,
#                      colour=accuracy)) + geom_point() + opts(title="Random Forest model, 30% holdback, rsq=0.3575")
# dev.off()
# sample <- sample(1:length(alldata[[1]]), size=(round(length(alldata[,1])*0.7)))
# trainData <- alldata[sample,]
# holdback <- alldata[!(1:length(alldata[[1]]) %in% sample),]
# 
# 
# TNTP_model<- lm(data=trainData, logChl ~ lognitrogen + logPhosphorus.as.P)
# summary(TNTP_model)
# 
# summary(lm(holdback$logChl ~ predict(TNTP_model, holdback)))
# plot(holdback$logChl ~ predict(TNTP_model, holdback))
#      
# # ###linear model###
# # registerDoSEQ()
# # ctrl2 <- rfeControl(functions = lmFuncs, method = "cv",verbose = FALSE, returnResamp = "all")
# # 
# # set.seed(10)
# # caretLMResults <- rfe(normal[,c(17:19, 21:22, 49:50, 52)], normal[,53], sizes=c(2, 4, 6), rfeControl=ctrl2)
# # summary(caretLMResults$fit)
# # 
# # ####
# # plot(alldata[,33], 10^alldata2[,53])
