##LDA###
library(MASS)
library(ggplot2)
alldata$standardQualDiff <- rep("good", nrow(alldata))
alldata$standardQualDiff[which((alldata$standardQual2k_BenthicChlora /
  alldata$observedChla.mg.m2.) > 2)] <- "over"
alldata$standardQualDiff[which((alldata$standardQual2k_BenthicChlora /
  alldata$observedChla.mg.m2.) < .5)] <- "under"
alldata$standardQualDiff <- as.factor(alldata$standardQualDiff)

fit <- lda(standardQualDiff ~  CanopyClosure + nitrogen + OrthoPhosphate.as.P + WaterDepth +
  WaterTemperature + Turbidity + ratio_in_M, alldata, na.action="na.omit")
fit
fitpredict <- predict(fit, alldata)
plotfit <- as.data.frame(fitpredict$x)
plotfit$class <- fitpredict$class
ggplot(plotfit, aes(x=LD1, y=LD2, colour=class)) + geom_point()

variables <- c("Latitude", "CanopyClosure", "nitrogen", "OrthoPhosphate.as.P", "WaterDepth", "WaterTemperature", "Turbidity", "ratio_in_M")
manova <- manova(as.matrix(alldata[,variables]) ~ alldata$standardQualDiff)
summary(manova, test="Wilks")]
summary.aov(manova)

for(i in 32:40){
select <- which(alldata$Latitude > i & alldata$Latitude <= i + 2)
print(paste(i, "through", i+2))
print(table(alldata$standardQualDiff[select])/sum(table(alldata$standardQualDiff[select])))
}

###logistic regression

lregress <- glm(standardQualDiff ~ Latitude + CanopyClosure + nitrogen + OrthoPhosphate.as.P + WaterDepth +
  WaterTemperature + Turbidity + ratio_in_M + Latitude*CanopyClosure, alldata, family=binomial)
summary(lregress)


###PCA###
pfit <- prcomp(~ Latitude + CanopyClosure + nitrogen + OrthoPhosphate.as.P + WaterDepth +
  WaterTemperature + Turbidity + ratio_in_M, data=alldata, na.action="na.omit")

pfit

plot(predict(pfit, alldata)[,1]~ predict(pfit, alldata)[,2])

index <- -500:500
total <- sapply(index, function(i){
  table <- table(alldata$standardQualDiff[which(predict(pfit, alldata)[,1] > i)])
  sum(table)
}
)

over <- sapply(index, function(i){
  table <- table(alldata$standardQualDiff[which(predict(pfit, alldata)[,1] > i)])
  table[which(names(table) %in% "over")]
}
)
over <- over/total
plot(over ~ index)

good <- sapply(index, function(i){
  table <- table(alldata$standardQualDiff[which(predict(pfit, alldata)[,1] > i)])
  table[which(names(table) %in% "good")]
}
)
good <- good/total
plot(good ~ index)

under<- sapply(index, function(i){
  table <- table(alldata$standardQualDiff[which(predict(pfit, alldata)[,1] > i)])
  table[which(names(table) %in% "under")]
}
)
under <- under/total
plot(under ~ index)

###Exclusion plots###
pdf(file="Removed_observations.pdf")

remove <- c(100, 150, 200)
for(i in 1:3){
plot(standardQual2k_BenthicChlora ~ observedChla.mg.m2., 
     data=alldata[which(alldata$observedChla.mg.m2 < remove[i]),], main=paste("Observed above", remove[i], "removed"),
     sub=paste("rsq=", qualr[i]))
}
for(i in 1:3){
  plot(log10(standardQual2k_BenthicChlora) ~ log10(observedChla.mg.m2.), 
       data=alldata[which(alldata$observedChla.mg.m2 < remove[i]),], main=paste("Observed above", remove[i], "removed"),
       sub=paste("rsq=", qualrl[i]))
}

for(i in 1:3){
  plot(MeanBenthicChlA_Dodds97 ~ observedChla.mg.m2., 
       data=alldata[which(alldata$observedChla.mg.m2 < remove[i]),], main=paste("Observed above", remove[i], "removed"),
       sub=paste("rsq=", doddsr[i]))
}
for(i in 1:3){
  plot(log10(MeanBenthicChlA_Dodds97) ~ log10(observedChla.mg.m2.), 
       data=alldata[which(alldata$observedChla.mg.m2 < remove[i]),], main=paste("Observed above", remove[i], "removed"),
       sub=paste("rsq=", doddsrl[i]))
}

dev.off()

for(i in c(100, 150, 200)){
  print(summary(lm(log10(standardQual2k_BenthicChlora) ~ log10(observedChla.mg.m2.), 
       data=alldata[which(alldata$observedChla.mg.m2 < i),])))
}

qualr <- c(0.05955, 0.1282, 0.1781)
qualrl <- c(0.0405, 0.07263, 0.09406)
doddsr <- c(0.04713, 0.08983, 0.1154)
doddsrl <- c(0.02332, 0.04207, 0.05599)

###Latitude plots###
library(ggplot2)
pdf(file="hist.pdf")
ggplot(alldata, aes(x=Latitude, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=CanopyClosure, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=nitrogen, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=OrthoPhosphate.as.P, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=WaterDepth, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=WaterTemperature, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=Turbidity, fill=standardQualDiff)) + geom_bar()
ggplot(alldata, aes(x=ratio_in_M, fill=standardQualDiff)) + geom_bar()
}
dev.off()