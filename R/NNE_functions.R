Qual2k <- function (testlabdata, testphabdata) {
  ###Error check and format lab data###
  testdataf <- NNEformat(testlabdata)
  
  ###Merge data###
  calcdata <- merge(testdataf, testphabdata, by=c("StationCode", "SampleDate"),
                                                  all=F, all.x=T, all.y=F)
  calcdata <- na.omit(calcdata)
  ###Fix column names###
  colnames(calcdata)[grep("Canopy", colnames(calcdata))] <- "CanopyClosure"
  colnames(calcdata)[grep("Depth", colnames(calcdata))] <- "WaterDepth"
  colnames(calcdata)[grep("Temperature", colnames(calcdata))] <- "WaterTemperature"
  
  ###Calc misc data###
  
  calcdata$SolarRadiation <- SolarRadiation(calcdata)
  
  calcdata$LightFactor <- LightFactor(calcdata)
  
  calcdata$LightExtinction <- LightExtinction(calcdata)
  
  ###Use calculators###
  
  density_qual2k <- MaxAlgaeDensity_standardQual2k(calcdata)
  standardQual2k_BenthicChlora <- BenthicChlora(density_qual2k)
  
  density_qual2krevised <- MaxAlgaeDensity_revisedQual2k(calcdata)
  RevisedQual2k_BenthicChlora <- RevisedQual2k_BenthicChlora(density_qual2krevised)
  
  density_qual2kaccural <- MaxAlgaeDensity_accrual(calcdata, density_qual2krevised)
  benthic_qual2kaccural <- BenthicChlora_accrual(calcdata, RevisedQual2k_BenthicChlora)
  
  ###Bind results###
  qual2k_results <- data.frame(density_qual2k, standardQual2k_BenthicChlora, density_qual2krevised, RevisedQual2k_BenthicChlora,
                               density_qual2kaccural, benthic_qual2kaccural)
  qual2k_results <- cbind(
    calcdata[, c("StationCode", "SampleDate", "CollectionTime", "Replicate")],
    qual2k_results)
  return(merge(calcdata, qual2k_results, all=F, all.y=T, all.x=T))
}


outliercheck <- function(data){
  Result <- vector()
  AnalyteName <- vector()
  CONCATENATE <- vector()
  columns <- which(colnames(data) %in% c("Ammonia.as.N", "Nitrate...Nitrite.as.N",
                                         "Nitrate.as.N", "Nitrite.as.N", "Nitrogen..Total.Kjeldahl"))
  ###Calculate z scores for all observations in selected columns###
  for(j in columns){
    cmean <- .colMeans(data[, j], m=length(data[, j]), n=1, na.rm=T)
    csd <- sd(data[, j], na.rm=T)
    problems <- which(sapply(1:length(data[, j]), function(i){
      (data[i, j]-cmean)/csd})>2)
    ###Concatenate the results from each loop together###
    CONCATENATE <- c(CONCATENATE, data[problems, "CONCATENATE"])
    Result <- c(Result, data[problems, j])
    AnalyteName <- c(AnalyteName, rep(colnames(data)[j], length(data[problems, j])))
  }
  ###Bind columns together## 
  return(cbind(CONCATENATE, AnalyteName, Result))
}

NNEformat <- function (data) {
  ###Clean up SWAMP data table###
  
  data <- data[which(data$SampleTypeName %in% c("Grab", "Integrate")),]
  data_aggregate <- aggregate(Result~StationCode+SampleDate+CollectionTime+
    AnalyteName, data=data, FUN=mean, na.rm=T)
  data_aggregate$CONCATENATE <- paste(data_aggregate$StationCode, data_aggregate$SampleDate,
                                      data_aggregate$CollectionTime)
  data$CONCATENATE <- paste(data$StationCode, data$SampleDate,
                            data$CollectionTime)
  data_aggregate$DWC_Month <- data$DWC_Month[match(data_aggregate$CONCATENATE, data$CONCATENATE)]
  data_aggregate$Replicate <- data$Replicate[match(data_aggregate$CONCATENATE, data$CONCATENATE)]
  
  data <- data_aggregate
  data <- data[order(data$CONCATENATE),]
  
  ###Select Relevant Columns###
  
  workingdata<-data[, c("CONCATENATE", "AnalyteName", "Result" )]
  
  ###Convert from flat to wide format###
  if(require(reshape)==F)
  {
    install.packages("reshape")
    require(reshape)
  }
  workingdata <- cast(workingdata, formula = CONCATENATE ~ AnalyteName, 
                      value="Result", fun.aggregate=sum)
  workingdata <- data.frame(workingdata)
  workingdata <- merge(workingdata, data[, c("CONCATENATE", "StationCode", "DWC_Month",
                                                 "SampleDate", "CollectionTime", "Replicate")],
                       all=F, all.x=T, all.y=F)
  values <- c("Ammonia.as.N", "Nitrate...Nitrite.as.N",
              "Nitrate.as.N", "Nitrite.as.N", "Nitrogen..Total.Kjeldahl")
  

  workingdata <- workingdata[!duplicated(workingdata[,"CONCATENATE"]),]

  ###Remove observations with missing values###
  for(i in which(colnames(workingdata) %in% c(values, 
                                              "OrthoPhosphate.as.P", "Phosphorus.as.P"))){
    if(length(which(workingdata[, i]<0))>0){
    workingdata[which(workingdata[, i]<0), i] <- NA }
  }
  
  ###Sum for total nitrogen###
  workingdata$nitrogen <- apply(workingdata[,values] , 1, sum)
  workingdata$nitrogen[which(workingdata$Nitrogen..Total>0)] <- 
    workingdata$Nitrogen..Total[which(workingdata$Nitrogen..Total>0)]
  workingdata <- workingdata[!is.na(workingdata$nitrogen),]
  
  ###Check for P errors and remove###
  Pcheck <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0)
  Pcheck2 <- which(workingdata$OrthoPhosphate.as.P>workingdata$Phosphorus.as.P)
  workingdata <- workingdata[!(1:length(workingdata[[1]]) %in% intersect(Pcheck, Pcheck2)),]
  
  ###Fill in Total P###
  pfill <- which(workingdata$Phosphorus.as.P==0 & workingdata$OrthoPhosphate.as.P>0)
  workingdata$Phosphorus.as.P[pfill] <- workingdata$OrthoPhosphate.as.P[pfill]
  
  ###Calculate Organic P###
  Pcheck <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0)
  workingdata$OrganicP <- rep(NA, length(workingdata[[1]]))
  workingdata$OrganicP[Pcheck] <- workingdata$Phosphorus.as.P[Pcheck] - 
    workingdata$OrthoPhosphate.as.P[Pcheck]
  workingdata$OrganicP[is.na(workingdata$OrganicP)] <- 0
  
  
  ###Check for outliers###
  
  if(require(RGtk2Extras)==F)
  {
    install.packages("RGtl2Extras")
    require(Gtk2Extras)
  }
  
  workingdata <- workingdata[which(workingdata$nitrogen<20),]
  workingdata <- workingdata[which(workingdata$Phosphorus.as.P<20),]

  outliers <- outliercheck(workingdata)
  fix <- which(workingdata$CONCATENATE %in% outliers[, 1])
  fixcolumns <- c("OrthoPhosphate.as.P", "Phosphorus.as.P", "Ammonia.as.N", "Nitrate...Nitrite.as.N",
                  "Nitrate.as.N", "Nitrite.as.N", "Nitrogen..Total.Kjeldahl", "Nitrogen..Total")
  workingdata[fix, fixcolumns] <- dfedit(workingdata[fix, fixcolumns])
  workingdata <- workingdata[intersect(intersect(which(!is.na(workingdata$OrthoPhosphate.as.P)), 
                                       which(!is.na(workingdata$OrganicP))), intersect(
                                       which(!is.na(workingdata$Phosphorus.as.P)),
                                       which(!is.na(workingdata$nitrogen)))),]

  return(workingdata)
}

SolarRadiation <- function (data) {
  load("data/radiationTable.RData")
  data$DWC_Month <- as.character(data$DWC_Month)
  lowerBound <- as.numeric(sapply(1:length(data$StationCode), function(i){
    radiationTable[radiationTable$Month==floor(data$Latitude[i]), data$DWC_Month[i]]
  }))
  upperBound <- as.numeric(sapply(1:length(data$StationCode), function(i){
    radiationTable[radiationTable$Month==ceiling(data$Latitude[i]), data$DWC_Month[i]]
  }))
  return(sapply(1:length(lowerBound), function(i){
    upperBound[i] + (upperBound[i]-lowerBound[i])*(data$Latitude[i]-floor(data$Latitude[i]))
  }))
}

LightFactor <- function(data){
  if(is.null(data$LightFraction)){
    fraction <- .9
  } else
  {fraction <- data$LightFraction}
  return(1 - fraction * (1 - (1 - data$CanopyClosure / 100) ^ 2) ^ 0.5)
}

LightExtinction <- function(data){
  if(is.null(data$Turbidity)){
    Turbidity <- rep(.6, times=length(data$StationCode))
  } else
  {Turbidity <- data$Turbidity}
  return((.1 * Turbidity) + .44)
}

MaxAlgaeDensity_standardQual2k <- function(data){
  load("data/parameters.RData")
  ###Define Phi_Nb###
  Nitrogen <- data$nitrogen
  phos <- data$OrthoPhosphate.as.P
  N <- Nitrogen/(Nitrogen+parameters["Inorg_N_Half_Sat"])
  P <- phos/(phos+parameters["Inorg_P_Half_Sat"])
  
  Phi_NB <- sapply(1:length(Nitrogen), function(i){min(c(N[i], P[i]), na.rm=T)})
  
  ###Define Phi_lb###
  numerator <- data$SolarRadiation * data$LightFactor * exp(-1 *
    data$LightExtinction * (data$WaterDepth))
  Phi_lb <- numerator/(numerator + parameters["Light_Half_Sat"])
  
  ###Calculate metrics###  
  standardQual2k_MaxAlgaeDensity <- (parameters["Max_Growth20_standard"]* Phi_NB * Phi_lb)  * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature - 20)) / 
    (parameters["Respiration"] + parameters["Natural_Death_standard"])
  return(data.frame(standardQual2k_MaxAlgaeDensity, Phi_lb))
}

###Benthic Chlor a
BenthicChlora <- function(algae_input){
  load("data/parameters.RData")
  standardQual2k_BenthicChlora <- algae_input[,1] * parameters["C_AFDW_ratio"]
  return(standardQual2k_BenthicChlora)
}

###Revised QUAL2K Method###
MaxAlgaeDensity_revisedQual2k <- function(data){
  load("data/parameters.RData")
  ###Define Phi_Nb###
  Nitrogen <- data$nitrogen
  phos <- data$Phosphorus.as.P
  TNdenom <- parameters["gamma_TN"] / 
    (1+(exp(-log10(Nitrogen*1000)*parameters["beta_TN"] + parameters["alpha_TN"])))
  TNlimit <- 1 / (1-TNdenom)
  
  TPdenom <- parameters["gamma_TP"] / 
    (1+(exp(-log10(phos*1000)*parameters["beta_TP"] + parameters["alpha_TP"])))
  TPlimit <- 1 / (1-TPdenom)
  
  N <- Nitrogen/(Nitrogen+(parameters["TN_Half_Sat"]*TNlimit))
  P <- phos/(phos+(parameters["TP_Half_Sat"]*TPlimit))
  
  Phi_NB <- sapply(1:length(Nitrogen), function(i){min(c(N[i], P[i]), na.rm=T)})
  
  ###Define Phi_lb###
  numerator <- data$SolarRadiation * data$LightFactor * exp(-1 * 
    data$LightExtinction * (data$WaterDepth))
  Phi_lb <- numerator/(numerator + parameters["Light_Half_Sat"])
  
  ###Calculate metrics###  
  revisedQual2k_MaxAlgaeDensity <- (parameters["Max_Growth20_revised"]* Phi_NB * Phi_lb)  * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature-20)) / 
    (parameters["Respiration"] + parameters["Natural_Death_revised"])
  return(data.frame(revisedQual2k_MaxAlgaeDensity, Phi_lb))
}

###Revised Benthic Chlor a
RevisedQual2k_BenthicChlora <- function(algae_input){
  load("data/parameters.RData")
  RevisedQual2k_BenthicChlora <- algae_input[,1] * parameters["C_AFDW_ratio"]
  return(RevisedQual2k_BenthicChlora)
}

###RivsedQual2K Method, Accrual Adjustment###
MaxAlgaeDensity_accrual <- function (data, MaxAlgaeDensity) {
  load("data/parameters.RData")
  if(is.null(data$accrual)){
    accrual <- 120
  } else
  {accrual <- data$accrual}
  accrualadj <- 10^(parameters["Biggs_Coefficient1"]*(log10(accrual)+parameters["Biggs Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(accrual)^2+parameters["Biggs_Coefficient4"]))
  return(MaxAlgaeDensity * accrualadj)
}

BenthicChlora_accrual <- function (data, BenthicChlora) {
  load("data/parameters.RData")
  if(is.null(data$accrual)){
    accrual <- 120
  } else
  {accrual <- data$accrual}
  accrualadj <- 10^(parameters["Biggs_Coefficient1"]*(log10(accrual)+parameters["Biggs Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(accrual)^2+parameters["Biggs_Coefficient4"]))
  return(BenthicChlora * accrualadj)
}







