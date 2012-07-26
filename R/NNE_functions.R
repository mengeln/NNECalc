outliercheck <- function(data, columns){
  Result <- vector()
  AnalyteName <- vector()
  SampleID <- vector()
  ###Calculate z scores for all observations in selected columns###
  for(j in columns){
    cmean <- .colMeans(data[, j], m=length(data[, j]), n=1, na.rm=T)
    csd <- sd(data[, j], na.rm=T)
    problems <- which(sapply(1:length(data[, j]), function(i){
      (data[i, j]-cmean)/csd})>3)
    ###Concatenate the results from each loop together###
    SampleID <- c(SampleID, data[problems, 2])
    Result <- c(Result, data[problems, j])
    AnalyteName <- c(AnalyteName, rep(colnames(data)[j], length(data[problems, j])))
  }
  ###Bind columns together## 
  
  return(cbind(SampleID, AnalyteName, Result))
}



SolarRadiation <- function (data) {
  load("data/radiationTable.RData")
  data$Collection.Month <- as.character(data$Collection.Month)
  lowerBound <- sapply(1:length(data$SampleID), function(i){
    radiationTable[radiationTable$Month==floor(data$Latitude[i]), data$DWC_Month[i]]
  }).
  upperBound <- sapply(1:length(data$SampleID), function(i){
    radiationTable[radiationTable$Month==ceiling(data$Latitude[i]), data$DWC_Month[i]]
  })
  return(sapply(1:length(lowerBound), function(i){
    upperBound[i] + (upperBound[i]-lowerBound[i])*(data$Latitude[i]-floor(data$Latitude[i]))
  }))
}

LightFactor <- function(data){
  #1 - data$fraction * (1 - (1 - data$CanopyClosure / 100) ^ 2) ^ 0.5
  1 - .9 * (1 - (1 - data$CanopyClosure / 100) ^ 2) ^ 0.5
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
    (parameters["Light_Half_Sat"]/100) * (data$Depth.cm./100))
  Phi_lb <- numerator/(numerator + parameters["Light_Half_Sat"])
  
  ###Calculate metrics###  
  return((parameters["Max_Growth20_standard"]* Phi_NB * Phi_lb)  * (parameters["Arrhenius_Coefficient"]^(data$Temperature.degC-20)) / 
    (parameters["Respiration"] + parameters["Natural_Death_standard"]))
}

###Benthic Chlor a
BenthicChlora <- function(MaxAlgaeDensity){
  MaxAlgaeDensity * parameters["C_AFDW_ratio"]
}

###Revised QUAL2K Method###
MaxAlgaeDensity_revisedQual2k <- function(data){
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
    (parameters["Light_Half_Sat"]/100) * (data$Depth.cm./100))
  Phi_lb <- numerator/(numerator + parameters["Light_Half_Sat"])
  
  ###Calculate metrics###  
  return((parameters["Max_Growth20_revised"]* Phi_NB * Phi_lb)  * (parameters["Arrhenius_Coefficient"]^(data$Temperature.degC-20)) / 
    (parameters["Respiration"] + parameters["Natural_Death_revised"]))
}

###RivsedQual2K Method, Accrual Adjustment###
MaxAlgaeDensity_accrual <- function (data, MaxAlgaeDensity) {
  accrual <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accural)+parameters["Biggs_Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(data$accural)^2+parameters["Biggs_Coefficient4"]))
  return(MaxAlgaeDensity * accrual)
}

BenthicChlora <- function (data, BenthicChlora) {
  accrual <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accural)+parameters["Biggs_Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(data$accural)^2+parameters["Biggs_Coefficient4"]))
  return(BenthicChlora * accrual)
}







