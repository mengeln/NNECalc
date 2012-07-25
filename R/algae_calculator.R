###Read in the data###
testdata <- read.csv("tblChemCanopyDepthLatitude.csv")


###Check for outliers/bad data### 

#***SampleID MUST be in column 2***

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

outliercolumns <- c(3:8, 11, 13, 14)

outliercheck(testdata, outliercolumns)

###Check for missing critical data### 

#***SampleID MUST be in column 2***

missingcheck <- function(data, columns){
  toRemove <- unlist(sapply(importantcolumns, function(i)which(is.na(data[[i]]))))
  if(length(toRemove)==0){print("No critical missing data")}else{
    print("The following samples are missing critical data")
    return(unique(as.character(data[toRemove, 2])))
    }
}

importantcolumns <- 3:8


missingcheck(cdata, importantcolumns)
  

###Define Parameters###

load("data/parameters.RData")

###Solar Radiation Estimation###


SolarRadiation <- function (data) {
  load("data/radiationTable.RData")
  data$Collection.Month <- as.character(data$Collection.Month)
  lowerBound <- sapply(1:length(data$SampleID), function(i){
    radiationTable[radiationTable$Month==floor(data$Latitude[i]), data$Collection.Month[i]]
  })
  upperBound <- sapply(1:length(data$SampleID), function(i){
    radiationTable[radiationTable$Month==ceiling(data$Latitude[i]), data$Collection.Month[i]]
  })
  return(sapply(1:length(lowerBound), function(i){
    upperBound[i] + (upperBound[i]-lowerBound[i])*(data$Latitude[i]-floor(data$Latitude[i]))
  }))
}

testdata$SolarRadiation <- SolarRadiation(testdata)

###Calculate Light Factor

LightFactor <- function(data){
  1 - data$fraction * (1 - (1 - data$CanopyClosure / 100) ^ 2) ^ 0.5
}

testdata$LightFactor <- LightFactor(testdata)

###Standard Qual2k Method###

MaxAlgaeDensity_standardQual2k <- function(data){
  ###Define Phi_Nb###
  Nitrogen <- data$Ammonium.as.N.mg.L + data$"Nitrate...Nitrite.as.N.mg.L" + data$Nitrite.as.N.mg.L
  phos <- data$OrthoPhosphate.as.P.mg.L
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

BenthicChlora <- function(MaxAlgaeDensity){
  MaxAlgaeDensity * parameters["C_AFDW_ratio"]
}


testdata$MaxAlgaeDensity_standardQual2k <- MaxAlgaeDensity_standardQual2k(testdata)

testdata$Benthic_Chlor_a_standardQual2k <- BenthicChlora(testdata$MaxAlgaeDensity)

###Revised QUAL2K Method###

MaxAlgaeDensity_revisedQual2k <- function(data){
  ###Define Phi_Nb###
  Nitrogen <- data$Nitrogen..Total.Dissolved.mg.L
  phos <- data$Phosphorus..Total.Dissolved.mg.L + data$OrthoPhosphate.as.P.mg.L
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

MaxAlgaeDensity_revisedQual2k(testdata)