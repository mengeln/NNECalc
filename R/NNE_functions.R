####NNE functions: suite of functions for the NNE calculator###

###Functions SWAMPformat, Qual2k, Dodds, and NNEcalc are meant to be functions used by users;
###all other functions are low level components not meant for individual use###

##########################
###High Level Functions###
##########################

###Calculates all models starting from raw SWAMP data input###
NNEcalc <- function(testlabdata, testphabdata){
  data <- SWAMPformat(testlabdata, testphabdata)
  
  Qual2k <- Qual2k(data)
  Dodds <- Dodds(data)
  
  results <- merge(Qual2k, Dodds)
  return(results)
}

###Creates usable data frame from a SWAMP table and a phab data frame
SWAMPformat <- function(testlabdata, testphabdata){

  ###Format SWAMP data and merge it with phab data###
  calcdata <- convert_merge(testlabdata, testphabdata)

  ###Add in Misc data
  
  calcdata <- miscData(calcdata)
  
  calcdata
}

###Takes a data frame in wide format already merged with phab data and 
###returns Qual2k method predictions merged with all input data###
Qual2k <- function (calcdata) {
  ###Use calculators###
  
  density_qual2k <- MaxAlgaeDensity_standardQual2k(calcdata)
  standardQual2k_BenthicChlora <- BenthicChlora(density_qual2k)
  
  density_qual2krevised <- MaxAlgaeDensity_revisedQual2k(calcdata)
  RevisedQual2k_BenthicChlora <- BenthicChlora(density_qual2krevised)
  
  density_qual2kaccrual <- MaxAlgaeDensity_accrual(calcdata, density_qual2krevised)
  benthic_qual2kaccrual <- BenthicChlora_accrual(calcdata, RevisedQual2k_BenthicChlora)
  
  ###Target value calculator needs to be added here###
  
  ###Bind results###
  qual2k_results <- data.frame(density_qual2k, standardQual2k_BenthicChlora, density_qual2krevised, RevisedQual2k_BenthicChlora,
                               density_qual2kaccrual, benthic_qual2kaccrual)
  qual2k_results <- cbind(calcdata, qual2k_results)
  return(qual2k_results)
}

###Takes a data frame in wide format already merged with phab data and 
###returns Qual2k method predictions
Dodds <- function(calcdata){
  Dodds97 <- MaxAlgaeDensity_Dodds97(calcdata)
  Dodds02 <- MaxAlgaeDensity_Dodds02(calcdata)
  
  results <- cbind(calcdata, Dodds97, Dodds02)
  return(results)
}

###Rename columns for data being entered already in wide format###
colRename <- function(data){
  colnames(data)[grep("SiteCode", colnames(data))]  <- "StationCode"
  colnames(data)[grep("Nitrogen..Total", colnames(data))] <- "Nitrogen..Total"
  colnames(data)[grep("Ammoni", colnames(data))] <- "Ammonia.as.N"
  #colnames(data)[grep("Organic.N", colnames(data))] <- "nitrate...nitrite.as.N"
  #colnames(data)[grep("Nitrate/Nitrite", colnames(data))] <- "nitrate...nitrite.as.N"
  colnames(data)[grep("Nitrate", colnames(data))] <- "Nitrate.as.N"
  colnames(data)[grep("Nitrite", colnames(data))] <- "Nitrite.as.N"
  #colnames(data)[grep("nitrate...nitrite.as.N", colnames(data))] <- "Nitrate...Nitrite.as.N"
  colnames(data)[grep("Organic.N", colnames(data))] <-  "Nitrogen..Total.Kjeldahl"
  
  
  colnames(data)[grep("OrthoPhosphate", colnames(data))] <- "OrthoPhosphate.as.P"
  colnames(data)[grep("Phosphorus", colnames(data))] <- "Phosphorus.as.P"
  colnames(data)[grep("Organic.P", colnames(data))] <- "OrganicP"
  
  colnames(data)[grep("Depth", colnames(data))] <- "WaterDepth"
  colnames(data)[grep("Temperature", colnames(data))] <- "WaterTemperature"
  colnames(data)[grep("Turbidity", colnames(data))] <- "Turbidity"
  colnames(data)[grep("Canopy", colnames(data))] <- "CanopyClosure"
  colnames(data)[grep("Month", colnames(data))] <- "DWC_Month"
  
  if(!is.integer(data$DWC_Month)){
    month <- c("January", "February", "March", "April", "May", "June", "July", "August", "Sempter",
               "October", "November", "December")
    data$DWC_Month <- sapply(data$DWC_Month, function(d){which(month==as.character(d))
    })
  }
  
  data <- nitrogen(data)
  
  return(data)
}

#########################
###Low Level Functions###
#########################
###Convert SWAMP data and merge it with PHAB data###
convert_merge <- function(testlabdata, testphabdata){
  ###Error check and format lab data###
  testdataf <- NNEformat(testlabdata)

  ###Change IDs to characters###
  testphabdata$StationCode <- as.character(testphabdata$StationCode)
  testphabdata$SampleDate <- as.character(testphabdata$SampleDate)
  testlabdata$StationCode <- as.character(testlabdata$StationCode)
  testlabdata$SampleDate <- as.character(testlabdata$SampleDate)
  ###Merge data###
  calcdata <- merge(testdataf, testphabdata, by=c("StationCode", "SampleDate"),
                    all=F)
  calcdata
}

###Checks and returns data points with a zscore greater than 2###
outliercheck <- function(data){
  ###Create objects to store results from loop###
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

###Functions for formatting a SWAMP style data table; also does QA/QC###
NNEformat <- function (data) {
  ###Clean up SWAMP data table###
  
  data <- data[which(data$SampleTypeName %in% c("Grab", "Integrate")),]
  
  
  ###Fill in ND and DNQ##
  
  if(!is.null(data$ResQualCode) & !is.null(data$MDL) & !is.null(data$RL)){
    MDL <- which(data$MDL>0)
    
    ND <- intersect(which(data$ResQualCode=="ND"), MDL)
    data$Result[ND] <- data$MDL[ND]/2
    
    RL <- which(data$RL>0)
    DNQ <- intersect(which(data$ResQualCode=="DNQ"), RL)
    data$Result[DNQ] <- (data$RL[DNQ] - data$MDL[DNQ])/2
  }
  
  data$Result[which(data$Result<0)] <- NA
  data$Result[which(data$Result=="")] <- NA
  
  data$Result <- as.numeric(as.character(data$Result))
  data_aggregate <- aggregate(Result~StationCode+SampleDate+CollectionTime+
    AnalyteName, data=data, FUN=mean, na.rm=T)
  data_aggregate$Result <- round(data_aggregate$Result, digits=6)
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
  workingdata <- merge(workingdata, data[, c("CONCATENATE", "StationCode", "DWC_Month", ### merge with other useful data
                                                 "SampleDate", "CollectionTime", "Replicate")],
                       all=F, all.x=T, all.y=F)
  
  

  workingdata <- workingdata[!duplicated(workingdata[,"CONCATENATE"]),]###Remove duplicated rows

  ###Sum for total nitrogen###
  workingdata <- nitrogen(workingdata)
  
  ###Check for P errors and remove###
  Pcheck1 <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0) ###Non zero data points
  Pcheck2 <- which(workingdata$OrthoPhosphate.as.P>(workingdata$Phosphorus.as.P*1.1)) ###Orthophosphate > than total phosporous + 10%
  Pcheck <- intersect(Pcheck1, Pcheck2) ###Problem rows
  ###Replace total phosphorus with OrthoPhosphate on problem rows
  workingdata$Phosphorus.as.P[Pcheck] <- workingdata$OrthoPhosphate.as.P[Pcheck]
  
  ###Fill in Total P where Ortho is given but total P is not###
  pfill <- which(workingdata$Phosphorus.as.P==0 & workingdata$OrthoPhosphate.as.P>0)
  workingdata$Phosphorus.as.P[pfill] <- workingdata$OrthoPhosphate.as.P[pfill]
  
  ###Calculate Organic P###
  Pcheck <- which(workingdata$OrthoPhosphate.as.P>0 & workingdata$Phosphorus.as.P>0)
  workingdata$OrganicP <- rep(NA, length(workingdata[[1]]))
  workingdata$OrganicP[Pcheck] <- workingdata$Phosphorus.as.P[Pcheck] - 
    workingdata$OrthoPhosphate.as.P[Pcheck]
  workingdata$OrganicP[is.na(workingdata$OrganicP)] <- 0
  
  
  ###Check for outliers### Currently only throwing out data points outside limits; z score test 
  ##mechanism in comments
  
#   if(require(RGtk2Extras)==F)
#   {
#     install.packages("RGtl2Extras")
#     require(Gtk2Extras)
#   }
  
  workingdata <- workingdata[which(workingdata$nitrogen<100),]
  workingdata <- workingdata[which(workingdata$Phosphorus.as.P<15),]

#   outliers <- outliercheck(workingdata)
#   fix <- which(workingdata$CONCATENATE %in% outliers[, 1])
#   fixcolumns <- c("OrthoPhosphate.as.P", "Phosphorus.as.P", "Ammonia.as.N", "Nitrate...Nitrite.as.N",
#                   "Nitrate.as.N", "Nitrite.as.N", "Nitrogen..Total.Kjeldahl", "Nitrogen..Total")
#   workingdata[fix, fixcolumns] <- dfedit(workingdata[fix, fixcolumns])
#   workingdata <- workingdata[intersect(intersect(which(!is.na(workingdata$OrthoPhosphate.as.P)), 
#                                        which(!is.na(workingdata$OrganicP))), intersect(
#                                        which(!is.na(workingdata$Phosphorus.as.P)),
#                                        which(!is.na(workingdata$nitrogen)))),]

  return(workingdata)
}

miscData <- function (calcdata) {
  ###Fix column names###
  colnames(calcdata)[grep("Canopy", colnames(calcdata))] <- "CanopyClosure"
  colnames(calcdata)[grep("Depth", colnames(calcdata))] <- "WaterDepth"
  colnames(calcdata)[grep("Temperature", colnames(calcdata))] <- "WaterTemperature"
  
  calcdata <- calcdata[which(!is.na(calcdata$Latitude)),]
  calcdata <- calcdata[which(!is.na(calcdata$DWC_Month)),]
  
  calcdata$Turbidity[which(calcdata$Turbidity>300)] <- 300
  calcdata$Turbidity[which(calcdata$Turbidity<0)] <- 0
  
  ###Calc misc data###
  
  calcdata$SolarRadiation <- SolarRadiation(calcdata)
  
  calcdata$LightFactor <- LightFactor(calcdata)
  
  calcdata$LightExtinction <- LightExtinction(calcdata)
  
  calcdata$accrual <- accrual(calcdata)
  
  data.frame(calcdata)
}

###Calculate total nitrogen###
nitrogen <- function (workingdata) {
  values <- c("Ammonia.as.N", "Nitrate...Nitrite.as.N", ###Names of nitrogen columns except total
              "Nitrate.as.N", "Nitrite.as.N", "Nitrogen..Total.Kjeldahl")
  
  ###Remove observations with missing values###
  for(i in which(colnames(workingdata) %in% c(values, 
                                              "OrthoPhosphate.as.P", "Phosphorus.as.P"))){
    if(length(which(workingdata[, i]<0))>0){
      workingdata[which(workingdata[, i]<0), i] <- NA } ###Negative values replaced with NA
  }
  
  ###Sum for total nitrogen###
  workingdata$nitrogen <- apply(workingdata[,which(colnames(workingdata) %in% values)] , 1, sum) ###Sum nitrogen constituents 
  workingdata$nitrogen[which(workingdata$Nitrogen..Total>0)] <- 
    workingdata$Nitrogen..Total[which(workingdata$Nitrogen..Total>0)]  ###Replace summed total with reported total, if given
  workingdata <- workingdata[!is.na(workingdata$nitrogen),] ###Throw out missing data
  
  return(workingdata)
}

###Calculates SolarRadiation using month, latitude data with a lookup table
SolarRadiation <- function (data) {
  load("data/radiationTable.RData")
  data$DWC_Month <- as.character(data$DWC_Month)
  ###Round down##
  lbound <- Vectorize(function(L, M){
    radiationTable[which(radiationTable$Month==floor(as.numeric(L))), M]
  })
  lowerBound <- unlist(lbound(data$Latitude, data$DWC_Month))
  ###Round up###
  ubound <- Vectorize(function(L, M){
    radiationTable[which(radiationTable$Month==ceiling(as.numeric(L))), M]
  })
  upperBound <- unlist(ubound(data$Latitude, data$DWC_Month))
  ###Return result###
  result <- upperBound + ((lowerBound-upperBound)*(data$Latitude-floor(data$Latitude)))
  result
}


###Calculates light factor as a function of the canopy closure###
LightFactor <- function(data){
  if(is.null(data$LightFraction)){
    fraction <- .9 ###light fraction defaults to .9 if no value is given
  } else
  {fraction <- data$LightFraction
   fraction[is.na(fraction)] <- .9} ###Fills in NAs from user input with default value
  return(1 - fraction * (1 - (1 - data$CanopyClosure / 100) ^ 2) ^ 0.5)
}

###Calculates light extinction coefficient based on turbidity###
LightExtinction <- function(data){
  if(is.null(data$Turbidity)){
    Turbidity <- rep(.6, times=length(data$StationCode)) ###Defaults to .6 if no user input
  } else
  {Turbidity <- data$Turbidity
   Turbidity[is.na(Turbidity)] <- .6} ###Fills in NAs from user input with default value
  
  (.1 * Turbidity) + .44
}

###Creates column accrual if it doesn't yet exist###
accrual <- function(data) {
  if(is.null(data$accrual)){
    accrual <- rep(120, length(data$StationCode)) ###if accrual data not given, default to 120 days
  } else
  {accrual <- data$accrual
   data$accrual[is.na(accrual)] <- 120} ###Fills in missing user data with default value
  return(accrual)
}  


###Predicts AFDM Chl-a with standard Qual2k method###
###Function takes a data frame in wide format with columns for nutrients and phab data###
MaxAlgaeDensity_standardQual2k <- function(data){
  load("data/parameters.RData")
  ###Define Phi_Nb###
  #Nitrogen <- data$nitrogen
  Nitrogen <- apply(data[, c("Ammonia.as.N", "Nitrite.as.N", "Nitrate.as.N", 
                             "Nitrate...Nitrite.as.N")], 1, sum, na.rm=T)
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
  return(data.frame(standardQual2k_MaxAlgaeDensity, Phi_lb)) ###Phi_lb used in calculations of TN/TP limits later
}

###Benthic Chlor-a; uses results from MaxAlgaeDensity as the argument; works for both
### standard and revised methods###
BenthicChlora <- function(algae_input){
  load("data/parameters.RData")
  standardQual2k_BenthicChlora <- algae_input[,1] * parameters["C_AFDW_ratio"]
  return(standardQual2k_BenthicChlora)
}

###Predicts AFDM Chl-a with revised Qual2k method###
###Function takes a data frame in wide format with columns for nutrients and phab data###
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

###Predicts AFDM Chl-a with revised Qual2k method with accrual adjustment###
###Takes both the original data frame and results from revised MaxAlgaeDensity as arguments###
MaxAlgaeDensity_accrual <- function (data, MaxAlgaeDensity) {
  load("data/parameters.RData")
  accrualadj <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accrual)+parameters["Biggs Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(data$accrual)^2+parameters["Biggs_Coefficient4"]))
  MaxAlgaeDensity * accrualadj
}

###Revised Benthic Chlor-a with accrual adjustment; uses results from revised Benthic Chlor-a 
###and whole data frame as arguments
BenthicChlora_accrual <- function (data, BenthicChlora) {
  load("data/parameters.RData")
  accrualadj <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accrual)+parameters["Biggs Coefficient2"])+
    parameters["Biggs_Coefficient3"]*(log10(data$accrual)^2+parameters["Biggs_Coefficient4"]))
  return(BenthicChlora * accrualadj)
}

##Dodd's 1997##
MaxAlgaeDensity_Dodds97 <- function(data){
  ###Define Nutrients###
  Nitrogen <- data$nitrogen
  phos <- data$Phosphorus.as.P
  N <- Nitrogen*1000
  P <- phos*1000
  MeanBenthicChlA_Dodds97 <- 10^((-3.2236)+(2.8263*log10(N))-0.431247*(log10(N))^2+0.2564*log10(P))
  MaxBenthicChlA_Dodds97 <- 10^((-2.70217)+(2.7857*log10(N))-0.43340*(log10(N))^2+0.30568*log10(P))
  load("data/parameters.RData")
  MeanAlgalDen_Dodds97 <- MeanBenthicChlA_Dodds97/parameters["C_AFDW_ratio"]
  MaxAlgalDen_Dodds97  <- MaxBenthicChlA_Dodds97/parameters["C_AFDW_ratio"]
  return(cbind(MeanAlgalDen_Dodds97,MeanBenthicChlA_Dodds97,MaxBenthicChlA_Dodds97,MaxAlgalDen_Dodds97))
  
}

##Dodd's 2002##
MaxAlgaeDensity_Dodds02 <- function(data){
  ###Define Nutrients###
  Nitrogen <- data$nitrogen
  phos <- data$Phosphorus.as.P
  N <- Nitrogen*1000
  P <- phos*1000
  MeanBenthicChlA_Dodds02 <- 10^((-0.408)+(0.593*log10(N))+0.204*log10(P))
  MaxBenthicChlA_Dodds02 <- 10^((0.722)+(0.349*log10(N))+0.256*log10(P))
  load("data/parameters.RData")
  MeanAlgalDen_Dodds02 <- MeanBenthicChlA_Dodds02/parameters["C_AFDW_ratio"]
  MaxAlgalDen_Dodds02  <- MaxBenthicChlA_Dodds02/parameters["C_AFDW_ratio"]
  return(cbind(MeanAlgalDen_Dodds02,MeanBenthicChlA_Dodds02,MaxBenthicChlA_Dodds02,MaxAlgalDen_Dodds02))
  
}


####target calculator; in progress###
Qual2k_targets <- function (data, target) {
  load("Data/parameters.RData") 
  load("Data/plotadj.RData")
  load("Data/TNlookup.RData")
  load("Data/TPlookup.RData")
  TPlookup <- TPlookup[order(TPlookup$PhiNb),]
  
  model <- c("Standard QUAL2K, max algae density", "Standard QUAL2K, benthic chl a",
             "Revised QUAL2K, max algae density", "Revised QUAL2K, benthic chl a",
             "Revised QUAL2K with accrual adj, max algae density", 
             "Revised QUAL2K with accrual adj, benthic chl a")
  
  for(j in 1:6){
    adjratio <- plotadj[j,] * target
    if(j >= 5){
      accrual <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accrual)+parameters["Biggs Coefficient2"])+
        parameters["Biggs_Coefficient3"]*(log10(data$accrual)^2+parameters["Biggs_Coefficient4"]))
      adjratio <- adjratio * accrual
    }
    if(j < 3){
      Phi_NB <- sapply(1:length(data$Phi_lb), function(i){
        min(.999, (adjratio * (parameters["Respiration"] + parameters["Natural_Death_standard"]) /
          (parameters["Max_Growth20_standard"] * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature[i]-20)) *
          data[i, "Phi_lb"])))
      }
      )                 
      data[, paste(model[j],"_targetTP", sep="")] <- parameters["Inorg_P_Half_Sat"]*(Phi_NB/(1-Phi_NB)) * 
        (data$Phosphorus.as.P/data$OrthoPhosphate.as.P)
      
      orgN <- rowSums(data[, c("Ammonia.as.N", "Nitrate...Nitrite.as.N", 
                               "Nitrate.as.N", "Nitrite.as.N")], na.rm=T)
      data[, paste(model[j],"_targetTN", sep="")] <- parameters["Inorg_N_Half_Sat"]*(Phi_NB/(1-Phi_NB)) *
        (data$nitrogen/orgN)
    } else
      if(j >= 3){
        Phi_NB <- sapply(1:length(data$Phi_lb), function(i){
          min(.999, (adjratio * (parameters["Respiration"] + parameters["Natural_Death_revised"]) /
            (parameters["Max_Growth20_revised"] * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature[i]-20)) *
            data[i, grep(paste("Phi_lb.", round(j/1.9 - 1), sep=""), colnames(data))])))
        }
        )
        data[, paste(model[j],"_targetTN", sep="")] <- TNlookup$TN[findInterval(Phi_NB, TNlookup$Phi_Nb)]
        
        data[, paste(model[j],"_targetTP", sep="")] <- TPlookup$TP[findInterval(Phi_NB, TPlookup$PhiNb)]
      }
  }
  return(data)
}