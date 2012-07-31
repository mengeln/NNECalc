# 
# ###Dodds Method###
# 
# tn <- c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7, 1.8,1.9, 2.0)
# load("doddsdata/benthic.Rdoddsdata")
# load("doddsdata/density.Rdoddsdata")
# doddsdata <- cbind(tn, benthic, density)
# 
# 
# colnames(doddsdata) <- c("tn", "Dodds '97, mean benthic chl a", "Dodds '97, max benthic chl a", "Dodds '02/'06, mean benthic chl a",
#                     "Dodds '02/'06, max benthic chl a", "Dodds '97, mean algal density", "Dodds '97, max algal density", 
#                     "Dodds '02/'06, mean algal density", "Dodds '02/'06, max algal density"
# )
# d <- doddsdata.frame(cbind(rep(tn, 8), unlist(c(doddsdata[, 2:9]))))
# d$names <- rep(colnames(doddsdata)[2:9], each=20)
# colnames(d) <- c("tn", "P", "names")
# d$type <- rep(c("Benthic Chlorophyll a", "Algal Density"), each=80)
# d$model <- rep(c("Dodds '97", "Dodds '02/'06"), each=40, times=2)
# d$parameter <- rep(c("mean Chl a", "max Chl a"), each=20, times=4)
# 
# benthicChlor_dodds <- d[d$type=="Benthic Chlorophyll a",]
# 
# library(ggplot2)
# benthic_dodds <- ggplot(benthicChlor_dodds , aes(x=tn, y=P, color=names))+geom_line()+ 
#   coord_cartesian(xlim = c(0, 2),ylim = c(0, .35)) + 
#   scale_colour_brewer(palette="Set1", name="Model") +
#   scale_x_discrete(name="Total Nitrogen") + scale_y_continuous(name="Total Phosporous") +
#   facet_wrap(model ~ parameter, scales="free")
# 
# 
# algalDensity_dodds <- d[d$type=="Algal Density",]
#   
# algaldensity_dodds <- ggplot(algalDensity_dodds, aes(x=tn, y=P, group=names, color=names))+geom_line()+ 
#   coord_cartesian(xlim = c(0, 2),ylim = c(0, 25)) + 
#   scale_colour_brewer(palette="Set1", name="Model") +
#   scale_x_discrete(name="Total Nitrogen") + scale_y_continuous(name="Total Phosporous")  +
#   facet_wrap(model ~ parameter, scales="free")
# 
# pdf(file="dodds_test.pdf")
# benthic_dodds
# algaldensity_dodds
# dev.off()


###Qual2k Method###

qual2kgraph <- function (data) {
  library(ggplot2)
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
    adjratio <- plotadj[j,] * 100
    if(j >= 5){
      accural <- 10^(parameters["Biggs_Coefficient1"]*(log10(data$accural)+parameters["Biggs Coefficient2"])+
        parameters["Biggs_Coefficient3"]*(log10(data$accural)^2+parameters["Biggs_Coefficient4"]))
      adjratio <- adjratio * accural
      }
    if(j < 3){
      Phi_NB <- sapply(1:length(data$Phi_lb), function(i){
        min(.999, (adjratio * (parameters["Respiration"] + parameters["Natural_Death_standard"]) /
          (parameters["Max_Growth20_standard"] * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature[i]-20)) *
          data$Phi_lb[i])))
      }
      )                 
      TP <- parameters["Inorg_P_Half_Sat"]*(Phi_NB/(1-Phi_NB)) * 
        (data$Phosphorus.as.P/data$OrthoPhosphate.as.P)
      
      orgN <- rowSums(data[, c("Ammonia.as.N", "Nitrate...Nitrite.as.N", 
                               "Nitrate.as.N", "Nitrite.as.N")], na.rm=T)
      TN <- parameters["Inorg_N_Half_Sat"]*(Phi_NB/(1-Phi_NB)) *
        (data$nitrogen/orgN)
      targetlimits <- data.frame(data$StationCode, data$LabSampleID, data$SampleDate, TP, TN)
    } else
    if(j >= 3){
      Phi_NB <- sapply(1:length(data$Phi_lb), function(i){
        min(.999, (adjratio * (parameters["Respiration"] + parameters["Natural_Death_revised"]) /
          (parameters["Max_Growth20_revised"] * (parameters["Arrhenius_Coefficient"]^(data$WaterTemperature[i]-20)) *
          data$Phi_lb[i])))
      }
      )
      TN <- TNlookup$TN[findInterval(Phi_NB, TNlookup$Phi_Nb)]
      
      TP <- TPlookup$TP[findInterval(Phi_NB, TPlookup$PhiNb)]
      
      targetlimits <- data.frame(data$StationCode, data$LabSampleID, data$SampleDate, TP, TN)
    }
    data$within_boundries <- sapply(1:length(data$StationCode), function(i){
      if((targetlimits$TP[i] <= data$Phosphorus.as.P[i]) | 
        (targetlimits$TN[i] <= data$nitrogen[i])){"Over allowable"}else{"Within allowable"}})
    print(
      ggplot(data, aes(x=nitrogen, y=Phosphorus.as.P, color=within_boundries))+geom_point() +
        scale_x_continuous(name="Total Nitrogen (mg/L)") + scale_y_continuous(name="Total Phosphorus (mg/L)") +
        scale_colour_hue(name= paste("Nutrient Limit for model:\n", model[j]))
          )
  }
}
pdf(file="C:\\Documents and Settings\\gisuser\\Desktop\\test.pdf")
qual2kgraph(results[which(results$standardQual2k_MaxAlgaeDensity>0),])
dev.off()