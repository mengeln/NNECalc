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
  load("data/parameters.RData")
  data <- data[which(data$Phi_NB != 0), ]
  data <- data[which(data$Phi_NB != 0), ]
  TP <- parameters["Inorg_P_Half_Sat"]*(data$Phi_NB/(1-data$Phi_NB)) * 
    (data$Phosphorus.as.P/data$OrthoPhosphate.as.P)
  
  orgN <- rowSums(data[, c("Ammonia.as.N", "Nitrate...Nitrite.as.N", 
                           "Nitrate.as.N", "Nitrite.as.N")])
  TN <- parameters["Inorg_N_Half_Sat"]*(data$Phi_NB/(1-data$Phi_NB)) *
    (data$nitrogen/orgN)
  TN[is.na(TN)] <- 0
  TN[TN==Inf] <- 0
  TP[is.na(TP)] <- 0
  TP[TP==Inf] <- 0
  targetlimits <- data.frame(data$StationCode, data$LabSampleID, data$SampleDate, TP, TN)
  data$within_boundries <- sapply(1:length(data$StationCode), function(i){
    if((targetlimits$TP[i] > data$Phosphorus.as.P[i]) & (targetlimits$TN[i] > data$nitrogen[i])){T}else{F}})
  return(ggplot(data, aes(x=nitrogen, y=Phosphorus.as.P, color=within_boundries))+geom_point())

}

qual2kgraph(results)