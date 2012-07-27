tn <- c(0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7, 1.8,1.9, 2.0)
load("data/benthic.RData")
load("data/density.RData")
data <- cbind(tn, benthic, density)


colnames(data) <- c("tn", "Dodds '97, mean benthic chl a", "Dodds '97, max benthic chl a", "Dodds '02/'06, mean benthic chl a",
                    "Dodds '02/'06, max benthic chl a", "Dodds '97, mean algal density", "Dodds '97, max algal density", 
                    "Dodds '02/'06, mean algal density", "Dodds '02/'06, max algal density"
)
d <- data.frame(cbind(rep(tn, 8), unlist(c(data[, 2:9]))))
d$names <- rep(colnames(data)[2:9], each=20)
colnames(d) <- c("tn", "P", "names")
d$type <- rep(c("Benthic Chlorophyll a", "Algal Density"), each=80)
d$model <- rep(c("Dodds '97", "Dodds '02/'06"), each=40, times=2)
d$parameter <- rep(c("mean Chl a", "max Chl a"), each=20, times=4)

benthicChlor_dodds <- d[d$type=="Benthic Chlorophyll a",]

library(ggplot2)
benthic_dodds <- ggplot(benthicChlor_dodds , aes(x=tn, y=P, color=names))+geom_line()+ 
  coord_cartesian(xlim = c(0, 2),ylim = c(0, .35)) + 
  scale_colour_brewer(palette="Set1", name="Model") +
  scale_x_discrete(name="Total Nitrogen") + scale_y_continuous(name="Total Phosporous") +
  facet_wrap(model ~ parameter, scales="free")


algalDensity_dodds <- d[d$type=="Algal Density",]
  
algaldensity_dodds <- ggplot(algalDensity_dodds, aes(x=tn, y=P, group=names, color=names))+geom_line()+ 
  coord_cartesian(xlim = c(0, 2),ylim = c(0, 25)) + 
  scale_colour_brewer(palette="Set1", name="Model") +
  scale_x_discrete(name="Total Nitrogen") + scale_y_continuous(name="Total Phosporous")  +
  facet_wrap(model ~ parameter, scales="free")

pdf(file="test.pdf")
benthic_dodds
algaldensity_dodds
dev.off()