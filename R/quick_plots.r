
variables <- list("standardQual2k_BenthicChlora", "revisedQual2k_MaxAlgaeDensity",
               "MeanAlgalDen_Dodds97", "MeanAlgalDen_Dodds02")
pdf(file="validation_graphs.pdf")
nl <- lapply(variables, function(v)
  plot(alldata[, "observedChla.mg.m2."] ~ alldata[, v], main=v)
  )

l <- lapply(variables, function(v)plot(log10(alldata[, "observedChla.mg.m2."])
                                       ~ log10(alldata[, v]), main=paste(v, "log trans")))

plot(holdback$logChl ~ predict(TNTP_model, holdback), main="California TN + TP, log trans")
Chl <- 10^holdback$logChl
pred <- 10^predict(TNTP_model, holdback)
plot(Chl ~ pred, main="California TN + TP")
dev.off()