paramFrame <- read.csv("~/FlatheadPhotos/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)

dateLabelSequence <- seq(as.Date("2018-03-01"), as.Date("2020-02-01"), by = "month")
dateLabels <- dateLabelSequence[seq(1, length(dateLabelSequence), 3)]

###Mini par data
selectedPAR <- read.csv("~/FlatheadPPs/Tron/PAR_and_Depth.csv", 
                        stringsAsFactors = FALSE)
colnames(selectedPAR)[1] <- "Date"
selectedPAR$Date <- as.Date(selectedPAR$Date, format = "%m/%d/%Y")
selectedPAR <- selectedPAR[, -which(names(selectedPAR) %in% c("Site"))]
par_5m <- selectedPAR[selectedPAR$Depth_m == 5,]
par_5m <- par_5m[, -which(names(par_5m) %in% c("Depth_m"))]
colnames(par_5m)[2] <- "PAR_5m"
par_DCM <- selectedPAR[selectedPAR$Depth_m != 5,]
par_DCM <- par_DCM[, -which(names(par_DCM) %in% c("Depth_m"))]
colnames(par_DCM)[2] <- "PAR_DCM"
newParamdf <- merge(paramFrame, par_5m, by = "Date", all = TRUE)
newParamdf <- merge(newParamdf, par_DCM, by = "Date", all = TRUE)

cex.small = 0.7

tiff("FigBin/fullHierTronImPARBoth.tiff", width = 3.5, height = 4, pointsize = 12, 
     units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(0.19,0.99,0.65,0.98),las = 1, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IM_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 600),
     pch = 1, col = "#DC3220", cex = cex.small)
arrows(newParamdf$Date, newParamdf$IM_5m+newParamdf$IM_Sd_5m, newParamdf$Date, 
       newParamdf$IM_5m-newParamdf$IM_Sd_5m, length=0.03, angle=90, code=3, 
       col = "#DC3220")
points(newParamdf$Date, newParamdf$IM_5m, pch = 1, col = "#DC3220", cex = cex.small)
points(newParamdf$Date, newParamdf$PAR_5m, pch = 2, col = "#005AB5", cex = cex.small)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75")
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75")
axis(1, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, tck = -0.05, hadj = 0.6,
     at = seq(0, 600, by = 200), labels = NA)
mtext(seq(0, 600, by = 200), side = 2, padj = seq(5.5, -4.5, length.out = 4),line = 0.4)
text(as.Date("2018-04-01"), 550, "5 m", adj = 0)

par(new = "TRUE",plt = c(0.19,0.99,0.27,0.6),las = 1, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IM_DCM, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 600),
     pch = 1, col = "#DC3220", cex = cex.small)
arrows(newParamdf$Date, newParamdf$IM_DCM+newParamdf$IM_Sd_DCM, newParamdf$Date, 
       newParamdf$IM_DCM-newParamdf$IM_Sd_DCM, length=0.03, angle=90, code=3, 
       col = "#DC3220")
points(newParamdf$Date, newParamdf$IM_DCM, pch = 1, col = "#DC3220", cex = cex.small)
points(newParamdf$Date, newParamdf$PAR_DCM, pch = 2, col = "#005AB5", cex = cex.small)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75")
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75")
axis(1, tck = -0.05, hadj = 0.8,
     at = dateLabels, labels = format(dateLabels, "%b"), las = 2)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, tck = -0.05, hadj = 0.6,
     at = seq(0, 600, by = 200), labels = NA)
mtext(seq(0, 600, by = 200), side = 2, padj = seq(5.5, -4.5, length.out = 4),line = 0.4)
mtext("Date", 1, line = 3)
mtext("2018", 1, line = 2.3, adj = 0.15)
mtext("2019", 1, line = 2.3, adj = 0.7)
mtext(expression("Light (μmol photons m"^{"-2"}~" s"^{"-1"}~")"), 2, line = 2, 
      las = 0, adj = -0.5)
text(as.Date("2018-04-01"), 550, "Chl max", adj = 0)

par(xpd = TRUE)
legend("bottom", legend = c(expression("E"["m"]), 
                            expression(paste(italic("in situ"), " PAR"))), 
       horiz = TRUE, pch = c(1, 2), col = c("#DC3220", "#005AB5"), lty = 0,
       inset = c(0, -0.85), bty = "n")
par(xpd = FALSE)
dev.off()


miniPAR <- newParamdf[seq(1, 27), c("Date", "PAR_5m", "PAR_DCM")]
miniPAR$Month <- as.numeric(as.character(miniPAR$Date, "%m"))
miniPAR$Date <- NULL
m5Safe <- miniPAR[, c("PAR_5m", "Month")]

savedDatm5 <- do.call(data.frame, aggregate(.~Month, m5Safe, 
                                            FUN = function(x) c(mean(x), sd(x))))
savedDat <- do.call(data.frame, aggregate(.~Month, miniPAR, 
                                            FUN = function(x) c(mean(x), sd(x))))
