paramFrame <- read.csv("~/FlatheadPhotos/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)
paramFrame$X <- NULL

mixingRef <- read.csv("~/FlatheadPhotos/threeDepths.csv",
                      stringsAsFactors = FALSE)
mixingRef$Date <- as.Date(mixingRef$Date)
mixingRef$X <- NULL

##Three conditions:
##5 and chl max separated
##Both in epilimnion, but epilimnion shallower than photic
##Mixing below photic
colorList <- c("#1b9e77", "#d95f02", "#7570b3")
cols <- rep(colorList[2], nrow(paramFrame))
stratDates <- as.Date(mixingRef[mixingRef$DCM > mixingRef$td05, "Date"])
cols[which(paramFrame$Date %in% stratDates)] <- colorList[1]
mixedDates <- as.Date(mixingRef[mixingRef$PhoticDepth < mixingRef$td05, "Date"])
cols[which(paramFrame$Date %in% mixedDates)] <- colorList[3]

mean(paramFrame$IM_5m[which(!(paramFrame$Date %in% stratDates))])
mean(paramFrame$IM_DCM[which(!(paramFrame$Date %in% stratDates))], na.rm = TRUE)


###Stats tests
paramFrame$MixingStat <- "Strat"
paramFrame$MixingStat[which(paramFrame$Date %in% mixedDates)] <- "Mixed"

stratOnly <- paramFrame[paramFrame$MixingStat == "Strat" & !is.na(paramFrame$Alph_DCM),]
mixOnly <- paramFrame[paramFrame$MixingStat == "Mixed" & !is.na(paramFrame$Alph_DCM),]

fixedDf <- data.frame("Date" <- c(stratOnly$Date, stratOnly$Date),
                      "IM" <- c(stratOnly$IM_5m, stratOnly$IM_DCM),
                      "PBM" <- c(stratOnly$PBM_5m, stratOnly$PBM_DCM),
                      "Alph" <- c(stratOnly$Alph_5m, stratOnly$Alph_DCM),
                      "Beta" <- c(stratOnly$Bet_5m, stratOnly$Bet_DCM),
                      "Depth" <- c(rep("m5", nrow(stratOnly)), 
                                   rep("DCM", nrow(stratOnly))))
colnames(fixedDf) <- c("Date", "IM", "PBM", "Alph", "Beta", "Depth")

t.test(mixOnly$Alph_5m, mixOnly$Alph_DCM, paired = TRUE, alternative = "two.sided")

##Mixed are not normally distributed
t.test(mixOnly$Alph_5m, mixOnly$Alph_DCM, paired = FALSE)
shapiro.test(mixOnly$PBM_5m-mixOnly$PBM_DCM) #Alpha iffy
##

cellwidths <- 0.38
x0a <- 0.1
x0b <- 0.6
y0a <- 0.11
y0b <- 0.6
arrow.length <- 0.03
cex.small <- 1

tiff("FigBin/fullPaneledComp_Four_noUnits_threeMixes.tiff", width = 7, height = 7, 
     pointsize = 12, units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(paramFrame$Alph_5m, paramFrame$Alph_DCM, col = cols,
     ylim = c(0.01, 0.08), xlim = c(0.01, 0.08),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small)
arrows(paramFrame$Alph_5m, paramFrame$Alph_DCM+paramFrame$Alph_Sd_DCM, 
       paramFrame$Alph_5m, 
       paramFrame$Alph_DCM-paramFrame$Alph_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols)
arrows(paramFrame$Alph_5m+paramFrame$Alph_Sd_5m, paramFrame$Alph_DCM, 
       paramFrame$Alph_5m-paramFrame$Alph_Sd_5m, paramFrame$Alph_DCM,
       length=arrow.length, angle=90, code=3, 
       col = cols)
mtext(expression("Chl max "~alpha^{"Chl"}), 2, line = 1.95, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 2, 
#      line = 1.8, las = 0, adj = 1.1, cex = cex.small)
mtext(expression("5 m " ~ alpha^{"Chl"}), 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 1, 
#      line = 2.5, las = 0, adj = 1.1, cex = cex.small)
axis(1, tck = -0.03, padj = -1, cex.axis = cex.small,
     at = seq(0.01, 0.07, by = 0.02), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0.01, 0.07, by = 0.02), las = 2)
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(paramFrame$Bet_5m, paramFrame$Bet_DCM, 
     ylim = c(0, 0.00078), xlim = c(0, 0.00078), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small)
arrows(paramFrame$Bet_5m, paramFrame$Bet_DCM+paramFrame$Bet_Sd_DCM, paramFrame$Bet_5m, 
       paramFrame$Bet_DCM-paramFrame$Bet_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols)
arrows(paramFrame$Bet_5m+paramFrame$Bet_Sd_5m, paramFrame$Bet_DCM, 
       paramFrame$Bet_5m-paramFrame$Bet_Sd_5m, paramFrame$Bet_DCM,
       length=arrow.length, angle=90, code=3, 
       col = cols)
mtext(expression("Chl max "~beta^{"Chl"}), 2, line = 2.4, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 2, 
#      line = 2.4, las = 0, adj = 1.1, cex = cex.small)
mtext(expression("5 m "~beta^{"Chl"}), 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 1, 
#      line = 2.5, las = 0, adj = 1.1, cex = cex.small)
axis(1, tck = -0.03, padj = -0.75, at = seq(0.0001, 0.0007, by = 0.0002), las = 1,
     labels = c(expression("1x10"^{"-4"}), expression("3x10"^{"-4"}), 
                expression("5x10"^{"-4"}), expression("7x10"^{"-4"})), 
     cex.axis = cex.small)
axis(2, tck = -0.02, hadj = 0.7, cex.axis = cex.small, 
     at = seq(0.0001, 0.0007, by = 0.0002), las = 2,
     labels = c(expression("1x10"^{"-4"}), expression("3x10"^{"-4"}), 
                expression("5x10"^{"-4"}), expression("7x10"^{"-4"})))
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(paramFrame$PBM_5m, paramFrame$PBM_DCM, 
     ylim = c(0.5, 3), xlim = c(0.5,3), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small)
arrows(paramFrame$PBM_5m, paramFrame$PBM_DCM+paramFrame$PBM_Sd_DCM, paramFrame$PBM_5m, 
       paramFrame$PBM_DCM-paramFrame$PBM_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols)
arrows(paramFrame$PBM_5m+paramFrame$PBM_Sd_5m, paramFrame$PBM_DCM, 
       paramFrame$PBM_5m-paramFrame$PBM_Sd_5m, paramFrame$PBM_DCM,
       length=arrow.length, angle=90, code=3, 
       col = cols)
mtext(expression("Chl max P"["m"]^{"Chl"}), 2, line = 1.9, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}), 2, 
#      line = 1.8, las = 0, cex = cex.small)
mtext(expression("5 m P"["m"]^{"Chl"}), 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}), 1, 
#      line = 2.5, las = 0, cex = cex.small)
axis(1, tck = -0.03, padj = -1.4, cex.axis = cex.small,
     at = seq(0.5, 3, by = 0.5), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0.5, 3, by = 0.5), las = 2)
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(paramFrame$IM_5m, paramFrame$IM_DCM, 
     ylim = c(100, 600), xlim = c(100, 600), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small)
arrows(paramFrame$IM_5m, paramFrame$IM_DCM+paramFrame$IM_Sd_DCM, paramFrame$IM_5m, 
       paramFrame$IM_DCM-paramFrame$IM_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols)
arrows(paramFrame$IM_5m+paramFrame$IM_Sd_5m, paramFrame$IM_DCM, 
       paramFrame$IM_5m-paramFrame$IM_Sd_5m, paramFrame$IM_DCM,
       length=arrow.length, angle=90, code=3, 
       col = cols)
mtext(expression("Chl max E"["m"]), 2, line = 1.8, las = 0, 
      cex = cex.small)
#mtext(expression("μmol photons m"^{"-2"}~"s"^{"-1"}), 2, 
#      line = 1.8, las = 0, cex = cex.small)
mtext(expression("5 m E"["m"]), 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("μmol photons m"^{"-2"}~"s"^{"-1"}), 1, 
#      line = 2.5, las = 0, cex = cex.small)
axis(1, tck = -0.03, padj = -1.4, cex.axis = cex.small,
     at = seq(100, 500, by = 200), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(100, 500, by = 200), las = 2)
abline(a = 0, b = 1, lty = 2)
mtext("5 m", 1, line = 6, las = 0, adj = -0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

par(xpd = TRUE)
legend(-0.65, -0.98, horiz = TRUE,
       legend = c("Stratified", "Intermediate", "Deeply mixed"),
       fill = colorList, box.lwd = 0, bg = "transparent")
par(xpd = FALSE)

dev.off()