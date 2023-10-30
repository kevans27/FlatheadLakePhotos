paramFrame <- read.csv("~/FlatheadPhotos/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)

dateLabelSequence <- seq(as.Date("2018-03-01"), as.Date("2020-02-01"), by = "month")
dateLabels <- dateLabelSequence[seq(1, length(dateLabelSequence), 3)]

spaceL <- 0.15
spaceR <- 0.14
spaceT <- 0.01
spaceB <- 0.1
gapY <- 0.05
yWidth <- (1-spaceT-spaceB-gapY*2)/3

pchs <- c(1, 2)
paramLine <- 3.8

cols <- c("#B18500", "#0C7BDC")

cex.small = 1
arrow.length = 0.02

####Multi-panel trends
tiff("FigBin/fullHierTronPooled.tiff", width = 7, height = 6, 
     units = "in", res = 1200)
plot.new()

par(new = "TRUE",plt = c(spaceL,1-spaceR,1-spaceT-yWidth,1-spaceT),las = 1, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Alph_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.1),
     pch = pchs[1], cex = cex.small, col = cols[1])
arrows(paramFrame$Date, paramFrame$Alph_5m+paramFrame$Alph_Sd_5m, paramFrame$Date, 
       paramFrame$Alph_5m-paramFrame$Alph_Sd_5m, length=arrow.length, angle=90, code=3,
       col = cols[1])
points(paramFrame$Date, paramFrame$Alph_DCM, pch = pchs[2], cex = cex.small, col = cols[2])
arrows(paramFrame$Date, paramFrame$Alph_DCM+paramFrame$Alph_Sd_DCM, paramFrame$Date, 
       paramFrame$Alph_DCM-paramFrame$Alph_Sd_DCM, length=arrow.length, angle=90, code=3,
       col = cols[2])
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 1.5)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 1.5)

axis(1, tck = -0.04, padj = 1,
     at = dateLabels, labels = NA)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 1, 
     at = seq(0.0125, 0.1125, by = 0.025), labels = NA)
axis(2, tck = -0.03, hadj = 0.8, at = seq(0, 0.1, by = 0.025), 
     labels = c("0", expression("2.5x10"^{"-2"}), expression("5.0x10"^{"-2"}),
                expression("7.5x10"^{"-2"}), expression("1.0x10"^{"-1"})))
mtext(expression(alpha^{"Chl"}), 2, line = paramLine, las = 0)
#mtext(expression(alpha^{"B"}~"(gC (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol quanta m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}~")"), 2, 
#      line = 3.6, las = 0, adj = 1.1)



par(new = "TRUE",plt = c(spaceL,1-spaceR,1-spaceT-yWidth*2-gapY,1-spaceT-gapY-yWidth),
    las = 1, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Bet_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.001),
     pch = pchs[1], cex = cex.small, col = cols[1])
arrows(paramFrame$Date, paramFrame$Bet_5m+paramFrame$Bet_Sd_5m, paramFrame$Date, 
       paramFrame$Bet_5m-paramFrame$Bet_Sd_5m, length=arrow.length, angle=90, code=3, 
       col = cols[1])
points(paramFrame$Date, paramFrame$Bet_DCM, pch = pchs[2], cex = cex.small, col = cols[2])
arrows(paramFrame$Date, paramFrame$Bet_DCM+paramFrame$Bet_Sd_DCM, paramFrame$Date, 
       paramFrame$Bet_DCM-paramFrame$Bet_Sd_DCM, length=arrow.length, angle=90, code=3,
       col = cols[2])
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 1.5)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 1.5)

axis(1, tck = -0.04, padj = 1,
     at = dateLabels, labels = NA)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 1, 
     at = seq(0.000125, 0.001125, by = 0.00025), labels = NA)
axis(2, tck = -0.03, hadj = 0.8, at = seq(0, 0.001, by = 0.00025), 
     labels = c("0", expression("2.5x10"^{"-4"}), expression("5.0x10"^{"-4"}),
                expression("7.5x10"^{"-4"}), expression("1.0x10"^{"-3"})))
#mtext(expression(beta^{"B"}~"(gC (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol quanta m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}~")"), 2, 
#      line = 3.6, las = 0, adj = 1.1)
mtext(expression(beta^{"Chl"}), 2, line = paramLine, las = 0)

legend(as.Date("2020-03-01"), 0.00075, legend = c("5 m", "Chl max"), pch = pchs,
       box.lwd = 0, xpd = NA, col = cols)

par(new = "TRUE",plt = c(spaceL,1-spaceR,1-spaceT-yWidth*3-gapY*2,1-spaceT-gapY*2-yWidth*2),
    las = 1, xpd = FALSE)
plot(paramFrame$Date, paramFrame$PBM_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 3),
     pch = pchs[1], cex = cex.small, col = cols[1])
arrows(paramFrame$Date, paramFrame$PBM_5m+paramFrame$PBM_Sd_5m, paramFrame$Date, 
       paramFrame$PBM_5m-paramFrame$PBM_Sd_5m, length=arrow.length, angle=90, code=3,
       col = cols[1])
points(paramFrame$Date, paramFrame$PBM_DCM, pch = pchs[2], cex = cex.small, col = cols[2])
arrows(paramFrame$Date, paramFrame$PBM_DCM+paramFrame$PBM_Sd_DCM, paramFrame$Date, 
       paramFrame$PBM_DCM-paramFrame$PBM_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols[2])
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 1.5)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 1.5)

axis(1, tck = -0.04, padj = -1,
     at = dateLabels, labels = format(dateLabels, "%b"), las = 1)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3.5, by = 1), labels = NA)
axis(2, tck = -0.03, hadj = 0,
     at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1))
text(as.Date("2018-07-01"), -0.84, "2018", xpd = NA)
text(as.Date("2019-07-01"), -0.84, "2019", xpd = NA)
text(as.Date("2019-01-01"), -1.2, "Date", xpd = NA)
mtext(expression("P"["M"]^{"Chl"}), 2, line = paramLine, las = 0)
#mtext(expression("P"["M"]^{"B"}~"(gC (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~")"), 
#      2, line = 1, las = 0, adj = -10)

dev.off()