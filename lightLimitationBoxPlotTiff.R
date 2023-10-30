paramFrame <- read.csv("~/FlatheadPhotos/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)
paramFrame$X <- NULL
allDates <- paramFrame$Date

dateLabelSequence <- seq(as.Date("2018-03-01"), as.Date("2020-02-01"), by = "month")
dateLabels <- dateLabelSequence[seq(1, length(dateLabelSequence), 3)]

##K
kVals <- read.csv("~/FlatheadPublic/lightAttenuation_Jan2023_lm.csv", 
                  stringsAsFactors = FALSE)
kVals$Date <- as.Date(kVals$Date)
kVals$N <- NULL

##Incident light
YBHill <- read.csv("~/FlatheadPublic/Weather/FLBSHill_PAR_20180301to20200301.csv",
                    stringsAsFactors = FALSE)
YBHill$DateTime <- as.POSIXlt(YBHill$timestamp, format = "%m/%d/%Y %H:%M")
YBHill$id <- NULL
YBHill$parameterID <- NULL
YBHill$parameterName <- NULL
YBHill$parameterUnits <- NULL #µmol photons/s/m^2
YBHill$timestamp <- NULL
YBHill$Date <- as.Date.POSIXlt(YBHill$DateTime)

lightDf <- merge(kVals, YBHill, by = "Date")
lightDf <- lightDf[order(lightDf$DateTime),]

##DCM
dcmVals <- read.csv("~/FlatheadPublic/DCMData.csv", stringsAsFactors = FALSE)
dcmVals$Date <- as.Date(dcmVals$Date)
dcmVals$X <- NULL
dcmValsTrim <- dcmVals[, c("Date", "DCM")]
dcmValsTrim <- aggregate(DCM ~ Date, dcmValsTrim, FUN = mean)

allDf <- merge(lightDf, dcmValsTrim, by = "Date", all.x = TRUE)
allDf <- allDf[order(allDf$DateTime),]
allDf <- allDf[!is.na(allDf$DCM),]

allDf$PAR_5m <- allDf$parameterValue*exp(-allDf$kVal*5)
allDf[allDf$PAR_5m < 0, "PAR_5m"] <- 0
allDf$PAR_DCM <- allDf$parameterValue*exp(-allDf$kVal*allDf$DCM)
allDf[allDf$PAR_DCM < 0, "PAR_DCM"] <- 0

allDf0 <- merge(allDf, paramFrame, by = "Date", all.x = TRUE)
allDf0 <- allDf0[!is.na(allDf0$Alph_5m),]
allDf0 <- allDf0[allDf0$Date != as.Date("2018-08-13"),]
allDf0 <- allDf0[order(allDf0$DateTime),]

daysAll <- unique(allDf0$Date)
dateDiff <- as.numeric(daysAll - min(daysAll))
dcmDays <- unique(allDf0[!is.na(allDf0$Alph_DCM), "Date"])

dateLabelSequence <- seq.Date(as.Date("2018-03-01"), as.Date("2020-01-01"), by = "month")
dateLabels <- dateLabelSequence[seq(1, length(dateLabelSequence), by = 3)]

plattFunc <- function(par, alph, bet, pbs){
  pb <- pbs*(1-exp(-alph*par/pbs))*(exp(-bet*par/pbs))
  return(pb)
}

dois <- as.Date(c("2019-01-09", "2019-06-11"))
subdf1 <- allDf0[allDf0$Date == dois[1],]
subdf2 <- allDf0[allDf0$Date == dois[2],]

xl0 <- 0.09
xl1 <- 0.43
xGap1 <- 0.14
xGap2 <- 0.06
xrGap <- 0.01
xrWidth <- (1-xl1-xrGap-xGap1-xGap2)/2

tl0 <- 0.62
tl1 <- 0.98
bl0 <- 0.21
bl1 <- 0.57

tr0 <- 0.65
tr1 <- 0.92
br0 <- 0.3
br1 <- 0.57

timeSeq0 <- seq.POSIXt(as.POSIXlt("2019-01-09 00:00:00"), 
                       as.POSIXlt("2019-01-10 00:00:00"), by = "hour")
timeLabels0 <- timeSeq0[seq(1, 28, by = 6)]
timeSeq1 <- seq.POSIXt(as.POSIXlt("2019-06-11 00:00:00"), 
                       as.POSIXlt("2019-06-12 00:00:00"), by = "hour")
timeLabels1 <- timeSeq1[seq(1, 28, by = 6)]


cex.small = 0.6

tiff("FigBin/EmPARDiff_Diels.tiff", width = 7, height = 4, pointsize = 12, 
     units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(xl0,xl1,tl0,tl1),las = 1, xpd = FALSE)
boxplot(PAR_5m/IM_5m ~ Date, data = allDf0, at = daysAll, yaxt = "n",
        xaxt = "n", boxwex = 13, xlim = as.Date(c("2018-03-01", "2020-02-01")),
        ylim = c(0, 5.3), xlab = "", ylab = "", cex = cex.small)
abline(h = 1, lty = 3)
abline(v = as.Date(c("2019-01-01", "2020-01-01")), lty = 2, col = "darkgray")
axis(1, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(-1, 10, by = 2), labels = NA)
axis(2, tck = -0.05, hadj = 0,
     at = seq(-2, 10, by = 2), labels = seq(-2, 10, by = 2))
text(as.Date("2018-02-15"), 5, "5 m (a)", adj = 0)

par(new = "TRUE",plt = c(xl1+xGap1,xl1+xGap1+xrWidth,tr0,tr1),las = 1, 
    xpd = FALSE)
plot(subdf1$DateTime, 
     plattFunc(subdf1$PAR_5m, subdf1$Alph_5m, 
               subdf1$Bet_5m, subdf1$PBS_5m)/subdf1$PBM_5m,
     xlab = "", ylab = "", main = "", ylim = c(0, 1), type = "l", xaxt = "n", yaxt = "n")
points(subdf1$DateTime, plattFunc(subdf1$PAR_5m, subdf1$Alph_5m, 
                                  subdf1$Bet_5m, subdf1$PBS_5m)/subdf1$PBM_5m, 
       cex = cex.small)
axis(1, tck = -0.05, hadj = 0.8,
     at = timeLabels0, labels = NA, las = 2)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(0.25, 1.25, by = 0.5), labels = NA)
axis(2, tck = -0.05, hadj = 0.8,
     at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
abline(h = 1, lty = 3)
mtext("Jan 9, 2019", side = 3, line = 0.6)
mtext(expression("Modeled  "^{"14"}~"C primary production : P"["m"]^{"Chl"}),
      2, line = 2, las = 0, adj = 1)
mtext("c", 3, adj = 0.05, padj = 2)

par(new = "TRUE",plt = c(xl1+xGap1+xGap2+xrWidth,xl1+xGap1+xGap2+2*xrWidth,tr0,tr1),
    las = 1, xpd = FALSE)
plot(subdf2$DateTime, 
     plattFunc(subdf2$PAR_5m, subdf2$Alph_5m, 
               subdf2$Bet_5m, subdf2$PBS_5m)/subdf2$PBM_5m,
     xlab = "", ylab = "", main = "", ylim = c(0, 1), type = "l", xaxt = "n", yaxt = "n")
points(subdf2$DateTime, cex = cex.small, 
     plattFunc(subdf2$PAR_5m, subdf2$Alph_5m, 
               subdf2$Bet_5m, subdf2$PBS_5m)/subdf2$PBM_5m)
axis(1, tck = -0.05, hadj = 0.8,
     at = timeLabels1, labels = NA, las = 2)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(0.25, 1.25, by = 0.5), labels = NA)
axis(2, tck = -0.05, hadj = 0.8,
     at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
abline(h = 1, lty = 3)
mtext("Jun 11, 2019", side = 3, line = 0.6)
mtext("d", 3, adj = 0.05, padj = 2)


par(new = "TRUE",plt = c(xl0,xl1,bl0,bl1),las = 1, xpd = FALSE)
boxplot(PAR_DCM/IM_DCM ~ Date, data = allDf0, 
        at = daysAll[seq(4, length(daysAll))], 
        xaxt = "n", boxwex = 13, xlim = as.Date(c("2018-03-01", "2020-02-01")),
        yaxt = "n", ylim = c(0, 5.3), xlab = "", ylab = "", cex = cex.small)
abline(h = 1, lty = 3)
abline(v = as.Date(c("2019-01-01", "2020-01-01")), lty = 2, col = "darkgray")
axis(1, tck = -0.05, hadj = 0.8,
     at = dateLabels, labels = format(dateLabels, "%b"), las = 2)
axis(1, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(-1, 10, by = 2), labels = NA)
axis(2, tck = -0.05, hadj = 0,
     at = seq(-2, 10, by = 2), labels = seq(-2, 10, by = 2))
mtext("Date", 1, line = 3)
mtext(expression("in situ PAR : E"["m"]),
      2, line = 1.2, las = 0, adj = 2.8)
mtext("2018", 1, line = 2.3, adj = 0.2)
mtext("2019", 1, line = 2.3, adj = 0.75)
text(as.Date("2018-02-15"), 5, "Chl max (b)", adj = 0)


par(new = "TRUE",plt = c(xl1+xGap1,xl1+xGap1+xrWidth,br0,br1),las = 1, 
    xpd = FALSE)
plot(subdf1$DateTime,
     plattFunc(subdf1$PAR_DCM, subdf1$Alph_DCM, 
               subdf1$Bet_DCM, subdf1$PBS_DCM)/subdf1$PBM_DCM,
     xlab = "", ylab = "", main = "", type = "l", ylim = c(0, 1), xaxt = "n", yaxt = "n")
points(subdf1$DateTime, plattFunc(subdf1$PAR_DCM, subdf1$Alph_DCM, 
                                  subdf1$Bet_DCM, subdf1$PBS_DCM)/subdf1$PBM_DCM,
       cex = cex.small)
abline(h = 1, lty = 3)
axis(1, tck = -0.05, hadj = 0.8,
     at = timeLabels0, labels = format(timeLabels0, "%H:%M"), las = 2)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(0.25, 1.25, by = 0.5), labels = NA)
axis(2, tck = -0.05, hadj = 0.8,
     at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
mtext("e", 3, adj = 0.05, padj = 2)

par(new = "TRUE",plt = c(xl1+xGap1+xGap2+xrWidth,xl1+xGap1+xGap2+2*xrWidth,br0,br1),
    las = 1, xpd = FALSE)
plot(subdf2$DateTime,
     plattFunc(subdf2$PAR_DCM, subdf2$Alph_DCM, 
               subdf2$Bet_DCM, subdf2$PBS_DCM)/subdf2$PBM_DCM,
     xlab = "", ylab = "", main = "", type = "l", ylim = c(0, 1), xaxt = "n", yaxt = "n")
points(subdf2$DateTime, plattFunc(subdf2$PAR_DCM, subdf2$Alph_DCM, 
                                  subdf2$Bet_DCM, subdf2$PBS_DCM)/subdf2$PBM_DCM,
       cex = cex.small)
abline(h = 1, lty = 3)
axis(1, tck = -0.05, hadj = 0.8,
     at = timeLabels1, labels = format(timeLabels1, "%H:%M"), las = 2)
axis(2, tck = -0.02, hadj = 0.8, 
     at = seq(0.25, 1.25, by = 0.5), labels = NA)
axis(2, tck = -0.05, hadj = 0.8,
     at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
mtext("f", 3, adj = 0.05, padj = 2)

dev.off()