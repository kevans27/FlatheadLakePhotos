library(viridis)

chlDat <- read.csv("~/FlatheadPhotos/Data/chlCount.csv", stringsAsFactors = FALSE)
chlDat$Date <- as.Date(chlDat$Date, format = "%m/%d/%Y")

DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
DCMDat <- DCMDat[order(DCMDat$Date),]
DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)

DCMDat <- DCMDat[seq(99, nrow(DCMDat)),]

mldData <- read.csv("~/TMNT/MLDAll.csv", stringsAsFactors = FALSE)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

datSeq <- seq(as.Date("2018-03-01"), as.Date("2020-01-31"), by = "month")

fullset <- read.csv('~/FlatheadPublic/Hydrolab_Dec2020.csv', na.strings="")
fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
units <- fullset[1,]
fullset <- fullset[-1,]

main <- "Temp"
val <- main

df <- fullset[complete.cases(fullset[,main]),]
df <- df[seq(13969, nrow(df)),]

df$Date <- as.Date(df$Date, format = "%m/%d/%y")
suppressWarnings(df[, main] <- as.numeric(as.character(df[, main])))
df <- df[(!is.na(df[, main])),]
df <- df[(!is.na(df$Depth)),]
df$Depth <- round(as.numeric(as.character(df$Depth)))
df$Month <- as.numeric(format(df$Date, "%m"))

#df <- df[df$Depth %in% c(1, 5, 10, 15, 20, 50, 90),]

tempDF <- df[!is.na(df$Time),]

tempDF <- tempDF[tempDF$Depth < 91 & tempDF$Date > as.Date("2017-12-01") & tempDF$Date < as.Date("2020-02-01"),]

fluorDF <- tempDF
fluorDF$CHLa <- as.numeric(fluorDF$CHLa)
fluorDF <- fluorDF[!is.na(fluorDF$CHLa), c("Date", "Depth", "CHLa")]

m5Fluor <- fluorDF[fluorDF$Depth == 5, c("Date", "CHLa")]
m5Chl <- chlDat[chlDat$Depth == "5", c("Date", "Chl")]

calibrateChl <- merge(m5Fluor, m5Chl, by = "Date", all = TRUE)
calibrateChl$Ratio <- calibrateChl$Chl/calibrateChl$CHLa
calibrateChl[calibrateChl$Date == as.Date("2018-01-08"), "Ratio"] <- 1.19
calibrateChl[calibrateChl$Date == as.Date("2018-03-06"), "Ratio"] <- 6.71
calibrateChl[calibrateChl$Date == as.Date("2018-03-19"), "Ratio"] <- 3.03
calibrateChl[calibrateChl$Date == as.Date("2018-09-17"), "Ratio"] <- 5.3
calibrateChl[calibrateChl$Date == as.Date("2019-05-23"), "Ratio"] <- 11.6

fluorDF$adjChl <- fluorDF$CHLa * 
  calibrateChl[match(fluorDF$Date, calibrateChl$Date), "Ratio"]

fluorDF <- fluorDF[!is.na(fluorDF$adjChl),]

library(tgp)

allDatSeq <- seq(min(fluorDF$Date), max(fluorDF$Date), by = "day")
fldFluorAllDate <- interp.loess(as.numeric(fluorDF$Date), fluorDF$Depth, fluorDF$adjChl, 
                               gridlen = c(400,400), span = 0.12)

###Fluor
png("FigBinExtras/fluorProfileSub.png", width = 1600, height = 800)
plot.new()

#####Fluor
p <- filled.contour(x = fldFluorAllDate$x, bty = "n",
                    y = fldFluorAllDate$y, z = fldFluorAllDate$z,
                    color.palette = function(x)viridis(x), nlevels = 12, 
                    plot.title={
                      mtext("",1,line=5,las=1,cex=2.5)
                      mtext("Depth (m)",2,line=5.5, las = 0,cex=2.5, font = 2)
                    },
                    plot.axes={
                      axis(1, cex.axis=2.5, labels = NA, 
                           at = as.numeric(seq(as.Date("2018-01-01"), 
                                               as.Date("2020-01-01"), 
                                               by = "month")), lwd.ticks = 3, 
                           tck = -0.015);
                      #lines(mldData$Date, mldData$td05, col = "gray", lwd = 5); 
                      points(DCMDat$Date, DCMDat$DCM, col = "#FFF8EF",
                             pch = 18, cex = 2.5); 
                      abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", 
                             lwd = 5)
                      axis(2, cex.axis=2.5, hadj = 1.2, lwd.ticks = 3, 
                           tck = -0.03, font = 2, at = seq(30, 0, by = -10),
                           labels = NA);
                      axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, 
                           padj = 1, font = 2,
                           at = datSeq[seq(1, length(datSeq), 3)],
                           labels = NA);
                      box(lwd = 3)
                    },
                    key.axes={
                      axis(4, cex.axis = 0.001) #####################
                    },
                    ylim = c(30, 0),
                    xlim = c(as.Date("2018-01-01"), as.Date("2020-01-01")),
                    zlim = c(0, 6.5)
)
#mtext("B", 3, line = 1, font = 2, adj = 0, cex = 2.5)


dev.off()



MLDPar <- read.csv("~/FlatheadPublic/NPPIncubationPARRatios_2011-2022.csv", 
                   stringsAsFactors = FALSE)

dailyTotals <- data.frame("Date" = as.Date(MLDPar$Date), 
                          "DailyPar" = MLDPar$DailySum)
dailyTotals$DailyPar <- dailyTotals$DailyPar/1000
dailyTotals <- dailyTotals[dailyTotals$Date != as.Date("2018-08-13"),]

coll14C <- read.csv("~/FlatheadPublic/NPP/calculatedNPP_Jul2023.csv", 
                    stringsAsFactors = FALSE)
coll14C$X <- NULL
mean14C <-aggregate(NPP~Date+Depth, coll14C, FUN = mean)
mean14C$Date <- as.Date(mean14C$Date)
mean14C$NPP <- mean14C$NPP*12.01

mean14C <- mean14C[mean14C$Date > as.Date("2017-01-01") & 
                     mean14C$Date < as.Date("2020-06-01"),]

m1mean14C <- mean14C[mean14C$Depth == 1,]
m1mean14C$Depth <- 0

mean14C <- rbind(mean14C, m1mean14C)
mean14C[mean14C$NPP < 0, "NPP"] <- 0

fldNPPAllDate <- interp.loess(as.numeric(mean14C$Date), mean14C$Depth, mean14C$NPP, 
                                gridlen = c(400,400), span = 0.09)

###NPP
png("FigBinExtras/NPPProfileSub.png", width = 1600, height = 800)
plot.new()

#####NPP
p <- filled.contour(x = fldNPPAllDate$x, bty = "n",
                    y = fldNPPAllDate$y, z = fldNPPAllDate$z,
                    color.palette = function(x)viridis(x), nlevels = 12, 
                    plot.title={
                      mtext("",1,line=5,las=1,cex=2.5)
                      mtext("Depth (m)",2,line=5.5, las = 0,cex=2.5, font = 2)
                    },
                    plot.axes={
                      axis(1, cex.axis=2.5, labels = NA, 
                           at = as.numeric(seq(as.Date("2018-01-01"), 
                                               as.Date("2020-01-01"), 
                                               by = "month")), lwd.ticks = 3, 
                           tck = -0.015);
                      #lines(mldData$Date, mldData$td05, col = "gray", lwd = 5); 
                      #points(DCMDat$Date, DCMDat$DCM, col = "#FFF8EF",
                      #       pch = 18, cex = 2.5); 
                      abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", 
                             lwd = 5)
                      axis(2, cex.axis=2.5, hadj = 1.2, lwd.ticks = 3, 
                           tck = -0.03, font = 2, at = seq(30, 0, by = -10),
                           labels = NA);
                      axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, 
                           padj = 1, font = 2,
                           at = datSeq[seq(1, length(datSeq), 3)],
                           labels = NA);
                      box(lwd = 3)
                    },
                    key.axes={
                      axis(4, cex.axis = 0.001) #####################
                    },
                    ylim = c(30, 0),
                    xlim = c(as.Date("2018-01-01"), as.Date("2020-01-01")),
                    zlim = c(0, 36)
)

dev.off()




library(png)
mainCol <- "#DC3220"
col2 <- "#005AB5"

cexSmall <- 0.9

tiff("FigBinExtras/Flour_NPP_Intro_Duo.tiff", width = 7, height = 6, 
     pointsize = 12, units = "in", res = 1200)

par(mar = c(1, 3, 1, 0), xpd = NA)
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source="FigBinExtras/fluorProfileSub.png"), -0.2, 1.05, 1.9, 2.1)
rasterImage(readPNG(source="FigBinExtras/NPPProfileSub.png"), -0.2, 0.04, 1.9, 1.09)
mtext("Depth (m)", 2, line = 2, adj = 0.5)
text(1.93, 2.1, expression("μg Chl L"^{"-1"}))
text(1.9, 1.07, expression("mg C m"^{"-3"}~"day"^{"-1"}))
mtext("2018", 1, line = -1, adj = 0.2, cexSmall)
mtext("2019", 1, line = -1, adj = 0.65, cexSmall)
parSeqb <- seq(0, 36, by = 9)
parAdjb <- seq(19, 0.4, length.out = length(parSeqb))
mtext(parSeqb, 2, line = -28.5, padj = parAdjb, las = 1, cexSmall, adj = 0)
parSeqt <- seq(0, 6, by = 1)
parAdjt <- seq(-2.7, -20.2, length.out = length(parSeqt))
mtext(parSeqt, 2, line = -28.5, padj = parAdjt, las = 1, cexSmall, adj = 0)
depthSeq <- seq(0, 30, by = 10)
depthAdjb <- seq(0.4, 19, length.out = length(depthSeq))
depthAdjt <- seq(-21.8, -2.7, length.out = length(depthSeq))
mtext(depthSeq, 2, line = 1, padj = depthAdjb, las = 1, cexSmall)
mtext(depthSeq, 2, line = 1, padj = depthAdjt, las = 1, cexSmall)
datSeq0 <- format(datSeq[seq(1, length(datSeq), 3)], "%b")
datAdj0l <- seq(0.03, 0.82, length.out = length(datSeq0))
mtext(datSeq0, 1, line = -2.5, adj = datAdj0l, cexSmall)
mtext("Date", 1, line = -0.2,adj = 0.48)

dev.off()