allModeledFull <- read.csv("allValuesModeledFull_m3_withChl_5mOnly.csv", 
                           stringsAsFactors = FALSE)

allModeledFull$Date <- as.Date(allModeledFull$Date)

allDates <- seq.Date(min(allModeledFull$Date), max(allModeledFull$Date), by = "day")


##Daily PAR
YBPoint <- read.csv("~/Meetings/ASLO23/YBPoints_PAR.csv", stringsAsFactors = FALSE)
YBPoint$DateTime <- as.POSIXlt(YBPoint$timestamp, format = "%m/%d/%Y %H:%M")
YBPoint$id <- NULL
YBPoint$parameterID <- NULL
YBPoint$parameterName <- NULL
YBPoint$parameterUnits <- NULL #µmol photons/s/m^2
YBPoint$timestamp <- NULL
YBPoint$Date <- as.Date.POSIXlt(YBPoint$DateTime)
##Per sec, to per 5 min (*60*5/1000) (mmol/5 min)
YBPoint$parameterValue <- YBPoint$parameterValue *60*5/1000000

dailySumPAR <- aggregate(parameterValue ~ as.Date(YBPoint$Date), 
                         YBPoint, FUN = sum)
colnames(dailySumPAR) <- c("Date", "PAR")

calcNPP <- read.csv("~/FlatheadPublic/NPP/calculatedNPP_Jul2023.csv", 
                    stringsAsFactors = FALSE)
calcNPP$X <- NULL
calcNPP$Date <- as.Date(calcNPP$Date)

library(aqeco)

trim14 <- calcNPP[which(calcNPP$Date %in% allDates),  
                  c("Date", "Depth", "NPP")]
subdfMean <- aggregate(.~Date+Depth, trim14, FUN = mean)
subdfMean <- subdfMean[subdfMean$Depth == 5,]

subdfMean <- merge(subdfMean, allModeledFull, all = TRUE)
subdfMean$X <- NULL
subdfMean$Depth <- NULL
subdfMean$NPP <- subdfMean$NPP*12.01

monthSeq <- seq(as.Date("2018-01-01"), as.Date("2020-01-01"), by = "month")

trimSub <- subdfMean[!is.na(subdfMean$NPP),]

#v <- prcomp(cbind(trimSub$NPP,trimSub$mgCm.3))$rotation
#beta <- v[2,1]/v[1,1]
#interc <- mean(trimSub$mgCm.3) - beta*mean(trimSub$NPP)

library(smatr)
com.test <- sma(trimSub$mgCm.3~trimSub$NPP)
intercSMA <- com.test$coef[[1]]$`coef(SMA)`[1]
slopeSMA <- com.test$coef[[1]]$`coef(SMA)`[2]
r2SMA <- com.test$r2[[1]]

col2 <- "#D55E00"

cexSmall <- 0.8
tiff("FigBin/ModeledvsInSitu_5mTS_noPAR_11_withChl.tiff", width = 7, height = 5, 
     pointsize = 12, units = "in", res = 1200)

plot.new()
par(new = "TRUE",plt = c(0.125, 0.84, 0.6, 0.99), xpd = FALSE)
plot(subdfMean$Date, subdfMean$NPP, main = "", xaxt = "n", pch = 15,
     yaxt = "n", xlab = "", ylab = "", ylim = c(0, 42), 
     xlim = c(as.Date("2018-03-01"), as.Date("2020-01-01")))
points(subdfMean$Date, subdfMean$mgCm.3, pch = 24, bg = "white")
axis(2, at = seq(0, 40, by = 10), las = 1, hadj = 0.8)
axis(1, at = monthSeq, labels = NA, las = 1, tck = -0.03)
axis(1, at = monthSeq[seq(3, length(monthSeq), by = 3)], 
     labels = format(monthSeq[seq(3, length(monthSeq), by = 3)], "%b"), las = 1,
     tck = -0.05, padj = -1)
mtext(expression("Productivity"), 2, line = 3.2)
mtext(expression("(mg C m"^{"-3"}~"day"^{"-1"}~")"), 2, line = 2.1)
mtext("Date", 1, line = 1.5, adj = 0.48)
mtext("2018", 1, line = 1.25, adj = 0.21)
mtext("2019", 1, line = 1.25, adj = 0.75)
abline(v = as.Date("2019-01-01"), lty = 2, col = "gray")
abline(v = as.Date("2020-01-01"), lty = 2, col = "gray")
legend(as.Date("2020-01-20"), 30, xpd = NA, 
       legend = c("Measured", "Modeled"),
       #col = c("#E66100", "#5D3A9B"), 
       pch = c(15, 2), bg = "transparent", box.lwd = 0)

par(new = "TRUE",plt = c(0.25, 0.25+(5/7*0.4)*(3.5/2.5), 0.09, 0.49), xpd = FALSE)
plot(subdfMean$NPP, subdfMean$mgCm.3, main = "", xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", ylim = c(0, 30.025), xlim = c(0, 42.038))
abline(0, 1, lty = 2)
abline(intercSMA, slopeSMA, lty = 4, col = col2)
axis(1, at = seq(0, 40, by = 10), las = 1, padj = -1, tck = -0.05)
axis(1, at = seq(5, 45, by = 10), labels = NA, las = 1, hadj = 0.7, tck = -0.03)
axis(2, at = seq(0, 40, by = 10), las = 1, hadj = 0.7, tck = -0.05)
axis(2, at = seq(5, 45, by = 10), labels = NA, las = 1, hadj = 0.7, tck = -0.03)
mtext(expression("Modeled Productivity"), 2, line = 3.2)
mtext(expression("(mg C m"^{"-3"}~"day"^{"-1"}~")"), 2, line = 1.8)
mtext(expression("Measured Productivity (mg C m"^{"-3"}~"day"^{"-1"}~")"), 
      1, line = 1.5)
mtext(expression("R"^{"2"}~" = 0.55"), 4, line = 0.5, las = 1, padj = -4, col = col2)

legend(43, 20, legend = c("1:1 line", "Standardized major axis"), col = c("black", col2), 
       lty = c(2, 4), box.lwd = 0, xpd = NA, bg = "transparent")

dev.off()


