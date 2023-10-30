library(car)

col5m <- "#B18500"
colChl <- "#0C7BDC"
cols <- c(col5m, colChl)
pch5m <- 1
pchChl <- 2
pchs <- c(pch5m, pchChl)

paramFrame <- read.csv("~/FlatheadPhotos/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)
paramFrame$X <- NULL

mixingRef <- read.csv("~/FlatheadPhotos/MLDvsDCMDepthComparison.csv",
                      stringsAsFactors = FALSE)
mixingRef$Date <- as.Date(mixingRef$Date)

mixParam <- merge(paramFrame, mixingRef, by = "Date", all = TRUE)

fullset <- read.csv('~/FlatheadPublic/Hydrolab_Dec2020.csv', na.strings="")
fullset$Depth <- round(as.numeric(as.character(fullset$Depth)))
fullset$Temp <- as.numeric(as.character(fullset$Temp))
fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep' & 
                     fullset$Depth == 1, c("Date", "Temp")]
fullset <- fullset[400:560,]
fullset$Date <- as.Date(fullset$Date, format = "%m/%d/%y")

mixParam <- merge(mixParam, fullset, by = "Date", all.x = TRUE)

flbsHill <- read.csv("~/FlatheadPPs/FLBSHill_PARTotal_full.csv", stringsAsFactors = FALSE)
flbsHill$Date <- as.Date.POSIXct(as.POSIXct(flbsHill$timestamp, 
                                            format = "%m/%d/%Y %H:%M"))
goodDates <- names(table(flbsHill$Date)[table(flbsHill$Date) == 96])
dailyPARHill <- data.frame("Date" = as.Date(unique(flbsHill$Date)))
summedStuff <- aggregate(parameterValue~Date, flbsHill, FUN = sum)
colnames(summedStuff)[2] <- "dailyPAR"

mixParam <- merge(mixParam, summedStuff, by = "Date", all.x = TRUE)
mixParam[mixParam$dailyPAR > 100000, "dailyPAR"] <- NA

# trimParam <- mixParam[!is.na(mixParam$dailyPAR),]
# 
# cor(trimParam$td05, trimParam$dailyPAR)
# 
# model <- lm(trimParam$IM_5m + trimParam$IM_DCM ~ trimParam$td05 + trimParam$dailyPAR)
# summary(model)
# 
# vif(model)

targDeps <- c("Alph", "Bet", "PBM", "IM")
modelOutput <- data.frame(Dependent = targDeps, td05Fit = NA, depthFit = NA,
                           dailyPARFit = NA, intercFit = NA, td05P = NA, depthP = NA,
                           dailyPARP = NA, intercP = NA, R2 = NA)

miniDF <- data.frame("Alph" = c(mixParam$Alph_5m, mixParam$Alph_DCM),
                     "Bet" = c(mixParam$Bet_5m, mixParam$Bet_DCM),
                     "PBM" = c(mixParam$PBM_5m, mixParam$PBM_DCM),
                     "IM" = c(mixParam$IM_5m, mixParam$IM_DCM),
                     "td05" = c(mixParam$td05, mixParam$td05),
                     "dailyPAR" = c(mixParam$dailyPAR, mixParam$dailyPAR),
                     "depth" = c(rep("5m", nrow(mixParam)), rep("DCM", nrow(mixParam))),
                     "temp1m" = c(mixParam$Temp, mixParam$Temp))
for (i in 1:length(targDeps)){
  targDep <- targDeps[i]
  miniDf <- data.frame(td05 = miniDF$td05, dailyPAR = miniDF$dailyPAR,
                       temp1m = miniDF$temp1m,
                       param = miniDF[[targDep]], depth = miniDF$depth)
  miniDf <- miniDf[complete.cases(miniDf),]
  
  fit <- lm(param ~ td05 + dailyPAR + depth, data = miniDf)
  summary(fit)
  
  #fit2 <- lm(param ~ td05 + dailyPAR + depth + temp1m, data = miniDf)
  #summary(fit2)
  #anova_model <- aov(param ~ depth, data = miniDf)
  #summary(anova_model)
  
  #ancova_model <- aov(param ~ depth + td05 + dailyPAR, data = miniDf)
  
  #view summary of model
  #Anova(ancova_model, type="III")
  
  modelOutput[modelOutput$Dependent == targDep, "intercFit"] <-
    fit$coefficients[[1]]
  modelOutput[modelOutput$Dependent == targDep, "td05Fit"] <-
    fit$coefficients[[2]]
  modelOutput[modelOutput$Dependent == targDep, "dailyPARFit"] <-
    fit$coefficients[[3]]
  modelOutput[modelOutput$Dependent == targDep, "depthFit"] <-
    fit$coefficients[[4]]
  
  modelOutput[modelOutput$Dependent == targDep, "intercP"] <-
    summary(fit)$coefficients[,4][1]
  modelOutput[modelOutput$Dependent == targDep, "td05P"] <-
    summary(fit)$coefficients[,4][2]
  modelOutput[modelOutput$Dependent == targDep, "dailyPARP"] <-
    summary(fit)$coefficients[,4][3]
  modelOutput[modelOutput$Dependent == targDep, "depthP"] <-
    summary(fit)$coefficients[,4][4]
  
  modelOutput[modelOutput $Dependent == targDep, "R2"] <- summary(fit)$r.squared
}

miniDF$Alph_Modeled = miniDF$td05 * 
  modelOutput[modelOutput$Dependent == "Alph", "td05Fit"] + 
  miniDF$dailyPAR * 
  modelOutput[modelOutput$Dependent == "Alph", "dailyPARFit"] + 
  modelOutput[modelOutput$Dependent == "Alph", "intercFit"]
miniDF$Bet_Modeled = miniDF$td05 * 
  modelOutput[modelOutput$Dependent == "Bet", "td05Fit"] + 
  miniDF$dailyPAR * 
  modelOutput[modelOutput$Dependent == "Bet", "dailyPARFit"] + 
  modelOutput[modelOutput$Dependent == "Bet", "intercFit"]
miniDF$PBM_Modeled = miniDF$td05 * 
  modelOutput[modelOutput$Dependent == "PBM", "td05Fit"] + 
  miniDF$dailyPAR * 
  modelOutput[modelOutput$Dependent == "PBM", "dailyPARFit"] + 
  modelOutput[modelOutput$Dependent == "PBM", "intercFit"]
miniDF$IM_Modeled = miniDF$td05 * 
  modelOutput[modelOutput$Dependent == "IM", "td05Fit"] + 
  miniDF$dailyPAR * 
  modelOutput[modelOutput$Dependent == "IM", "dailyPARFit"] + 
  modelOutput[modelOutput$Dependent == "IM", "intercFit"]

miniDF[miniDF$depth == "DCM", "Alph_Modeled"] <- 
  miniDF[miniDF$depth == "DCM", "Alph_Modeled"] + 
  modelOutput[modelOutput$Dependent == "Alph", "depthFit"]
miniDF[miniDF$depth == "DCM", "Bet_Modeled"] <- 
  miniDF[miniDF$depth == "DCM", "Bet_Modeled"] + 
  modelOutput[modelOutput$Dependent == "Bet", "depthFit"]
miniDF[miniDF$depth == "DCM", "PBM_Modeled"] <- 
  miniDF[miniDF$depth == "DCM", "PBM_Modeled"] + 
  modelOutput[modelOutput$Dependent == "PBM", "depthFit"]
miniDF[miniDF$depth == "DCM", "IM_Modeled"] <- 
  miniDF[miniDF$depth == "DCM", "IM_Modeled"] + 
  modelOutput[modelOutput$Dependent == "IM", "depthFit"]
  

#vif(ancova_model)

#test for homogeneity of variance
#leveneTest(param ~ depth, data = miniDf)


#mixedDates <- as.Date(mixingRef[mixingRef$DCM < mixingRef$td05, "Date"])

#mean(paramFrame$IM_5m[which(!(paramFrame$Date %in% mixedDates))])
#mean(paramFrame$IM_DCM[which(!(paramFrame$Date %in% mixedDates))], na.rm = TRUE)
#
trimParamDF <- mixParam

cellwidths <- 0.38
x0a <- 0.1
x0b <- 0.6
y0a <- 0.11
y0b <- 0.6
arrow.length <- 0.03
cex.small <- 1

xVar <- "dailyPAR"

tiff("FigBin/fullPaneledEnvirFits_light.tiff", width = 7, height = 7, pointsize = 12, 
     units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(mixParam[, xVar], mixParam$Alph_5m, ylim = c(0, 0.08),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[1])
arrows(mixParam[, xVar], mixParam$Alph_5m+mixParam$Alph_Sd_5m, mixParam[, xVar], 
       mixParam$Alph_5m-mixParam$Alph_Sd_5m, length=arrow.length, angle=90, code=3, 
       col = cols[1])
points(mixParam[, xVar], mixParam$Alph_DCM,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[2])
arrows(mixParam[, xVar], mixParam$Alph_DCM+mixParam$Alph_Sd_DCM, mixParam[, xVar], 
       mixParam$Alph_DCM-mixParam$Alph_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols[2])


mtext(expression(alpha^{"Chl"}), 2, line = 1.95, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 2, 
#      line = 1.8, las = 0, adj = 1.1, cex = cex.small)
mtext(xVar, 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 1, 
#      line = 2.5, las = 0, adj = 1.1, cex = cex.small)
axis(1, tck = -0.03, padj = -1, cex.axis = cex.small,
     at = seq(0.01, 0.07, by = 0.02), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0.01, 0.07, by = 0.02), las = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(mixParam[, xVar], mixParam$Bet_5m, ylim = c(0, 0.0008),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[1])
arrows(mixParam[, xVar], mixParam$Bet_5m+mixParam$Bet_Sd_5m, mixParam[, xVar], 
       mixParam$Bet_5m-mixParam$Bet_Sd_5m, length=arrow.length, angle=90, code=3, 
       col = cols[1])
points(mixParam[, xVar], mixParam$Bet_DCM,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[2])
arrows(mixParam[, xVar], mixParam$Bet_DCM+mixParam$Bet_Sd_DCM, mixParam[, xVar], 
       mixParam$Bet_DCM-mixParam$Bet_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols[2])
mtext(expression(beta^{"Chl"}), 2, line = 2.4, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}~
#                   "(μmol photons m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}), 2, 
#      line = 2.4, las = 0, adj = 1.1, cex = cex.small)
mtext(xVar, 1, line = 1.6, las = 0, cex = cex.small)
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


par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(mixParam[, xVar], mixParam$PBM_5m, ylim = c(0, 3),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[1])
arrows(mixParam[, xVar], mixParam$PBM_5m+mixParam$PBM_Sd_5m, mixParam[, xVar], 
       mixParam$PBM_5m-mixParam$PBM_Sd_5m, length=arrow.length, angle=90, code=3, 
       col = cols[1])
points(mixParam[, xVar], mixParam$PBM_DCM,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[2])
arrows(mixParam[, xVar], mixParam$PBM_DCM+mixParam$PBM_Sd_DCM, mixParam[, xVar], 
       mixParam$PBM_DCM-mixParam$PBM_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols[2])
mtext(expression("P"["m"]^{"Chl"}), 2, line = 1.9, las = 0, 
      cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}), 2, 
#      line = 1.8, las = 0, cex = cex.small)
mtext(xVar, 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("gC  (gChl "~italic("a")~")"^{"-1"}~"h"^{"-1"}), 1, 
#      line = 2.5, las = 0, cex = cex.small)
axis(1, tck = -0.03, padj = -1.4, cex.axis = cex.small,
     at = seq(0.5, 3, by = 0.5), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0.5, 3, by = 0.5), las = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(mixParam[, xVar], mixParam$IM_5m, ylim = c(0, 600),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[1])
arrows(mixParam[, xVar], mixParam$IM_5m+mixParam$IM_Sd_5m, mixParam[, xVar], 
       mixParam$IM_5m-mixParam$IM_Sd_5m, length=arrow.length, angle=90, code=3, 
       col = cols[1])
points(mixParam[, xVar], mixParam$IM_DCM,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = cols[2])
arrows(mixParam[, xVar], mixParam$IM_DCM+mixParam$IM_Sd_DCM, mixParam[, xVar], 
       mixParam$IM_DCM-mixParam$IM_Sd_DCM, length=arrow.length, angle=90, code=3, 
       col = cols[2])
mtext(expression("E"["m"]), 2, line = 1.8, las = 0, 
      cex = cex.small)
#mtext(expression("μmol photons m"^{"-2"}~"s"^{"-1"}), 2, 
#      line = 1.8, las = 0, cex = cex.small)
mtext(xVar, 1, line = 1.6, las = 0, cex = cex.small)
#mtext(expression("μmol photons m"^{"-2"}~"s"^{"-1"}), 1, 
#      line = 2.5, las = 0, cex = cex.small)
axis(1, tck = -0.03, padj = -1.4, cex.axis = cex.small,
     at = seq(100, 500, by = 200), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(100, 500, by = 200), las = 2)
mtext("5 m", 1, line = 6, las = 0, adj = -0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

par(xpd = TRUE)
legend("bottom", horiz = TRUE,
       legend = c(expression("5 m"), 
                  expression("Chl max")),
       fill = cols, box.lwd = 0, bg = "transparent")
par(xpd = FALSE)

dev.off()





####Modeled v fit

maxAlph <- 0.08
maxBet <- 0.0008
maxPBM <- 3
maxIM <- 600

cellwidths <- 0.38
x0a <- 0.1
x0b <- 0.6
y0a <- 0.11
y0b <- 0.6
arrow.length <- 0.03
cex.small <- 1

r2Line <- -1.7
r2Adj = 0.1

miniDF$col <- cols[1]
miniDF[miniDF$depth == "DCM", "col"] <- cols[2]

tiff("FigBin/fullPaneledEnvirFits_modelVFit.tiff", width = 7, height = 7, pointsize = 12, 
     units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(miniDF$Alph, miniDF$Alph_Modeled, ylim = c(0, maxAlph), xlim = c(0, maxAlph),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = miniDF$col)

mtext(expression("R"^{"2"}~" = 0.53"), 3, line = r2Line, las = 0, cex = cex.small,
      adj = r2Adj)
mtext(expression("Modeled "~alpha^{"Chl"}), 2, line = 1.95, las = 0, 
      cex = cex.small)
mtext(expression("Measured " ~ alpha^{"Chl"}), 1, line = 1.6, las = 0, cex = cex.small)
axis(1, tck = -0.02, padj = -1, cex.axis = cex.small,
     at = seq(0.01, 0.07, by = 0.02), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0.01, 0.07, by = 0.02), las = 2)
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0b,y0b+cellwidths),las = 1, xpd = FALSE)
plot(miniDF$Bet, miniDF$Bet_Modeled, ylim = c(0, maxBet), xlim = c(0, maxBet),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = miniDF$col)
mtext(expression("R"^{"2"}~" = 0.36"), 3, line = r2Line, las = 0, cex = cex.small,
      adj = r2Adj)
mtext(expression("Modeled "~beta^{"Chl"}), 2, line = 2.6, las = 0, 
      cex = cex.small)
mtext(expression("Measured "~beta^{"Chl"}), 1, line = 1.8, las = 0, cex = cex.small)
axis(1, tck = -0.02, padj = -0.75, at = seq(0.0001, 0.0007, by = 0.0002), las = 1,
     labels = c(expression("1x10"^{"-4"}), expression("3x10"^{"-4"}), 
                expression("5x10"^{"-4"}), expression("7x10"^{"-4"})), 
     cex.axis = cex.small)
axis(2, tck = -0.02, hadj = 0.7, cex.axis = cex.small, 
     at = seq(0.0001, 0.0007, by = 0.0002), las = 2,
     labels = c(expression("1x10"^{"-4"}), expression("3x10"^{"-4"}), 
                expression("5x10"^{"-4"}), expression("7x10"^{"-4"})))
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0a,x0a+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(miniDF$PBM, miniDF$PBM_Modeled, ylim = c(0, maxPBM), xlim = c(0, maxPBM),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = miniDF$col)
mtext(expression("R"^{"2"}~" = 0.42"), 3, line = r2Line, las = 0, cex = cex.small,
      adj = r2Adj)
mtext(expression("Modeled P"["m"]^{"Chl"}), 2, line = 1.9, las = 0, 
      cex = cex.small)
mtext(expression("Measured P"["m"]^{"Chl"}), 1, line = 1.6, las = 0, cex = cex.small)
axis(1, tck = -0.02, padj = -1.4, cex.axis = cex.small,
     at = seq(0, 3, by = 1), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(0, 3, by = 1), las = 2)
abline(a = 0, b = 1, lty = 2)


par(new = "TRUE",plt = c(x0b,x0b+cellwidths,y0a,y0a+cellwidths),las = 1, xpd = FALSE)
plot(miniDF$IM, miniDF$IM_Modeled, ylim = c(0, maxIM), xlim = c(0, maxIM),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = cex.small, col = miniDF$col)
mtext(expression("R"^{"2"}~" = 0.48"), 3, line = r2Line, las = 0, cex = cex.small,
      adj = r2Adj)
mtext(expression("Modeled E"["m"]), 2, line = 1.8, las = 0, 
      cex = cex.small)
mtext(expression("Measured E"["m"]), 1, line = 1.6, las = 0, cex = cex.small)
axis(1, tck = -0.02, padj = -1.4, cex.axis = cex.small,
     at = seq(100, 500, by = 200), las = 1)
axis(2, tck = -0.02, hadj = 0.6, cex.axis = cex.small, 
     at = seq(100, 500, by = 200), las = 2)
abline(a = 0, b = 1, lty = 2)
mtext("5 m", 1, line = 6, las = 0, adj = -0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

par(xpd = TRUE)
legend("bottom", horiz = TRUE,
       legend = c(expression("5 m"), 
                  expression("Chl max")),
       fill = cols, box.lwd = 0, bg = "transparent")
par(xpd = FALSE)

dev.off()


write.csv(modelOutput, file = "~/FlatheadPhotos/paramModelFits.csv", quote = FALSE)
