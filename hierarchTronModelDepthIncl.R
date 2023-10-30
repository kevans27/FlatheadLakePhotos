##Despite the name, depth is in fact not included

library(rstan)
library(readxl)
library(ggplot2)
library(GGally)
options(mc.cores = parallel::detectCores())


#In mL
SAVol <- 0.25
sampVol <- 15

IMFun <- function(pbsVar, alphVar, betVar){
  x = pbsVar/alphVar * log((alphVar + betVar)/betVar)
  return(x)
}
PBMFun <- function(pbsVar, alphVar, betVar){
  x = pbsVar*(alphVar/(alphVar+betVar))*(betVar/(alphVar+betVar))**(betVar/alphVar)
  return(x)
}
IKFun <- function(pbmVar, alphVar){
  x = pbmVar/alphVar
  return(x)
}

umolC <- function(sampDPM, incTime, SADPM, SA, samp, DIC){
  umol_C <- ((sampDPM/samp)/as.numeric(((SADPM/SA)/DIC))*1.06)/incTime
  return(umol_C)
}
#PB <- PBS*(1-e^(-(alph*I)/PBS))*e^(-(bet*I)/PBS)

parDat <- read.csv("~/FlatheadPPs/FLBSHill_PARTotal.csv")
parDat$timestamp <- as.POSIXct(parDat$timestamp, format = "%m/%d/%Y %H:%M")
parDat$Date <- as.Date(parDat$timestamp)
parDat$parameterValue <- as.numeric(as.character(parDat$parameterValue))
dailyPAR <- aggregate(parameterValue~Date, parDat, FUN = "sum")
dailyPAR$parameterValue <- dailyPAR$parameterValue/1000

read_excel_allsheets_9 <- function(filename, tibble = FALSE){
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) 
    readxl::read_excel(filename, sheet = X, skip = 10, col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
mySheets <- read_excel_allsheets_9("~/FlatheadPPs/Tron/PvE2018FlatheadLake.xlsx")
allDates <- as.Date(names(mySheets), format = "%m_%d_%Y")

read_excel_allsheets_2 <- function(filename, tibble = FALSE){
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) 
    readxl::read_excel(filename, sheet = X, range = "H3:J5"))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
mySheetsMin <- read_excel_allsheets_2("~/FlatheadPPs/Tron/PvE2018FlatheadLake.xlsx")
allDates <- as.Date(names(mySheets), format = "%m_%d_%Y")

metaDatTab <- data.frame("Date" = allDates, "DIC_5m" = NA, "Chl_5m" = NA, 
                         "DIC_DCM" = NA, "Chl_DCM" = NA)

for(i in 1:length(allDates)){
  miniFram <- mySheetsMin[[i]]
  if (is.na(miniFram$Depth[2])){
    metaDatTab[i, "DIC_5m"] <- miniFram[1,2]
    metaDatTab[i, "Chl_5m"] <- miniFram[1,3]
  }
  else{
    metaDatTab[i, "DIC_5m"] <- miniFram[1,2]
    metaDatTab[i, "Chl_5m"] <- miniFram[1,3]
    metaDatTab[i, "DIC_DCM"] <- miniFram[2,2]
    metaDatTab[i, "Chl_DCM"] <- miniFram[2,3]
  }
}

umolC.df <- c()
chla.df <- c()
for(i in 1:length(allDates)){
  newFram <- mySheets[[i]]
  newFram <- newFram[complete.cases(newFram[ , 5:8]),]
  saRow <- newFram[which(newFram$Well == 'SA'),]
  saIncTime <- saRow$`Incubation time [h]`
  miniFram <- newFram[-which(newFram$Well == 'SA'),]
  miniFram <- miniFram[-which(miniFram$`Light Intensity (uE/m2/s)` == '0'),]
  lightList <- as.numeric(miniFram[, grep("ight", colnames(miniFram))])
  nFram <- data.frame("Light" = lightList, "umolC_5m" = NA, "umolC_DCM" = NA)
  chlN <- data.frame("Light" = lightList, "chl_5m" = NA, "chl_DCM" = NA)
  for(j in grep("DPM", colnames(miniFram))){
    dpmList <- miniFram[,j]
    if (colnames(miniFram)[j] == "DPM Chla Max"){
      metaDatDIC <- metaDatTab[i, "DIC_DCM"]
      metaDatChl <- metaDatTab[i, "Chl_DCM"]
      saDPMItm <- saRow[j]
      tcol = 3
    }else{
      metaDatDIC <- metaDatTab[i, "DIC_5m"]
      metaDatChl <- metaDatTab[i, "Chl_5m"]
      saDPMItm <- saRow[j]
      tcol = 2
    }
    nFram[, tcol] <- umolC(dpmList, saIncTime, saDPMItm, SAVol, sampVol, metaDatDIC)*12.01
    chlN[, tcol] <- nFram[, tcol]/metaDatChl
  }
  if (is.null(umolC.df)){
    umolC.df <- c(list(nFram))
    chla.df <- c(list(chlN))
  }else{
    umolC.df <- c(umolC.df, list(nFram))
    chla.df <- c(chla.df, list(chlN))
  }
}
names(umolC.df) <- allDates

names(chla.df) <- allDates
df <- do.call("rbind", chla.df)
df$Date <- rownames(df)
df$Date <- as.Date(gsub("\\..*","",df$Date))
just5m <- df[, c("Light", "chl_5m", "Date")]
colnames(just5m) <- c("Light", "PhRate", "Date")
just5m$Depth <- "5m"
just5m$DepthNum <- 5
justDCM <- df[!is.na(df$chl_DCM), c("Light", "chl_DCM", "Date")]
colnames(justDCM) <- c("Light", "PhRate", "Date")
justDCM$Depth <- "DCM"
justDCM$DepthNum <- 15

chlAllDf <- rbind(just5m, justDCM)
chlAllDf$DOY <- as.numeric(format(chlAllDf$Date, "%j"))
chlAllDf$DateNum <- 
  chlAllDf$Date-min(chlAllDf$Date)+as.numeric(format(min(chlAllDf$Date), "%j"))

cDF <- do.call("cbind", chla.df)
DF_5m <- cDF[, grep("Light|5m", colnames(cDF))]
photDF_5m <- DF_5m[, grep("5m", colnames(DF_5m))]
lightDF_5m <- DF_5m[, grep("Light", colnames(DF_5m))]
DF_DCM <- cDF[, grep("Light|DCM", colnames(cDF))]
DF_DCM[, seq(1,6)] <- NULL
photDF_DCM <- DF_DCM[, grep("DCM", colnames(DF_DCM))]
lightDF_DCM <- DF_DCM[, grep("Light", colnames(DF_DCM))]
colnames(lightDF_5m) <- gsub(".Light", "", colnames(lightDF_5m))
colnames(photDF_5m) <- gsub(".chl_5m", "", colnames(photDF_5m))
colnames(photDF_DCM) <- gsub(".chl_DCM", "", colnames(photDF_DCM))
colnames(lightDF_DCM) <- gsub(".Light", "", colnames(lightDF_DCM))

allDates_5m <- as.Date(gsub("\\.1", "", colnames(lightDF_5m)))
allDates_DCM <- as.Date(gsub("\\.1", "", colnames(lightDF_DCM)))
allPAR_5m <- rep(NA, length(allDates_5m))
allPAR_DCM <- rep(NA, length(allDates_DCM))
for(i in 1:length(allDates_5m)){
  firstDat <- allDates_5m[i] - 6
  
  weeklyPARdf <- dailyPAR[dailyPAR$Date >= firstDat & dailyPAR$Date <= allDates_5m[i],]
  avgWeeklyPAR <- mean(weeklyPARdf$parameterValue)
  allPAR_5m[i] <- avgWeeklyPAR
}
for(i in 1:length(allDates_DCM)){
  firstDat <- allDates_DCM[i] - 6
  
  weeklyPARdf <- dailyPAR[dailyPAR$Date >= firstDat & dailyPAR$Date <= allDates_DCM[i],]
  avgWeeklyPAR <- mean(weeklyPARdf$parameterValue)
  allPAR_DCM[i] <- avgWeeklyPAR
}

for(i in 1){
  fullset <- read.csv('~/FlatheadPublic/hydrolab.csv', na.strings="")
  fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
  units <- fullset[1,]
  fullset <- fullset[-1,]
  tempStart <- as.Date("2016-01-01")
  tempEnd  <- as.Date("2020-01-31")
  
  
  df <- fullset[complete.cases(fullset[,"Temp"]),]
  df <- df[seq(13969, nrow(df)), c("Date", "Depth", "Temp")]
  
  df$Date <- as.Date(df$Date, format = "%m/%d/%y")
  suppressWarnings(df[, "Temp"] <- as.numeric(as.character(df[, "Temp"])))
  df <- df[(!is.na(df[, "Temp"])),]
  df <- df[(!is.na(df$Depth)),]
  df$Depth <- round(as.numeric(as.character(df$Depth)))
  df <- df[df$Depth == 1,]
  
  allTemp_5m <- df[match(allDates_5m, df$Date), "Temp"]
  allTemp_DCM <- df[match(allDates_DCM, df$Date), "Temp"]
}


lightDF_5m <- lightDF_5m[,colSums(is.na(lightDF_5m)) != nrow(lightDF_5m)]
lightDF_DCM <- lightDF_DCM[,colSums(is.na(lightDF_DCM)) != nrow(lightDF_DCM)]
photDF_5m <- photDF_5m[,colSums(is.na(photDF_5m)) != nrow(photDF_5m)]
photDF_DCM <- photDF_DCM[,colSums(is.na(photDF_DCM)) != nrow(photDF_DCM)]

lightDF_5m <- lightDF_5m[, !is.na(allTemp_5m)]
photDF_5m <- photDF_5m[, !is.na(allTemp_5m)]
allPAR_5m <- allPAR_5m[!is.na(allTemp_5m)]
allTemp_5m <- allTemp_5m[!is.na(allTemp_5m)]

lightDF_DCM <- lightDF_DCM[, !is.na(allTemp_DCM)]
photDF_DCM <- photDF_DCM[, !is.na(allTemp_DCM)]
allPAR_DCM <- allPAR_DCM[!is.na(allTemp_DCM)]
allTemp_DCM <- allTemp_DCM[!is.na(allTemp_DCM)]



sink("chlAdjMol.stan")
cat("
      data {
      int <lower=1> N; //number of data points
      int <lower=1> J; //number of groups
      int <lower=1> K; //number of samples per group
      matrix[K,J] Light; //input model matrix
      matrix[K,J] Phot; //photosynthetic response model matrix
      }
      
      parameters {
      vector <lower=0> [J] alph;
      real <lower=0> mu_alph;
      real <lower=0> sd_alph;
      vector <lower=0> [J] bet;
      real <lower=0> mu_bet;
      real <lower=0> sd_bet;
      vector <lower=0> [J] PBS;
      real <lower=0> mu_PBS;
      real <lower=0> sd_PBS;
      vector <lower=0> [J] sigm;
      }
      
      model { 
      //priors. 
      alph ~ normal(mu_alph, sd_alph);
      mu_alph ~ normal(0.03, 0.02);
      sd_alph ~ normal(0.02, 0.01);
      bet ~ normal(mu_bet, sd_bet);
      mu_bet ~ normal(0.0003, 0.0003);
      sd_bet ~ normal(0.0003, 0.0003);
      PBS ~ normal(mu_PBS, sd_PBS);
      mu_PBS ~ normal(1.5, 1);
      sd_PBS ~ normal(1, 0.5);

      //likelihood    	
      for (i in 1:J) {
        for (j in 1:K) {
            Phot[j,i] ~ normal(PBS[i]*(1-exp(-(alph[i]*Light[j,i])/PBS[i]))*(exp(-(bet[i]*Light[j,i])/PBS[i])), sigm[i]);
        }
      }
      
      }",
    fill=TRUE)
sink()

##Now to stan it
chlAdjDat_5m <- list(N=ncol(lightDF_5m)*nrow(lightDF_5m), J=ncol(lightDF_5m), 
                     K=nrow(lightDF_5m), Light=lightDF_5m, Phot=photDF_5m)
chlAdjDat_DCM <- list(N=ncol(lightDF_DCM)*nrow(lightDF_DCM), J=ncol(lightDF_DCM), 
                      K=nrow(lightDF_DCM), Light=lightDF_DCM, Phot=photDF_DCM)

save(list = c("chlAdjDat_5m", "chlAdjDat_DCM"), file = "~/FlatheadPhotos/hierarchTronModDat.RData")

chlAdj_5m <- stan(file="chlAdjMol.stan", data = chlAdjDat_5m, 
                  iter = 10000, chains = 8, control = list(adapt_delta = 0.99, 
                                                           max_treedepth = 15))

#png("~/FlatheadPhotos/FigBinExtras/testingPairs.png", width = 16000, height = 16000)
#pairs(chlAdj_5m)
#dev.off()

params_5m <- as.data.frame(extract(chlAdj_5m, permuted=FALSE))
names(params_5m) <- gsub(":", ".", names(params_5m), fixed = TRUE)
names(params_5m) <- gsub("[", ".", names(params_5m), fixed = TRUE)
names(params_5m) <- gsub("]", "", names(params_5m), fixed = TRUE)
params_5m$iter <- 1:5000

divergent <- get_sampler_params(chlAdj_5m, inc_warmup=FALSE)[[1]][,'divergent__']
params_5m$divergent <- divergent

traceplot(chlAdj_5m)
stan_diag(chlAdj_5m, information = "sample")

div_params_ncp <- params_5m[params_5m$divergent == 1,]
nondiv_params_ncp <- params_5m[params_5m$divergent == 0,]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params_ncp$chain.8.PBM.1, log(nondiv_params_ncp$chain.8.alph.1),
     xlab="theta.1", ylab="log(tau)",
     col="#8F2727", pch=16, cex=0.8)
points(div_params_ncp$chain.8.PBM.1, log(div_params_ncp$chain.8.alph.1),
       col="green", pch=16, cex=0.8)


allVals <- extract(chlAdj_5m, permuted = TRUE)

trimDf <- data.frame("Alpha" = allVals$mu_alph, "Beta" = allVals$mu_bet, 
                     "PBS" = allVals$mu_PBS, 
                     "lp__" = allVals$lp__)


png("FigBinExtras/paired5mBasic.png", width = 600, height = 600)
par(mar = c(3, 3, 3, 3))
ggpairs(trimDf)
dev.off()

alphaVal <- colMeans(allVals$alph)
alphaValSd <- apply(allVals$alph, 2, sd)
alphaMean <- mean(allVals$mu_alph)
alphaSd <- mean(allVals$sd_alph)
betaVal <- colMeans(allVals$bet)
betaValSd <- apply(allVals$bet, 2, sd)
betaMean <- mean(allVals$mu_bet)
betaSd <- mean(allVals$sd_bet)
PBSVal <- colMeans(allVals$PBS)
PBSValSd <- apply(allVals$PBS, 2, sd)
PBSMean <- mean(allVals$mu_PBS)
PBSSd <- mean(allVals$sd_PBS)
tronSigm <- mean(allVals$sigm)

for(i in 1:ncol(lightDF_5m)){
  IMAll0 <- IMFun(allVals$PBS[,i], allVals$alph[,i], allVals$bet[,i])
  IMVal0 <- mean(IMAll0)
  IMSd0 <- sd(IMAll0)
  PBMAll0 <- PBMFun(allVals$PBS[,i], allVals$alph[,i], allVals$bet[,i])
  PBMVal0 <- mean(PBMAll0)
  PBMSd0 <- sd(PBMAll0)
  IKAll0 <- IKFun(PBMAll0, allVals$alph[,i])
  IKVal0 <- mean(IKAll0)
  IKSd0 <- sd(IKAll0)
  if (i == 1){
    IMVal <- c(IMVal0)
    IMSd <- c(IMSd0)
    PBMVal <- c(PBMVal0)
    PBMSd <- c(PBMSd0)
    IKVal <- c(IKVal0)
    IKSd <- c(IKSd0)
  }
  else{
    IMVal <- append(IMVal, IMVal0)
    IMSd <- append(IMSd, IMSd0)
    PBMVal <- append(PBMVal, PBMVal0)
    PBMSd <- append(PBMSd, PBMSd0)
    IKVal <- append(IKVal, IKVal0)
    IKSd <- append(IKSd, IKSd0)
  }
}
sampleIDs <- colnames(photDF_5m)

sampDF5m <- data.frame("Alph_5m" = alphaVal,
                       "Alph_Sd_5m" = alphaValSd,
                       "Bet_5m" = betaVal, 
                       "Bet_Sd_5m" = betaValSd,
                       "PBS_5m" = PBSVal, 
                       "PBS_Sd_5m" = PBSValSd,
                       "PBM_5m" = PBMVal,
                       "PBM_Sd_5m" = PBMSd,
                       "IM_5m" = IMVal,
                       "IM_Sd_5m" = IMSd,
                       "IK_5m" = IKVal,
                       "IK_Sd_5m" = IKSd)
sampDF5m$Date <- sampleIDs



##Now to stan it
chlAdjDat_DCM <- list(N=ncol(lightDF_DCM)*nrow(lightDF_DCM), J=ncol(lightDF_DCM), 
                      K=nrow(lightDF_DCM), Light=lightDF_DCM, Phot=photDF_DCM, 
                      PAR=allPAR_DCM, Temp=allTemp_DCM)

chlAdj_DCM <- stan(file="chlAdjMol.stan", data = chlAdjDat_DCM, 
                   iter = 10000, chains = 8, control = list(adapt_delta = 0.99, 
                                                            max_treedepth = 15))


allVals <- extract(chlAdj_DCM, permuted = TRUE)

trimDf <- data.frame("Alpha" = allVals$mu_alph, "Beta" = allVals$mu_bet, 
                     "PBS" = allVals$mu_PBS, 
                     "lp__" = allVals$lp__)


png("FigBinExtras/pairedDCMBasic.png", width = 600, height = 600)
par(mar = c(3, 3, 3, 3))
ggpairs(trimDf)
dev.off()

alphaVal <- colMeans(allVals$alph)
alphaValSd <- apply(allVals$alph, 2, sd)
alphaMean <- mean(allVals$mu_alph)
alphaSd <- mean(allVals$sd_alph)
betaVal <- colMeans(allVals$bet)
betaValSd <- apply(allVals$bet, 2, sd)
betaMean <- mean(allVals$mu_bet)
betaSd <- mean(allVals$sd_bet)
PBSVal <- colMeans(allVals$PBS)
PBSValSd <- apply(allVals$PBS, 2, sd)
PBSMean <- mean(allVals$mu_PBS)
PBSSd <- mean(allVals$sd_PBS)
tronSigm <- mean(allVals$sigm)

for(i in 1:ncol(lightDF_DCM)){
  IMAll0 <- IMFun(allVals$PBS[,i], allVals$alph[,i], allVals$bet[,i])
  IMVal0 <- mean(IMAll0)
  IMSd0 <- sd(IMAll0)
  PBMAll0 <- PBMFun(allVals$PBS[,i], allVals$alph[,i], allVals$bet[,i])
  PBMVal0 <- mean(PBMAll0)
  PBMSd0 <- sd(PBMAll0)
  IKAll0 <- IKFun(PBMAll0, allVals$alph[,i])
  IKVal0 <- mean(IKAll0)
  IKSd0 <- sd(IKAll0)
  if (i == 1){
    IMVal <- c(IMVal0)
    IMSd <- c(IMSd0)
    PBMVal <- c(PBMVal0)
    PBMSd <- c(PBMSd0)
    IKVal <- c(IKVal0)
    IKSd <- c(IKSd0)
  }
  else{
    IMVal <- append(IMVal, IMVal0)
    IMSd <- append(IMSd, IMSd0)
    PBMVal <- append(PBMVal, PBMVal0)
    PBMSd <- append(PBMSd, PBMSd0)
    IKVal <- append(IKVal, IKVal0)
    IKSd <- append(IKSd, IKSd0)
  }
}
sampleIDs <- colnames(photDF_DCM)

sampDF_DCM <- data.frame("Alph_DCM" = alphaVal,
                         "Alph_Sd_DCM" = alphaValSd,
                         "Bet_DCM" = betaVal, 
                         "Bet_Sd_DCM" = betaValSd,
                         "PBS_DCM" = PBSVal, 
                         "PBS_Sd_DCM" = PBSValSd,
                         "PBM_DCM" = PBMVal,
                         "PBM_Sd_DCM" = PBMSd,
                         "IM_DCM" = IMVal,
                         "IM_Sd_DCM" = IMSd,
                         "IK_DCM" = IKVal,
                         "IK_Sd_DCM" = IKSd)

sampDF_DCM$Date <- sampleIDs

paramFrame <- merge(sampDF5m, sampDF_DCM, by = "Date", all = TRUE)

#write.csv(paramFrame, "~/FlatheadPhotos/fullHierTronFits.csv", quote = FALSE)


#write.csv(predictParamFrame, "~/FlatheadPhotos/fullHierTronFitsSeasonal.csv", 
#          quote = FALSE)







###Figures below



paramFrame <- read.csv("~/FlatheadPPs/fullHierTronFits.csv", stringsAsFactors = FALSE)
paramFrame$Date <- as.Date(paramFrame$Date)

plot(paramFrame$Date, paramFrame$Bet_5m)
allDates <- paramFrame$Date

png("fullHier5mSharedAxes.png", width = 1600, height = 1200)
par(mfrow = c(4,7))
for (i in 1:length(allDates)){
  plot(chla.df[[i]]$Light, chla.df[[i]]$chl_5m, xlab = "", ylab = "", xlim = c(0, 3000), 
       ylim = c(0, 3), main = allDates[i])
  forCurve <- seq(0, 3400, length.out = 1000)
  curveOutput <- paramFrame[i, "PBS_5m"]*
    (1-exp(-(paramFrame[i, "Alph_5m"]*forCurve)/paramFrame[i, "PBS_5m"]))*
    exp(-(paramFrame[i, "Bet_5m"]*forCurve)/paramFrame[i, "PBS_5m"])
  lines(forCurve, curveOutput, col = 'red')
}
dev.off()

png("fullHier5mUniqueAxes.png", width = 1600, height = 1200)
par(mfrow = c(4,7))
for (i in 1:length(allDates)){
  plot(chla.df[[i]]$Light, chla.df[[i]]$chl_5m, xlab = "", ylab = "", main = allDates[i])
  forCurve <- seq(0, 3400, length.out = 1000)
  curveOutput <- paramFrame[i, "PBS_5m"]*
    (1-exp(-(paramFrame[i, "Alph_5m"]*forCurve)/paramFrame[i, "PBS_5m"]))*
    exp(-(paramFrame[i, "Bet_5m"]*forCurve)/paramFrame[i, "PBS_5m"])
  lines(forCurve, curveOutput, col = 'red')
}
dev.off()


png("fullHierDCMSharedAxes.png", width = 1600, height = 1200)
par(mfrow = c(4,7))
for (i in 1:length(allDates)){
  if (anyNA(chla.df[[i]]$chl_DCM)){
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
         xaxt = "n", yaxt = "n")
  } else{
    plot(chla.df[[i]]$Light, chla.df[[i]]$chl_DCM, xlab = "", ylab = "", xlim = c(0, 3000), 
         ylim = c(0, 3), main = allDates[i])
    forCurve <- seq(0, 3400, length.out = 1000)
    curveOutput <- paramFrame[i, "PBS_DCM"]*
      (1-exp(-(paramFrame[i, "Alph_DCM"]*forCurve)/paramFrame[i, "PBS_DCM"]))*
      exp(-(paramFrame[i, "Bet_DCM"]*forCurve)/paramFrame[i, "PBS_DCM"])
    lines(forCurve, curveOutput, col = 'red')
  }
}
dev.off()

png("fullHierDCMUniqueAxes.png", width = 1600, height = 1200)
par(mfrow = c(4,7))
for (i in 1:length(allDates)){
  if (anyNA(chla.df[[i]]$chl_DCM)){
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
         xaxt = "n", yaxt = "n")
  } else{
    plot(chla.df[[i]]$Light, chla.df[[i]]$chl_DCM, xlab = "", ylab = "", main = allDates[i])
    forCurve <- seq(0, 3400, length.out = 1000)
    curveOutput <- paramFrame[i, "PBS_DCM"]*
      (1-exp(-(paramFrame[i, "Alph_DCM"]*forCurve)/paramFrame[i, "PBS_DCM"]))*
      exp(-(paramFrame[i, "Bet_DCM"]*forCurve)/paramFrame[i, "PBS_DCM"])
    lines(forCurve, curveOutput, col = 'red')
  }
}
dev.off()

png("fullHierTronSep2418_DCM_trim.png", width = 1600, height = 1200)
par(las = 1, mar = c(8, 8, 4, 1))
plot(chla.df[[10]]$Light, chla.df[[10]]$chl_DCM, main = allDates[10], xlab = "", ylab = "",
     cex = 2.5, cex.main = 2.5, yaxt = "n", xaxt = "n", lwd = 2.5, ylim = c(0, 1.5),
     xlim = c(0, 500))
box(lwd = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.015, padj = 0.8,
     at = seq(0, 3000, by = 500), labels = seq(0, 3000, by = 500), font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.01, padj = 1, 
     at = seq(250, 2750, by = 500), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.015, hadj = 1.2,
     at = seq(0, 3, by = 0.5), labels = seq(0, 3, by = 0.5), font = 2)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.01, hadj = 1.1,
     at = seq(0.25, 2.75, by = 0.5), labels = NA)
mtext(expression(bold("Irrandiance (μmol quanta m"^{"-2"}~"s"^{"-1"}~")")), 1, line = 5.5, 
      cex = 2.5, font = 2)
mtext(expression(bold(""^{"14"}~"C Assimilation (gC (gChl-a)"^{"-1"}~"h"^{"-1"}~")")), 2, 
      line = 5, cex = 2.5, font = 2, las = 0)
forCurve <- seq(0, 3400, length.out = 1000)
curveOutput <- paramFrame[10, "PBS_DCM"]*
  (1-exp(-(paramFrame[10, "Alph_DCM"]*forCurve)/paramFrame[10, "PBS_DCM"]))*
  exp(-(paramFrame[10, "Bet_DCM"]*forCurve)/paramFrame[10, "PBS_DCM"])
lines(forCurve, curveOutput, col = 'red', lwd = 2.5)
#lines(c(0, 3000), c(1.5, 1.5-3000*0.000603), col = 'purple', lwd = 2.5)
dev.off()



png("fullHierTronSep2418_DCM_annotated.png", width = 1600, height = 1200)
par(las = 1, mar = c(8, 8, 4, 1))
plot(chla.df[[11]]$Light, chla.df[[11]]$chl_DCM, main = allDates[11], xlab = "", ylab = "",
     cex = 2.5, cex.main = 2.5, yaxt = "n", xaxt = "n", lwd = 2.5, ylim = c(0, 2),
     xlim = c(0, 800))
box(lwd = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.015, padj = 0.8,
     at = seq(0, 3000, by = 250), labels = seq(0, 3000, by = 250), font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.01, padj = 1, 
     at = seq(125, 2875, by = 500), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.015, hadj = 1.2,
     at = seq(0, 3, by = 0.5), labels = seq(0, 3, by = 0.5), font = 2)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.01, hadj = 1.1,
     at = seq(0.25, 2.75, by = 0.5), labels = NA)
mtext(expression(bold("Irrandiance (μmol quanta m"^{"-2"}~"s"^{"-1"}~")")), 1, line = 5.5, 
      cex = 2.5, font = 2)
mtext(expression(bold(""^{"14"}~"C Assimilation (gC (gChl-a)"^{"-1"}~"h"^{"-1"}~")")), 2, 
      line = 5, cex = 2.5, font = 2, las = 0)
forCurve <- seq(0, 3400, length.out = 1100)
curveOutput <- paramFrame[11, "PBS_DCM"]*
  (1-exp(-(paramFrame[11, "Alph_DCM"]*forCurve)/paramFrame[11, "PBS_DCM"]))*
  exp(-(paramFrame[11, "Bet_DCM"]*forCurve)/paramFrame[11, "PBS_DCM"])
lines(forCurve, curveOutput, col = 'red', lwd = 4)
lines(c(-110, 3000), c(paramFrame[11, "PBM_DCM"], paramFrame[11, "PBM_DCM"]), 
      col = '#D55E00', lwd = 4, lty = 2)
text(x = 750, 1.95, expression(bold(P["M"]^{"B"})), col = '#D55E00', cex = 4)
lines(c(paramFrame[11, "IM_DCM"], paramFrame[11, "IM_DCM"]), c(-2,5), 
      col = '#0072B2', lwd = 4, lty = 2)
text(x = 275, 0.3, expression(bold(E["M"])), col = '#0072B2', cex = 4)
abline(a=0, b=paramFrame[11, "Alph_DCM"], c(-2,5), 
       col = '#009E73', lwd = 4, lty = 2)
text(x = 0, 1, expression(bold(alpha^{"B"})), col = '#009E73', cex = 4)
abline(a=2.05, b=-paramFrame[11, "Bet_DCM"], c(-2,5), 
       col = '#CC79A7', lwd = 4, lty = 2)
text(x = 750, 1.75, expression(bold(beta^{"B"})), col = '#CC79A7', cex = 4)
dev.off()



dateLabelSequence <- seq(as.Date("2018-03-01"), as.Date("2020-02-01"), by = "month")
dateLabels <- dateLabelSequence[seq(1, length(dateLabelSequence), 3)]


png("fullHierTronPBMboth.png", width = 600, height = 800)
plot.new()
par(new = "TRUE",plt = c(0.16,0.99,0.63,0.99),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$PBM_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 3),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$PBM_5m+paramFrame$PBM_Sd_5m, paramFrame$Date, 
       paramFrame$PBM_5m-paramFrame$PBM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1), font = 2)
text(as.Date("2018-04-01"), 0.22, "5 m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.16,0.99,0.22,0.58),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$PBM_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 3),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$PBM_DCM+paramFrame$PBM_Sd_DCM, paramFrame$Date, 
       paramFrame$PBM_DCM-paramFrame$PBM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%Y"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1), font = 2)
mtext("Date", 1, line = 11, cex = 2.5, font = 2)
mtext(expression(bold("P"["M"]^{"B"}~"(gC (gChl-")~bolditalic("a")~
                   bold(")"^{"-1"}~"m"^{"-2"}~"h"^{"-1"})), 
      2, line = 3, cex = 2.5, font = 2, adj = -2, las = 0)
text(as.Date("2018-04-01"), 0.22, "Chl max", font = 2, cex = 2.5, adj = 0)
dev.off()




png("fullHierTronAlphboth.png", width = 600, height = 800)
plot.new()
par(new = "TRUE",plt = c(0.22,0.99,0.63,0.99),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Alph_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.1),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Alph_5m+paramFrame$Alph_Sd_5m, paramFrame$Date, 
       paramFrame$Alph_5m-paramFrame$Alph_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0125, 0.1125, by = 0.025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 0.1, by = 0.025), labels = seq(0, 0.1, by = 0.025), font = 2)
text(as.Date("2018-04-01"), 0.007, "5 m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.22,0.99,0.22,0.58),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Alph_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 0.1),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Alph_DCM+paramFrame$Alph_Sd_DCM, paramFrame$Date, 
       paramFrame$Alph_DCM-paramFrame$Alph_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%Y"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0125, 0.1125, by = 0.025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 0.1, by = 0.025), labels = seq(0, 0.1, by = 0.025), font = 2)
mtext("Date", 1, line = 11, cex = 2.5, font = 2)
mtext(expression(bold(alpha^{"B"})), 
      2, line = 7, cex = 2.5, font = 2, las = 1, adj = 0.7, padj = -5)
text(as.Date("2018-04-01"), 0.007, "Chl max", font = 2, cex = 2.5, adj = 0)
dev.off()




png("fullHierTronBetboth.png", width = 600, height = 800)
plot.new()
par(new = "TRUE",plt = c(0.22,0.99,0.63,0.99),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Bet_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.001),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Bet_5m+paramFrame$Bet_Sd_5m, paramFrame$Date, 
       paramFrame$Bet_5m-paramFrame$Bet_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.000125, 0.001125, by = 0.00025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 0.001, by = 0.00025), 
     labels = c("0", "2.5e-4", "5e-4", "7.5e-4", "1e-3"), 
     font = 2)
text(as.Date("2018-04-01"), 0.00006, "5 m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.22,0.99,0.22,0.58),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Bet_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 0.001),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Bet_DCM+paramFrame$Bet_Sd_DCM, paramFrame$Date, 
       paramFrame$Bet_DCM-paramFrame$Bet_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%Y"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.000125, 0.001125, by = 0.00025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 0.001, by = 0.00025), 
     labels = c("0", "2.5e-4", "5e-4", "7.5e-4", "1e-3"), 
     font = 2)
mtext("Date", 1, line = 11, cex = 2.5, font = 2)
mtext(expression(bold(beta^{"B"})), 2, line = 7, cex = 2.5, font = 2, las = 1, 
      adj = 0.7, padj = -4.5)
text(as.Date("2018-04-01"), 0.00006, "Chl max", font = 2, cex = 2.5, adj = 0)
dev.off()








png("fullTronEKboth.png", width = 600, height = 800)
plot.new()
par(new = "TRUE",plt = c(0.18,0.99,0.63,0.99),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$IK_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 600),
     pch = 1, cex = 1, lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 600, by = 200), labels = seq(0, 600, by = 200), font = 2)
text(as.Date("2018-04-01"), 550, "5m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.18,0.99,0.22,0.58),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$IK_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 600),
     pch = 1, cex = 1, lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%Y"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 600, by = 200), labels = seq(0, 600, by = 200), font = 2)
mtext("Date", 1, line = 11, cex = 2.5, font = 2)
mtext(expression(bold("E"["K"])), 2, line = 5, cex = 2.5, font = 2, las = 0, adj = 1.1)
text(as.Date("2018-04-01"), 550, "Chl max", font = 2, cex = 2.5, adj = 0)
dev.off()



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

png("FigBin/fullHierTronImPARBoth.png", width = 600, height = 700)
plot.new()
par(new = "TRUE",plt = c(0.19,0.99,0.65,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IM_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 600),
     pch = 1, cex = 1, lwd = 2.5, col = "#DC3220")
arrows(newParamdf$Date, newParamdf$IM_5m+newParamdf$IM_Sd_5m, newParamdf$Date, 
       newParamdf$IM_5m-newParamdf$IM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$IM_5m, cex = 1, pch = 1, lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$PAR_5m, cex = 1, pch = 2, lwd = 2.5, col = "#005AB5")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 600, by = 200), labels = seq(0, 600, by = 200), font = 2)
text(as.Date("2018-04-01"), 550, "5 m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.19,0.99,0.27,0.6),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IM_DCM, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 600),
     pch = 1, cex = 1, lwd = 2.5, col = "#DC3220")
arrows(newParamdf$Date, newParamdf$IM_DCM+newParamdf$IM_Sd_DCM, newParamdf$Date, 
       newParamdf$IM_DCM-newParamdf$IM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$IM_DCM, cex = 1, pch = 1, lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$PAR_DCM, cex = 1, pch = 2, lwd = 2.5, col = "#005AB5")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%y"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 600, by = 200), labels = seq(0, 600, by = 200), font = 2)
mtext("Date", 1, line = 8, cex = 2.5, font = 2)
mtext(expression(bold("Light (μmol m"^{"-2"}~" s"^{"-1"}~")")), 2, line = 5, cex = 2.5, 
      font = 2, 
      las = 0, adj = -2)
text(as.Date("2018-04-01"), 550, "Chl max", font = 2, cex = 2.5, adj = 0)

par(xpd = TRUE)
legend("bottom", legend = c(expression(bold("E"["m"])), 
                            expression(bold(paste(bolditalic("in situ"), " PAR")))), 
       horiz = TRUE, pch = c(1, 2), col = c("#DC3220", "#005AB5"), lwd = 2.5, lty = 0,
       inset = c(0, -0.85), box.lwd = 3, cex = 2, bty = "n")
par(xpd = FALSE)
dev.off()


png("fullHierTronEKEMPARBoth.png", width = 600, height = 800)
plot.new()
par(new = "TRUE",plt = c(0.19,0.99,0.65,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IK_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 620),
     pch = 1, cex = 1, lwd = 2.5, col = "#DC3220")
arrows(newParamdf$Date, newParamdf$IK_5m+newParamdf$IK_Sd_5m, newParamdf$Date, 
       newParamdf$IK_5m-newParamdf$IK_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$PAR_5m, cex = 1, pch = 2, lwd = 2.5, col = "#005AB5")
points(newParamdf$Date, newParamdf$IM_5m, cex = 1, pch = 5, lwd = 2.5, col = "#FFB000")
arrows(newParamdf$Date, newParamdf$IM_5m+newParamdf$IM_Sd_5m, newParamdf$Date, 
       newParamdf$IM_5m-newParamdf$IM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#FFB000")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(50, 750, by = 50), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 600, by = 200), labels = seq(0, 600, by = 200), font = 2)
text(as.Date("2018-02-15"), 550, "5 m", font = 2, cex = 2.5, adj = 0)

par(new = "TRUE",plt = c(0.19,0.99,0.27,0.6),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(newParamdf$Date, newParamdf$IK_DCM, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 310),
     pch = 1, cex = 1, lwd = 2.5, col = "#DC3220")
arrows(newParamdf$Date, newParamdf$IK_DCM+newParamdf$IK_Sd_DCM, newParamdf$Date, 
       newParamdf$IK_DCM-newParamdf$IK_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
points(newParamdf$Date, newParamdf$PAR_DCM, cex = 1, pch = 2, lwd = 2.5, col = "#005AB5")
points(newParamdf$Date, newParamdf$IM_DCM, cex = 1, pch = 5, lwd = 2.5, col = "#FFB000")
arrows(newParamdf$Date, newParamdf$IM_DCM+newParamdf$IM_Sd_DCM, newParamdf$Date, 
       newParamdf$IM_DCM-newParamdf$IM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#FFB000")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%b"), font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(25, 375, by = 25), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, hadj = 1.1,
     at = seq(0, 300, by = 100), labels = seq(0, 300, by = 100), font = 2)
mtext("Date", 1, line = 9, cex = 2.5, font = 2)
mtext(expression(bold("Irradiance (μmol quanta m"^{"-2"}~" s"^{"-1"}~")")), 2, 
      line = 5, cex = 2.5, font = 2, las = 0, adj = -0.2)
text(as.Date("2018-02-15"), 275, "CM", font = 2, cex = 2.5, adj = 0)
mtext("2018", side = 1, line = 7, cex = 2.5, font = 2, adj = 0.18)
mtext("2019", side = 1, line = 7, cex = 2.5, font = 2, adj = 0.735)

par(xpd = TRUE)
legend("bottom", legend = c(expression(bold("E"["k"])), 
                            expression(bold("E"["m"])), 
                            expression(bolditalic("in situ"))), 
       horiz = TRUE, pch = c(1, 5, 2), col = c("#DC3220", "#FFB000", "#005AB5"), lwd = 2.5, 
       inset = c(0, -0.8), box.lwd = 3, cex = 2, lty = 0, x.intersp = 0)
par(xpd = FALSE)
dev.off()


#####Looking at light
parDat <- read.csv("~/FlatheadPublic/ShortedPARData.csv", stringsAsFactors = FALSE)
parDat$Date <- as.Date(parDat$Date)
parFits <- read.csv("~/FlatheadPPs/LightAttenuationFits.csv", stringsAsFactors = FALSE)
parFits$Date <- as.Date(parFits$Date)

png("fullRelKCurves.png", width = 1600, height = 1200)
par(mfrow = c(4,7))
for (i in 1:length(allDates)){
  dat <- allDates[i]
  
  subParDat <- parDat[parDat$Date == dat,]
  subParFits <- parFits[parFits$Date == dat,]
  
  plot(subParDat$AdjDeep, -subParDat$Depth, xlab = "", ylab = "", main = allDates[i])
  forCurve <- seq(0, 100, length.out = 1000)
  curveOutput <- subParFits$SurfaceLight * exp(-(subParFits$kVal * forCurve))
  lines(curveOutput, -forCurve, col = 'red')
}
dev.off()






paramFrame$Month <- as.numeric(as.character(format(paramFrame$Date, "%m")))
cols <- rep("#D41159", nrow(paramFrame))
cols[which(paramFrame$Month < 6 | paramFrame$Month > 9)] <- "#1A85FF"

png("~/FlatheadPPs/fullPBMComp.png", width = 800, height = 600)
plot.new()
par(mar = c(8, 9, 2, 2))
plot(paramFrame$PBM_5m, paramFrame$PBM_DCM, cex = 2, lwd = 2, 
     ylim = c(0.5, 3), xlim = c(0.5,3), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
mtext(expression(bold("P"["M"]^"B"~"(5 m)")), 1, line = 6, cex = 2.5, font = 2, las = 0, 
      adj = 0.5)
mtext(expression(bold("P"["M"]^"B"~"(FM)")), 2, line = 5, cex = 2.5, font = 2, las = 0, 
      adj = 0.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.5, 3, by = 0.5), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3, by = 0.5), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("bottomright", cex = 2.5,
       legend = c(expression(bold("Mixed")), expression(bold("Stratified"))),
       fill = c("#1A85FF", "#D41159"), box.lwd = 2)
box(lwd = 2)
dev.off()

png("~/FlatheadPPs/fullAlphComp.png", width = 800, height = 600)
plot.new()
par(mar = c(8, 9, 2, 2))
plot(paramFrame$Alph_5m, paramFrame$Alph_DCM, col = cols, cex = 2, lwd = 2,
     ylim = c(0.01, 0.075), xlim = c(0.01, 0.075),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
mtext(expression(bold(alpha^{"B"}~"(5 m)")), 1, line = 6, cex = 2.5, font = 2, las = 0, 
      adj = 0.5)
mtext(expression(bold(alpha^{"B"}~"(FM)")), 2, line = 5.5, cex = 2.5, font = 2, las = 0, 
      adj = 0.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.01, 0.07, by = 0.01), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.01, 0.07, by = 0.01), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("bottomright", cex = 2.5,
       legend = c(expression(bold("Mixed")), expression(bold("Stratified"))),
       fill = c("#1A85FF", "#D41159"), box.lwd = 2)
box(lwd=2)
dev.off()





png("FigBin/fullPaneledComp.png", width = 1728, height = 1200)
plot.new()
par(new = "TRUE",plt = c(0.08,0.33,0.59,0.95),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Alph_5m, paramFrame$Alph_DCM, col = cols, cex = 2, lwd = 2,
     ylim = c(0.01, 0.075), xlim = c(0.01, 0.075),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(0.015, 0.07, expression(bold(alpha^{"B"})), cex = 2.5)
mtext("Chl max", 2, line = 6, cex = 2.5, font = 2, las = 0, 
      adj = -0.2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.01, 0.07, by = 0.02), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.01, 0.07, by = 0.02), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.4,0.65,0.59,0.95),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Bet_5m, paramFrame$Bet_DCM, cex = 2, lwd = 2, 
     ylim = c(0, 0.0007), xlim = c(0, 0.0007), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(0.00005, 0.00065, expression(bold(beta^{"B"})), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.0001, 0.0007, by = 0.0002), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0001, 0.0007, by = 0.0002), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.72,0.97,0.59,0.95),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$PBM_5m, paramFrame$PBM_DCM, cex = 2, lwd = 2, 
     ylim = c(0.5, 3), xlim = c(0.5,3), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(0.7, 2.8, expression(bold("P"["M"]^"B")), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.5, 3, by = 0.5), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3, by = 0.5), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.08,0.33,0.16,0.52),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$IM_5m, paramFrame$IM_DCM, cex = 2, lwd = 2, 
     ylim = c(100, 550), xlim = c(100, 550), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(120, 500, expression(bold("I"["M"])), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(100, 500, by = 200), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(100, 500, by = 200), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.4,0.65,0.16,0.52),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$IK_5m, paramFrame$IK_DCM, cex = 2, lwd = 2, 
     ylim = c(10, 110), xlim = c(10, 110), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(20, 100, expression(bold("I"["K"])), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(20, 100, by = 20), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(20, 100, by = 20), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)
mtext("5 m", 1, line = 5.5, cex = 2.5, font = 2, las = 0)

par(xpd = TRUE)
legend("bottom", cex = 2.5, horiz = TRUE,
       legend = c(expression(bold("Mixed (Oct-May)")), 
                  expression(bold("Stratified (Jun-Sep)"))),
       fill = c("#1A85FF", "#D41159"), box.lwd = 0, inset = c(0, -0.4))
par(xpd = FALSE)

dev.off()




paramFrame$Month <- as.numeric(as.character(format(paramFrame$Date, "%m")))
cols <- rep("#D41159", nrow(paramFrame))
cols[which(paramFrame$Month < 6 | paramFrame$Month > 9)] <- "#1A85FF"

png("FigBin/fullPaneledComp_Four.png", width = 1400, height = 1510)
plot.new()
par(new = "TRUE",plt = c(0.09,0.5,0.58,0.96),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Alph_5m, paramFrame$Alph_DCM, col = cols, cex = 2, lwd = 2,
     ylim = c(0.01, 0.08), xlim = c(0.01, 0.08),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
arrows(paramFrame$Alph_5m, paramFrame$Alph_DCM+paramFrame$Alph_Sd_DCM, paramFrame$Alph_5m, 
       paramFrame$Alph_DCM-paramFrame$Alph_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
arrows(paramFrame$Alph_5m+paramFrame$Alph_Sd_5m, paramFrame$Alph_DCM, 
       paramFrame$Alph_5m-paramFrame$Alph_Sd_5m, paramFrame$Alph_DCM,
       length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
text(0.015, 0.075, expression(bold(alpha^{"B"})), cex = 2.5)
mtext("Chl max", 2, line = 6, cex = 2.5, font = 2, las = 0, 
      adj = -0.2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.01, 0.07, by = 0.02), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.01, 0.07, by = 0.02), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.57,0.98,0.58,0.96),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Bet_5m, paramFrame$Bet_DCM, cex = 2, lwd = 2, 
     ylim = c(0, 0.0007), xlim = c(0, 0.0007), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
arrows(paramFrame$Bet_5m, paramFrame$Bet_DCM+paramFrame$Bet_Sd_DCM, paramFrame$Bet_5m, 
       paramFrame$Bet_DCM-paramFrame$Bet_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
arrows(paramFrame$Bet_5m+paramFrame$Bet_Sd_5m, paramFrame$Bet_DCM, 
       paramFrame$Bet_5m-paramFrame$Bet_Sd_5m, paramFrame$Bet_DCM,
       length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
text(0.00005, 0.00065, expression(bold(beta^{"B"})), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.0001, 0.0007, by = 0.0002), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0001, 0.0007, by = 0.0002), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.09,0.5,0.13,0.51),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$PBM_5m, paramFrame$PBM_DCM, cex = 2, lwd = 2, 
     ylim = c(0.5, 3), xlim = c(0.5,3), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
arrows(paramFrame$PBM_5m, paramFrame$PBM_DCM+paramFrame$PBM_Sd_DCM, paramFrame$PBM_5m, 
       paramFrame$PBM_DCM-paramFrame$PBM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
arrows(paramFrame$PBM_5m+paramFrame$PBM_Sd_5m, paramFrame$PBM_DCM, 
       paramFrame$PBM_5m-paramFrame$PBM_Sd_5m, paramFrame$PBM_DCM,
       length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
text(0.7, 2.8, expression(bold("P"["M"]^"B")), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(0.5, 3, by = 0.5), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3, by = 0.5), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)


par(new = "TRUE",plt = c(0.57,0.98,0.13,0.51),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$IM_5m, paramFrame$IM_DCM, cex = 2, lwd = 2, 
     ylim = c(100, 550), xlim = c(100, 550), col = cols,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
arrows(paramFrame$IM_5m, paramFrame$IM_DCM+paramFrame$IM_Sd_DCM, paramFrame$IM_5m, 
       paramFrame$IM_DCM-paramFrame$IM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
arrows(paramFrame$IM_5m+paramFrame$IM_Sd_5m, paramFrame$IM_DCM, 
       paramFrame$IM_5m-paramFrame$IM_Sd_5m, paramFrame$IM_DCM,
       length=0.05, angle=90, code=3, 
       lwd = 2.5, col = cols)
text(120, 510, expression(bold("E"["M"])), cex = 2.5)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1,
     at = seq(100, 500, by = 200), font = 2, las = 1)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(100, 500, by = 200), font = 2, las = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
box(lwd=2)
mtext("5 m", 1, line = 6, cex = 2.5, font = 2, las = 0, adj = -0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

par(xpd = TRUE)
legend("bottom", cex = 2.5, horiz = TRUE,
       legend = c(expression(bold("Mixed (Oct-May)")), 
                  expression(bold("Stratified (Jun-Sep)"))),
       fill = c("#1A85FF", "#D41159"), box.lwd = 0)
par(xpd = FALSE)

dev.off()



#######CHL, PBM, PM Table
miniParamFrame <- paramFrame[, c("Date", "PBM_5m", "PBM_Sd_5m", "PBM_DCM", "PBM_Sd_DCM")]
PBMEtAlTab <- 
  merge(miniParamFrame, metaDatTab[, c("Date", "Chl_5m", "Chl_DCM")], 
        by = "Date", all = TRUE)
PBMEtAlTab$PM_5m <- PBMEtAlTab$PBM_5m*PBMEtAlTab$Chl_5m
PBMEtAlTab$PM_Sd_5m <- PBMEtAlTab$PBM_Sd_5m*PBMEtAlTab$Chl_5m
PBMEtAlTab$PM_DCM <- PBMEtAlTab$PBM_DCM*PBMEtAlTab$Chl_DCM
PBMEtAlTab$PM_Sd_DCM <- PBMEtAlTab$PBM_Sd_DCM*PBMEtAlTab$Chl_DCM
PBMEtAlTab[, 2:11] <- round(PBMEtAlTab[, 2:11], digits = 2)
PBMEtAlTab$PM_Comb_5m <- paste0(PBMEtAlTab$PM_5m, "+", PBMEtAlTab$PM_Sd_5m)
PBMEtAlTab$PM_Comb_DCM <- paste0(PBMEtAlTab$PM_DCM, "+", PBMEtAlTab$PM_Sd_DCM)
PBMEtAlTab$PBM_Comb_5m <- paste0(PBMEtAlTab$PBM_5m, "+", PBMEtAlTab$PBM_Sd_5m)
PBMEtAlTab$PBM_Comb_DCM <- paste0(PBMEtAlTab$PBM_DCM, "+", PBMEtAlTab$PBM_Sd_DCM)


#write.csv(PBMEtAlTab, "~/FlatheadPPs/PBMvsChlvsPM.csv", quote = FALSE)



####Multi-panel trends
png("FigBin/fullHierTronAlphBetaPBMBoth.png", width = 2800, height = 1000)
plot.new()
par(new = "TRUE",plt = c(0.08,0.33,0.56,0.94),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Alph_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.1),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Alph_5m+paramFrame$Alph_Sd_5m, paramFrame$Date, 
       paramFrame$Alph_5m-paramFrame$Alph_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.04, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0125, 0.1125, by = 0.025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 0.1, by = 0.025), labels = seq(0, 0.1, by = 0.025), font = 2)
mtext("5 m", 2, line = 9, font = 2, cex = 2.5, las = 0)
mtext(expression(bold(alpha^{"B"}~"(gC (gChl ")~bolditalic("a")~
                   bold(")"^{"-1"}~"h"^{"-1"}~"(μmol quanta m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}~")")), 3, 
      line = 1, cex = 2.5, font = 2, las = 1)

par(new = "TRUE",plt = c(0.08,0.33,0.13,0.51),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Alph_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 0.1),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Alph_DCM+paramFrame$Alph_Sd_DCM, paramFrame$Date, 
       paramFrame$Alph_DCM-paramFrame$Alph_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2, lwd.ticks = 3, tck = -0.04, padj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%y"), font = 2, las = 1)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.0125, 0.1125, by = 0.025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 0.1, by = 0.025), labels = seq(0, 0.1, by = 0.025), font = 2)
mtext("Chl max", 2, line = 9, font = 2, cex = 2.5, las = 0)



par(new = "TRUE",plt = c(0.4,0.65,0.56,0.94),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Bet_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 0.001),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Bet_5m+paramFrame$Bet_Sd_5m, paramFrame$Date, 
       paramFrame$Bet_5m-paramFrame$Bet_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.04, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.000125, 0.001125, by = 0.00025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 0.001, by = 0.00025), 
     labels = c("0", "2.5e-4", "5e-4", "7.5e-4", "1e-3"), 
     font = 2)
mtext(expression(bold(beta^{"B"}~"(gC (gChl ")~bolditalic("a")~
                   bold(")"^{"-1"}~"h"^{"-1"}~"(μmol quanta m"^{"-2"}~"s"^{"-1"}~")"^{"-1"}~")")), 3, 
      line = 1, cex = 2.5, font = 2, las = 1)

par(new = "TRUE",plt = c(0.4,0.65,0.13,0.51),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$Bet_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 0.001),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$Bet_DCM+paramFrame$Bet_Sd_DCM, paramFrame$Date, 
       paramFrame$Bet_DCM-paramFrame$Bet_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2, lwd.ticks = 3, tck = -0.04, padj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%y"), font = 2, las = 1)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.000125, 0.001125, by = 0.00025), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 0.001, by = 0.00025), 
     labels = c("0", "2.5e-4", "5e-4", "7.5e-4", "1e-3"), 
     font = 2)
mtext("Date", 1, line = 6, cex = 2.5, font = 2)

par(new = "TRUE",plt = c(0.72,0.97,0.56,0.94),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$PBM_5m, xaxt = "n", ylab = "", yaxt = "n",
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xlab = "", ylim = c(0, 3),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$PBM_5m+paramFrame$PBM_Sd_5m, paramFrame$Date, 
       paramFrame$PBM_5m-paramFrame$PBM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.04, padj = 1,
     at = dateLabels, labels = NA)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1), font = 2)
mtext(expression(bold("P"["M"]^{"B"}~"(gC (gChl ")~bolditalic("a")~
                   bold(")"^{"-1"}~"h"^{"-1"}~")")), 
      3, line = 1, cex = 2.5, font = 2)

par(new = "TRUE",plt = c(0.72,0.97,0.13,0.51),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(paramFrame$Date, paramFrame$PBM_DCM, xlab = "", ylab = "", yaxt = "n", 
     xlim = as.Date(c("2018-03-01", "2020-02-01")), xaxt = "n", ylim = c(0, 3),
     pch = 1, cex = 1, lwd = 2.5)
arrows(paramFrame$Date, paramFrame$PBM_DCM+paramFrame$PBM_Sd_DCM, paramFrame$Date, 
       paramFrame$PBM_DCM-paramFrame$PBM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2, lwd.ticks = 3, tck = -0.04, padj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%y"), font = 2, las = 1)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, hadj = 1, 
     at = seq(0.5, 3.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.1,
     at = seq(0, 3, by = 1), labels = seq(0, 3, by = 1), font = 2)

dev.off()










pbmDF <- read.csv("~/FlatheadPPs/PBMvsChlvsPM.csv", stringsAsFactors = FALSE)
pbmDF$Date <- as.Date(pbmDF$Date)
pbmDF$DOY <- as.numeric(format(pbmDF$Date, "%j"))
sortPBMdf <- pbmDF[order(pbmDF$DOY),]

pbmDF_5m <- sortPBMdf[!is.na(sortPBMdf$PBM_5m),]
pbmDF_DCM <- sortPBMdf[!is.na(sortPBMdf$PBM_DCM),]


png("FigBin/pbm_vs_pm.png", width = 1000, height = 750)
plot.new()

###5m###
par(new = "TRUE",plt = c(0.12,0.94,0.6,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(pbmDF_5m$DOY, pbmDF_5m$PBM_5m, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, col = "#DC3220",
     main = "", xlim = c(0, 365), ylim = c(0.5, 3))
arrows(pbmDF_5m$DOY, pbmDF_5m$PBM_5m+pbmDF_5m$PBM_Sd_5m, pbmDF_5m$DOY, 
       pbmDF_5m$PBM_5m-pbmDF_5m$PBM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
par(new = "TRUE")
plot(pbmDF_5m$DOY, pbmDF_5m$PM_5m, lty = 2, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, col = "#005AB5",
     main = "", xlim = c(0, 365), ylim = c(0.5, 3))
arrows(pbmDF_5m$DOY, pbmDF_5m$PM_5m+pbmDF_5m$PM_Sd_5m, pbmDF_5m$DOY, 
       pbmDF_5m$PM_5m-pbmDF_5m$PM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#005AB5")
box(lwd = 3)
axis(1, cex.axis=2.5, at = c(60, 152, 244, 335), tck = -0.04, 
     labels = NA,
     padj = 1, lwd.ticks = 3, font = 2)
axis(1, cex.axis=2.5, at = c(1, 32, 91, 121, 182, 213, 274, 305, 366), tck = -0.02, 
     labels = NA, padj = 1, lwd.ticks = 3)
text(-5, 2.9, "5 m", font = 2, cex = 2, pos = 4)

###DCM###
par(new = "TRUE",plt = c(0.12,0.94,0.18,0.56),las = 1, cex.axis = 2.5)
plot(pbmDF_DCM$DOY, pbmDF_DCM$PBM_DCM, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#DC3220", xlim = c(0, 365), ylim = c(0.5, 3.6))
arrows(pbmDF_DCM$DOY, pbmDF_DCM$PBM_DCM+pbmDF_DCM$PBM_Sd_DCM, pbmDF_DCM$DOY, 
       pbmDF_DCM$PBM_DCM-pbmDF_DCM$PBM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
par(new = "TRUE")
plot(pbmDF_DCM$DOY, pbmDF_DCM$PM_DCM, lty = 2, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#005AB5", xlim = c(0, 365), ylim = c(0.5, 3.6))
arrows(pbmDF_DCM$DOY, pbmDF_DCM$PM_DCM+pbmDF_DCM$PM_Sd_DCM, pbmDF_DCM$DOY, 
       pbmDF_DCM$PM_DCM-pbmDF_DCM$PM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#005AB5")
box(lwd = 3)
axis(1, cex.axis=2.5, at = c(60, 152, 244, 335), tck = -0.04, 
     labels = month.abb[seq(3,12, by = 3)],
     padj = 1, lwd.ticks = 3, font = 2)
axis(1, cex.axis=2.5, at = c(1, 32, 91, 121, 182, 213, 274, 305, 366), tck = -0.02, 
     labels = NA, padj = 1, lwd.ticks = 3)
mtext("Photosynthetic Rate", 
      2, line = 6, cex = 2.5, font = 2, las = 0, adj = -33)
text(-5, 3.4, "Chl max", font = 2, cex = 2, pos = 4)
par(xpd = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(-0.7, -0.85, 
       legend = c(expression(bold("P"["M"]^{"B"}~"(g C (g Chl ")~bolditalic("a")~
                               bold(")"^{"-1"}~"h"^{"-1"}~")")), 
                  expression(bold("P"["M"]~"(g C h"^{"-1"}~")"))), 
       fill = c("#DC3220", "#005AB5"), horiz = TRUE, cex = 2.5, box.lwd = 0)

dev.off()









pbmDF <- read.csv("~/FlatheadPPs/PBMvsChlvsPM.csv", stringsAsFactors = FALSE)
pbmDF$Date <- as.Date(pbmDF$Date)

pbmDF_5m <- pbmDF[!is.na(pbmDF$PBM_5m),]
pbmDF_DCM <- pbmDF[!is.na(pbmDF$PBM_DCM),]


png("FigBin/pbm_vs_pm_date.png", width = 1000, height = 750)
plot.new()

###5m###
par(new = "TRUE",plt = c(0.12,0.94,0.6,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(pbmDF_5m$Date, pbmDF_5m$PBM_5m, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxt = "n",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, col = "#DC3220",
     main = "", ylim = c(0.5, 3), xlim = c(as.Date("2018-03-01"), as.Date("2020-02-01")))
arrows(pbmDF_5m$Date, pbmDF_5m$PBM_5m+pbmDF_5m$PBM_Sd_5m, pbmDF_5m$Date, 
       pbmDF_5m$PBM_5m-pbmDF_5m$PBM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
par(new = "TRUE")
plot(pbmDF_5m$Date, pbmDF_5m$PM_5m, lty = 2, pch = 18,
     type = "b", yaxt = "n", xlab = "", ylab = "", xaxt = "n",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, col = "#005AB5",
     main = "", ylim = c(0.5, 3), xlim = c(as.Date("2018-03-01"), as.Date("2020-02-01")))
arrows(pbmDF_5m$Date, pbmDF_5m$PM_5m+pbmDF_5m$PM_Sd_5m, pbmDF_5m$Date, 
       pbmDF_5m$PM_5m-pbmDF_5m$PM_Sd_5m, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#005AB5")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.04, hadj = 1.1,
     at = dateLabels, labels = NA, font = 2, las = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
text(as.Date("2018-03-01"), 2.9, "5 m", font = 2, cex = 2, pos = 4)

###DCM###
par(new = "TRUE",plt = c(0.12,0.94,0.18,0.56),las = 1, cex.axis = 2.5)
plot(pbmDF_DCM$Date, pbmDF_DCM$PBM_DCM, pch = 18,
     type = "b", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, xaxt = "n", 
     main = "", col = "#DC3220", ylim = c(0.5, 3.6), 
     xlim = c(as.Date("2018-03-01"), as.Date("2020-02-01")))
arrows(pbmDF_DCM$Date, pbmDF_DCM$PBM_DCM+pbmDF_DCM$PBM_Sd_DCM, pbmDF_DCM$Date, 
       pbmDF_DCM$PBM_DCM-pbmDF_DCM$PBM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#DC3220")
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
par(new = "TRUE")
plot(pbmDF_DCM$Date, pbmDF_DCM$PM_DCM, lty = 2, pch = 18,
     type = "b", yaxt = "n", xlab = "", ylab = "", xaxt = "n",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#005AB5", ylim = c(0.5, 3.6), 
     xlim = c(as.Date("2018-03-01"), as.Date("2020-02-01")))
arrows(pbmDF_DCM$Date, pbmDF_DCM$PM_DCM+pbmDF_DCM$PM_Sd_DCM, pbmDF_DCM$Date, 
       pbmDF_DCM$PM_DCM-pbmDF_DCM$PM_Sd_DCM, length=0.05, angle=90, code=3, 
       lwd = 2.5, col = "#005AB5")
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2020-01-01"), lty = 2, col = "gray75", lwd = 3)
box(lwd = 3)
axis(1, cex.axis=2, lwd.ticks = 3, tck = -0.04, padj = 1.1,
     at = dateLabels, labels = format(dateLabels, "%m-%y"), font = 2, las = 1)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = dateLabelSequence, labels = NA)
mtext("Photosynthetic Rate", 
      2, line = 6, cex = 2.5, font = 2, las = 0, adj = -33)
text(as.Date("2018-03-01"), 3.4, "Chl max", font = 2, cex = 2, pos = 4)
par(xpd = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(-0.7, -0.85, 
       legend = c(expression(bold("P"["M"]^{"B"}~"(g C (g Chl ")~bolditalic("a")~
                               bold(")"^{"-1"}~"h"^{"-1"}~")")), 
                  expression(bold("P"["M"]~"(g C h"^{"-1"}~")"))), 
       fill = c("#DC3220", "#005AB5"), horiz = TRUE, cex = 2.5, box.lwd = 0)

dev.off()

