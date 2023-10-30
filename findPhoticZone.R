coll14C <- read.csv("~/FlatheadPublic/NPP/calculatedNPP_Jul2023.csv", 
                    stringsAsFactors = FALSE)
coll14C$X <- NULL
mean14C <-aggregate(NPP~Date+Depth, coll14C, FUN = mean)
mean14C$Date <- as.Date(mean14C$Date)
mean14C$Depth <- as.numeric(as.character(mean14C$Depth))

allDates <- unique(mean14C$Date)
photicDeps <- data.frame("Date" = allDates, "PhoticZone_01" = NA, "PhoticZone_05" = NA,
                         "PhoticZone_03" = NA, "PhoticZone_0" = NA, "N" = NA)

for (i in 1:length(allDates)){
  subDf <- mean14C[mean14C$Date == allDates[i],]
  subDf <- subDf[!is.na(subDf$NPP),]
  
  photicDeps[photicDeps$Date == allDates[i], "PhoticZone_01"] <- 
    round(approx(subDf$NPP, subDf$Depth, xout=0.1)$y)
  
  photicDeps[photicDeps$Date == allDates[i], "PhoticZone_03"] <- 
    round(approx(subDf$NPP, subDf$Depth, xout=0.3)$y)
  
  photicDeps[photicDeps$Date == allDates[i], "PhoticZone_05"] <- 
    round(approx(subDf$NPP, subDf$Depth, xout=0.5)$y)
  
  photicDeps[photicDeps$Date == allDates[i], "PhoticZone_0"] <- 
    round(approx(subDf$NPP, subDf$Depth, xout=0)$y)
  
  photicDeps[photicDeps$Date == allDates[i], "N"] <- nrow(subDf)
}

photicDepZero <- photicDeps[!is.na(photicDeps$PhoticZone_0),]


##K
kVals <- read.csv("~/FlatheadPublic/lightAttenuation_Jan2023_lm.csv", 
                  stringsAsFactors = FALSE)
kVals$Date <- as.Date(kVals$Date)
kVals$N <- NULL

##Incident light
##umol photon m^-2 s^-1
YBHill <- read.csv("~/FlatheadPublic/Weather/YBPoint_PAR_20110101to20191231.csv",
                   stringsAsFactors = FALSE)
YBHill$DateTime <- as.POSIXlt(YBHill$timestamp, format = "%m/%d/%Y %H:%M")
YBHill$id <- NULL
YBHill$parameterID <- NULL
YBHill$parameterName <- NULL
YBHill$parameterUnits <- NULL #µmol photons/s/m^2
YBHill$timestamp <- NULL
YBHill$Date <- as.Date.POSIXlt(YBHill$DateTime)

dailySumYBPoint <- aggregate(parameterValue ~ Date, data = YBHill, FUN = sum)
dailySumYBPoint$DailySum <- dailySumYBPoint$parameterValue * 60 * 5 /  1000
  
lightDf <- merge(kVals, YBHill, by = "Date")
lightDf <- lightDf[order(lightDf$DateTime),]

dailySum <- merge(kVals, dailySumYBPoint, by = "Date")
dailySum <- dailySum[order(dailySum$Date),]

##umol photons m2 s-1
photicDepZero$maxLight <- NA
for (i in 1:nrow(photicDepZero)){
  lightSub <- lightDf[lightDf$Date == photicDepZero$Date[i],]
  lightSub <- lightSub[!is.na(lightSub$parameterValue),]
  kval <- lightSub$kVal[6]
  maxLight <- max(lightSub$parameterValue)
  
  photicDepZero[i, "maxLight"] <- maxLight*exp(-kval*photicDepZero[i, "PhoticZone_0"])
}

#umol photons m-2 s-1
isolumeTarg <- mean(photicDepZero$maxLight, na.rm = TRUE)


##Daily sum
photicDepZero$dailyLight <- NA
for (i in 1:nrow(photicDepZero)){
  dailySumMini <- dailySum[dailySum$Date == photicDepZero$Date[i],]
  dailySumMini <- dailySumMini[!is.na(dailySumMini$DailySum),]
  kval <- dailySumMini$kVal
  maxLight <- dailySumMini$DailySum
  
  if (length(kval) > 0){
    photicDepZero[i, "dailyLight"] <- maxLight*exp(-kval*photicDepZero[i, "PhoticZone_0"])
  }
}

#mmol photons m-2 day-1
isolumeTarg <- mean(photicDepZero$dailyLight, na.rm = TRUE)



####max light of 3 umol photons m-2 s-1


allModeledFull <- read.csv("allValuesModeledFull_m3_withChl_5mOnly.csv", 
                           stringsAsFactors = FALSE)

allModeledFull$Date <- as.Date(allModeledFull$Date)

targSampDates <- data.frame("Date" = allModeledFull$Date, "PhoticDepth" = NA)
for (i in 1:nrow(targSampDates)){
  lightDat <- lightDf[lightDf$Date == targSampDates[i, "Date"],]
  maxLight <- max(lightDat$parameterValue)
  kVal <- mean(lightDat$kVal)
  
  targSampDates[i, "PhoticDepth"] <- abs(log(3.3/maxLight)/kVal)
}

targSampDates$Month <- as.numeric(as.character(format(targSampDates$Date, "%m")))

targSampDates$DailySum <- NA
for (i in 1:nrow(targSampDates)){
  dailySumLight <- dailySum[dailySum$Date == targSampDates[i, "Date"], "DailySum"]
  kVal <- mean(dailySum$kVal)
  
  targSampDates[i, "DailySum"] <- abs(log(99/dailySumLight)/kVal)
}

#write.csv(targSampDates, "~/FlatheadPhotos/PhoticZoneDepth.csv", quote = FALSE)

i <- 12
mean(targSampDates[targSampDates$Month == i, "PhoticDepth"])
sd(targSampDates[targSampDates$Month == i, "PhoticDepth"])


mixingRef <- read.csv("~/FlatheadPhotos/MLDvsDCMDepthComparison.csv",
                      stringsAsFactors = FALSE)
mixingRef$Date <- as.Date(mixingRef$Date)

targSampDatesTrim <- targSampDates[, c("Date", "PhoticDepth")]

allDepsData <- merge(mixingRef, targSampDatesTrim, by = "Date")

plot(allDepsData$Date, allDepsData$DCM, ylim = c(40, 0))
points(allDepsData$Date, allDepsData$PhoticDepth, col = "red")
points(allDepsData$Date, allDepsData$td05, col = "blue")

#write.csv(allDepsData, "~/FlatheadPhotos/threeDepths.csv", quote = FALSE)
