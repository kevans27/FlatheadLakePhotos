library(readxl)

#Read in NPP
fullNPP <- read.csv("~/TMNT/14C/Historic Data/NPP.csv")
wantNPP <- fullNPP[,c("CollectDate", "Param", "CorrectedReportedResult")]

saveWantNPP <- wantNPP
#wantNPP <- saveWantNPP

wantNPP$Date <- as.Date(wantNPP$CollectDate, format = "%m/%d/%y")
wantNPP$CollectDate <- NULL
names(wantNPP)[names(wantNPP) == 'CorrectedReportedResult'] <- 'DPM'
wantNPP$Status <- NA
wantNPP[grepl("Dark", wantNPP$Param), "Status"] <- "Dark"
wantNPP[grepl("Light", wantNPP$Param), "Status"] <- "Light"
wantNPP[grepl("Stock", wantNPP$Param), "Status"] <- "Stock"
wantNPP[grepl("Prekill", wantNPP$Param), "Status"] <- "Prekill"
subNPP <- wantNPP[!is.na(wantNPP$Status),]
subNPP$Depth <- NA
depths <- c("1m", "5m", "2_5m", "2.5m", "10m", "20m", "30m", "Stock", "Prekill")
subin <- c("1", "5", "2.5", "2.5", "10", "20", "30", NA, NA)
for (i in 1:length(depths)){
  dep <- depths[i]
  sub <- subin[i]
  subNPP[grepl(dep, subNPP$Param), "Depth"] <- sub
}

#Find the average Stock by date and add on to the wantNPP df
subStatus <- subNPP[subNPP$Status == "Stock",]
subStatus$DPM <- as.numeric(as.character(subStatus$DPM))
subStatus <- subStatus[!is.na(subStatus$DPM),]

prekillStatus <- subNPP[subNPP$Status == "Prekill",]
prekillStatus$DPM <- as.numeric(as.character(prekillStatus$DPM))
prekillStatus <- prekillStatus[!is.na(prekillStatus$DPM),]
prekillStatus <- prekillStatus[, c("DPM", "Date")]

monthlymeans <- aggregate(DPM ~ Date, subStatus, mean)
subNPP <- merge(subNPP, monthlymeans, by = "Date")
wantNPPAll <- merge(subNPP, prekillStatus, by = "Date")
names(wantNPPAll) <- c("Date", "Param", "DPM", "Status", "Depth", "Stock", "Prekill")
wantNPPAll <- wantNPPAll[!is.na(wantNPPAll$Depth),]

lights <- wantNPPAll[wantNPPAll$Status == "Light",]
lights$Status <- NULL
lights$DPM <- as.numeric(lights$DPM)

darks <- wantNPPAll[wantNPPAll$Status == "Dark",]
darks$Status <- NULL
darks$DPM <- as.numeric(darks$DPM)


parRatios <- read.csv("~/FlatheadPublic/NPPIncubationPARRatios_2011-2022.csv",
                      stringsAsFactors = FALSE)
goodPAR <- parRatios[, c("Date", "ParRat")]
goodPAR$Date <- as.Date(goodPAR$Date)
lights <- merge(lights, goodPAR, by = "Date")


vivianDIC <- read_excel("~/FlatheadPublic/ChurchLab/FMP_DIC.xlsx")
vivianDIC$Date <- as.Date(vivianDIC$`Date - alkalinity measurement`)
vivDIC <- vivianDIC[vivianDIC$Date > as.Date("2011-01-01"),]
miniDIC <- vivDIC[, c("Date", "Depth (m)", "Alkalinity conc. (mg CaCO3 / L)", "Temp (˚C)",
                      "pH", "Derived DIC (mmol C / L)")]
colnames(miniDIC) <- c("Date", "Depth", "Alk_mg_CaCO3_L", "Temp", "pH", "vivDIC")
miniDIC <- miniDIC[miniDIC$Depth == 5,]

avgdDIC <- aggregate(. ~ Date, data = miniDIC, FUN = mean)
avgdDIC$Salinity <- 0
avgdDIC$Phosphate <- 0
avgdDIC$Atm <- 0.9
avgdDIC$HydroP <- 0.491 #Hydrostatic pressure in bar at 5 m

library(readxl)
allSilica <- read_excel("~/FlatheadPublic/FMPSilica_2019.xlsx")
miniSilica <- allSilica[1121:nrow(allSilica), ] 
miniSilica <- miniSilica[as.numeric(as.character(miniSilica$`Start Depth`)) <= 30,]
miniSilica <- miniSilica[, c("Date", "Value")] 
miniSilica$Value <- as.numeric(as.character(miniSilica$Value))
miniSilica <- aggregate(. ~ Date, miniSilica, FUN = mean, na.rm = TRUE)
#In mg/L, divide by 28, / 1000
miniSilica$Value <- as.numeric(as.character(miniSilica$Value))/28/1000
miniSilica$Date <- as.Date(miniSilica$Date, format = "%m/%d/%Y")
colnames(miniSilica) <- c("Date", "Silica")

avgdDIC <- merge(avgdDIC, miniSilica, by = "Date", all.x = TRUE)


#Convert alk from mg/L caco2 to mol/L = /0.02/1000
avgdDIC$Alk_mol_L <- avgdDIC$Alk_mg_CaCO3_L/100/1000*2

#pH is probably NBS, in which case pH(total) = pH(NBS) - 0.15
avgdDIC$pH_Total <- avgdDIC$pH - 0.15

library("seacarb")

##Carb wants: pH, ALK (flag 8), temp, salinity, atmospheric pressure, phosphate concentration (0), silicate concentration

savedDIC <- carb(8, avgdDIC$pH_Total, avgdDIC$Alk_mol_L, S = avgdDIC$Salinity, 
                 T = avgdDIC$Temp, Patm = avgdDIC$Atm, Pt = avgdDIC$Phosphate, 
                 Sit = avgdDIC$Silica, P = avgdDIC$HydroP)
#DIC is in mol/kg, convert to mmol/L
savedDIC$DIC <- savedDIC$DIC*1000
savedDIC$Date <- avgdDIC$Date


savedDIC$VivDIC <- avgdDIC$vivDIC
miniDIC <- savedDIC[, c("Date", "DIC")]

darksModern <- darks[darks$Date >= as.Date("2010-01-01"),]
darksTrim <- aggregate(DPM ~ Date + Depth, data = darksModern, FUN = mean)

lightsNew <- merge(lights, darksTrim, by = c("Date", "Depth"), all.x = TRUE)
colnames(lightsNew)[ncol(lightsNew)] <- "DarkDPM"
colnames(lightsNew)[4] <- "LightDPM"

letsroll <- merge(lightsNew, miniDIC, by = "Date")

NPPFunc <- function(dpm, stock, prekill, dic, parrat){
  K1 <- 1.06
  toLiters <- 1000
  NPP <- (dpm-prekill)/stock*dic*K1*parrat*toLiters
}

letsroll$NPP <- NPPFunc(as.numeric(as.character(letsroll$LightDPM)), 
                        as.numeric(as.character(letsroll$Stock)), 
                        as.numeric(as.character(letsroll$DarkDPM)),
                        letsroll$DIC, letsroll$ParRat)

letsroll <- letsroll[!is.infinite(letsroll$NPP),]
letsroll <- letsroll[!is.na(letsroll$NPP),]

#write.csv(letsroll, "~/FlatheadPublic/NPP/calculatedNPP_Jul2023.csv", quote = FALSE)


m30s <- letsroll[letsroll$Depth == 30,]
which(m30s$LightDPM <= m30s$DarkDPM)

plot(m30s$Date, m30s$LightDPM)
points(m30s$Date, m30s$DarkDPM, col = "red")
