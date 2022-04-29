#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(tidyverse)
library(reshape2)
library(sf)

#-----------#
#-Load Data-#
#-----------#

Data <- NULL
files <- list.files(path = "~/HMSNO/DataFormat/RawData", pattern = "Data", full.names = TRUE)
for(i in 1:length(files)){
  sheets <- excel_sheets(files[i])
  for(j in 1:length(sheets)){
    Temp <- read_excel(files[i], sheet = sheets[j])
    Temp <- Temp %>% mutate(species = as.character(sheets[j]), park = as.factor(i)) %>% 
      select(park, `deployment ID`, Year, species, ls(Temp, pattern = "R"), density)
    Data <- bind_rows(Data, Temp)
  }
}

#UDZ ungulates 2018-2019
sheets <- excel_sheets("~/HMSNO/DataFormat/RawData/UDZ_ungulates_2018_2019.xlsx")
UDZdf <- NULL
for(j in 3:6){
  Temp <- read_excel("~/HMSNO/DataFormat/RawData/UDZ_ungulates_2018_2019.xlsx",
                     sheet = sheets[j])
  UDZdf <- bind_rows(UDZdf, Temp)
}

Data <- bind_rows(Data, UDZdf %>% select(Trap...2, Year, Species, ls(Temp, pattern = "R")) %>% select(-R7) %>%
  rename(`deployment ID` = Trap...2, species = Species) %>%
  mutate(species = recode(species, "Cephalophus harveyi" = "harveyi", "Cephalophus spadix" = "spadix", 
                          "Nesotragus moschatus" = "moschatus", "Tragelaphus scriptus" = "scriptus"),
         park = as.factor(5), density = 0) %>%
  relocate(park, `deployment ID`))

#Remaining ungulates
sheets <- excel_sheets("~/HMSNO/DataFormat/RawData/Remaining_ungulates.xlsx")
RUdf <- NULL
for(k in 1:length(sheets)){
  Temp <- read_excel("~/HMSNO/DataFormat/RawData/Remaining_ungulates.xlsx",
                     sheet = sheets[k])
  RUdf <- bind_rows(RUdf, Temp)
}

Data <- bind_rows(Data, RUdf %>% select(park, `deployment ID`, Year, Species, ls(Temp, pattern = "R"), density) %>%
  rename(species = Species) %>%
  mutate(species = recode(species, "Tragelaphus scriptus" = "scriptus",
                                   "Tragelaphus spekii" = "spekii"),
         park = as.factor(park)))

#Add missing UDZ data
sheets <- excel_sheets("~/HMSNO/DataFormat/RawData/Bushbuck_suni_UDZ.xlsx")
UDZdf <- NULL
for(k in 3:4){
  Temp <- read_excel("~/HMSNO/DataFormat/RawData/Bushbuck_suni_UDZ.xlsx",
                     sheet = sheets[k])
  Temp$park <- 5
  UDZdf <- bind_rows(UDZdf, Temp)
}

Data <- bind_rows(Data, UDZdf %>% select(park, `deployment ID`, Year, species, ls(Temp, pattern = "R"), density) %>%
        mutate(park = as.factor(park)))


  
#Site data
SiteData <- NULL
file <- list.files(path = "~/HMSNO/DataFormat/RawData", pattern = "metadata", full.names = TRUE)
sheets <- excel_sheets(file)
for(j in 1:length(sheets)){
  Temp <- read_excel(file, sheet = sheets[j])
  Temp <- Temp %>% select(`deployment ID`, Year, ls(Temp, pattern = "R"), 
                          elevation, edge, days)
  SiteData <- bind_rows(SiteData, Temp)
}

SiteData <- bind_rows(SiteData, read_excel("~/HMSNO/DataFormat/RawData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 1) %>%
                       select(`deployment ID`, days, elevation, edge, Year)) %>%
  bind_rows(., read_excel("~/HMSNO/DataFormat/RawData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 2) %>%
              select(`deployment ID`, days, elevation, edge, Year)) %>%
  distinct() %>%
  select(`deployment ID`, days, elevation, edge, Year)

SiteData <- bind_rows(SiteData, SiteData %>% filter(grepl("CT-UDZ", `deployment ID`) & Year >= 2016) %>%
        mutate(Year = ifelse(Year == 2016, 2018, 2019),
               days = read_excel("~/HMSNO/DataFormat/RawData/UDZ_ungulates_2018_2019.xlsx", sheet = 4) %>%
                      select(trapdays) %>% .$trapdays))

#Add GPS data
sites <- NULL
sheets <- excel_sheets("~/HMSNO/DataFormat/RawData/trap locations and dates.xlsx")
for(i in 1:length(sheets)){
  Temp <- read_excel(path = "~/HMSNO/DataFormat/RawData/trap locations and dates.xlsx", sheet = sheets[i])
  Temp <- Temp %>% mutate(latitude = as.numeric(as.character(latitude)),
                          longitude = as.numeric(as.character(longitude)),
                          start_date = as.Date(start_date),
                          end_date = as.Date(end_date)) 
  sites <- bind_rows(sites, Temp)
}

sites$deployment <- str_remove(sites$deployment, "CT-NNP-1-")
sites$deployment <- str_remove(sites$deployment, "NNP-2015-")
tmp <- sites$deployment[grep("TR", sites$deployment)]
tmp[-grep("-", tmp)] <- gsub("(....)(.*)","\\1-\\2",tmp[-grep("-", tmp)])
sites$deployment[grep("TR", sites$deployment)] <- tmp
tmp <- sites$deployment[grep("TR-", sites$deployment)]
tmp <- stringr::str_replace(tmp, "-", "")
sites$deployment[grep("TR-", sites$deployment)] <- tmp
sites <- sites %>% rename(`deployment ID` = deployment, Year = year)
SiteData <- left_join(SiteData, sites, by = c("deployment ID", "Year"))

#Set coordinate system
SiteData <- SiteData %>% group_by(`deployment ID`) %>%
  fill(longitude, latitude, start_date, end_date, .direction = "downup")

SiteData <- SiteData %>% 
  drop_na(longitude|latitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326)

#-------------#
#-Format Data-#
#-------------#

#Turn dashes to NAs
Data[,5:10] <- apply(Data[,5:10], 2, function(y) as.numeric(gsub("-", NA, y)))

#Convert CT-UDZ-1-14.1 to CT-UDZ-1-14
Data$`deployment ID`[Data$`deployment ID` == "CT-UDZ-1-14.1"] <- "CT-UDZ-1-14"

#Remove CT-UDZ-3-21
Data <- Data %>% filter(`deployment ID` != "CT-UDZ-3-21")

#Remove 2017 from VIR
Data <- Data %>% filter(!(park == 6 & Year == 2017))

#Reorder park levels
park.levels <- c("1" = "5", "2" = "6", "3" = "4", "4" = "1", "5" = "3", "6" = "2")

Data <- Data %>% mutate(`deployment ID` = ifelse(species == "sylvilocutor" & Year == "2010" & park == 1,
                                                 gsub("TEAM-001", "CT", `deployment ID`), `deployment ID`),
                        `deployment ID` = ifelse(species == "sylvilocutor" & Year == "2010" & park == 1,
                                                 gsub("-2010", "", `deployment ID`), `deployment ID`),
                        park = fct_recode(factor(park, levels = c("5","6","4","1","3","2")), !!!park.levels)) %>% 
                arrange(park, Year, `deployment ID`)


Data <- left_join(Data, SiteData, by = c("deployment ID", "Year"))

#Control for days sampled
Data$days[Data$days < 0] <- 0 #Set negative days to 0
Data$days[Data$days > 30] <- 30 #Truncate days at max 30 days
Data$days[Data$days == 0] <- NA #Set 0 days to NA (ie, no sampling)
Data[is.na(Data$days),5:10] <- NA #Set all occ data to NA for sites with NA days (ie, no sampling)
Data <- Data %>% mutate(R6 = replace(R6, which(days>=25 & is.na(R6)), 0), #Replace values with NAs or 0s
                        R6 = replace(R6, which(days<25), NA),
                        R5 = replace(R5, which(days>=20 & is.na(R5)), 0),
                        R5 = replace(R5, which(days<20), NA),
                        R4 = replace(R4, which(days>=15 & is.na(R4)), 0),
                        R4 = replace(R4, which(days<15), NA),
                        R3 = replace(R3, which(days>=10 & is.na(R3)), 0),
                        R3 = replace(R3, which(days<10), NA),
                        R2 = replace(R2, which(days>=5 & is.na(R2)), 0),
                        R2 = replace(R2, which(days<5), NA),
                        R1 = replace(R1, which(days>0 & is.na(R1)), 0))

Data$days[is.na(Data$days)] <- 0 #Reset no sampling to zero days

Data <- Data %>% drop_na(days) %>%
  mutate(parkID = as.numeric(park),
         yearID = as.numeric(as.factor(Year)),
         specID = as.numeric(as.factor(species))) %>%
  group_by(parkID) %>%
  mutate(siteID = as.numeric(factor(`deployment ID`, levels = unique(`deployment ID`)))) %>%
  ungroup(parkID) %>%   
  arrange(parkID, specID, siteID, yearID)

#---------#
#-Indices-#
#---------#

#Number of species
nspecs <- max(Data$specID)

#Number of parks
nparks <- max(Data$parkID)

#First park for each park
parkS <- as.numeric(Data %>% group_by(specID) %>%
                      summarize(parkS = min(parkID) - min(Data$parkID) + 1) %>%
                      select(parkS) %>% .$parkS)

#Last park for each park
parkE <- as.numeric(Data %>% group_by(specID) %>%
                      summarize(parkE = max(parkID) - min(Data$parkID) + 1) %>%
                      select(parkE) %>% .$parkE)
#Max number of sites
nsites <- max(Data$siteID)

#Number of sites/park
nsite <- as.numeric(Data %>% group_by(parkID) %>% 
                      summarize(nsite = n_distinct(siteID)) %>%
                      select(nsite) %>% .$nsite)

#Number of years
nyrs <- max(Data$yearID)

#First year for each park
nstart <- as.numeric(Data %>% group_by(parkID) %>%
                       summarize(nstart = min(yearID) - min(Data$yearID) + 1) %>%
                       select(nstart) %>% .$nstart)

#Last year for each park
nend <- as.numeric(Data %>% group_by(parkID) %>%
                     summarize(nend = max(yearID) - min(Data$yearID) + 1) %>%
                     select(nend) %>% .$nend)

#Number of replicates
nreps <- 6

#Nested indices
yr <- Data$yearID
site <- Data$siteID
park <- Data$parkID
spec <- Data$specID

#--------------#
#-Extract data-#
#--------------#

#Occupancy
y <- array(NA, dim = c(nreps, nyrs, nsites, nparks, nspecs))

for(i in 1:dim(Data)[1]){
  y[1:6,yr[i],site[i],park[i],spec[i]] <- as.numeric(Data[i,5:10])
}

#Covariates
Cov <- Data %>% group_by(parkID, siteID, yearID) %>%
  select(days, edge, elevation, geometry, start_date, end_date) %>% 
  distinct(parkID, siteID, yearID, .keep_all = TRUE)


#Fix start_date years that do not match year ID
for(i in 1:dim(Cov)[1]){
  if(as.numeric(as.factor(format(Cov$start_date[i], "%Y"))) != Cov$yearID[i]){
    lubridate::year(Cov$start_date[i]) <- Cov$yearID[i] + 2008
  }
}

Cov <- Cov %>% mutate(end_date_old = end_date,
                      end_date = start_date + 30)


# edge <- elevation <- array(NA, dim = c(nyrs, nsites, nparks))
# days <- array(NA, dim = c(nreps, nyrs, nsites, nparks))

for(i in 1:dim(Cov)[1]){
  
  if(Cov$days[i] > 25){
    tmp <- c(rep(5,5), (Cov$days[i] - 25))
  }else{
    if(Cov$days[i] <= 25 & Cov$days[i] > 20){
      tmp <- c(rep(5,4), (Cov$days[i] - 20), 0)
    }else{
      if(Cov$days[i] <= 20 & Cov$days[i] > 15){
        tmp <- c(rep(5,3), (Cov$days[i] - 15), rep(0, 2))
      }else{
        if(Cov$days[i] <= 15 & Cov$days[i] > 10){
          tmp <- c(rep(5,2), (Cov$days[i] - 10), rep(0, 3))
        }else{
          if(Cov$days[i] <= 10 & Cov$days[i] > 5){
            tmp <- c(1, (Cov$days[i] - 5), rep(0, 4))
          }else{
            if(Cov$days[i] <= 5 & Cov$days[i] > 0){
              tmp <- c(Cov$days[i], rep(0, 5))
            }else{
              if(Cov$days[i] == 0 | is.na(Cov$days[i])){
                tmp <- rep(0, 6)
              }
            }
          }
        }
      }
    }
  }
  
  days[ , Cov$yearID[i], Cov$siteID[i], Cov$parkID[i]] <- tmp
  # edge[Cov$yearID[i], Cov$siteID[i], Cov$parkID[i]] <- Cov$edge[i]
  # elevation[Cov$yearID[i], Cov$siteID[i], Cov$parkID[i]] <- Cov$elevation[i]
}


days.scaled <- (days - mean(days, na.rm = TRUE))/sd(days, na.rm = TRUE)
# edge.scaled <- (edge - mean(edge, na.rm = TRUE))/sd(edge, na.rm = TRUE)
# elevation.scaled <- (elevation - mean(elevation, na.rm = TRUE))/sd(elevation, na.rm = TRUE)
# density <- as.numeric(Data %>% group_by(park) %>% distinct(density) %>% select(density) %>% .$density)

#----------#
#-Rainfall-#
#----------#
# Geomdf <- st_as_sf(Cov)
# 
# Geomdf <- Geomdf %>% mutate(
#   annual_precip = 0,
# )
# 
# year <- as.character(rep(2009:2020, each = 12))
# month <- as.character(rep(sprintf('%0.2d', 1:12), length(2009:2020)))
# dates <- paste0(year, "-", month)
# year_ID <- as.numeric(as.factor(year))
# 
# 
# #Code to download CHIRPS data
# for(i in 1:length(year)){
#   url <- paste0("https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs/chirps-v2.0.", year[i], ".", month[i], ".tif.gz")
#   httr::GET(url = url,
#             httr::write_disk(path = paste0(getwd(), "/DataFormat/CHIRPS/", year[i], month[i], ".tif.gz"), overwrite = TRUE))
# }
# 
# filenames <- list.files("~/HMSNO/DataFormat/CHIRPS/", pattern = "chirps", full.names = TRUE)
# for(i in 1:length(filenames)){
#   rainfall <- raster::raster(filenames[i])
#   raster::values(rainfall)[raster::values(rainfall) < 0] <- NA
#   Geomdf$cell <- tabularaster::cellnumbers(rainfall, Geomdf)$cell_
#   Geomdf <- Geomdf %>% mutate(annual_precip = ifelse(yearID == year_ID[i], raster::extract(rainfall, cell) + annual_precip, annual_precip))
# }
#   
# precip <- array(NA, dim = c(22, nsites, nparks))
# 
# for(i in 1:dim(Geomdf)[1]){
#   precip[Geomdf$yearID[i], Geomdf$siteID[i], Geomdf$parkID[i]] <- Geomdf$annual_precip[i]
# }

#--------------#
#-Compile data-#
#--------------#

HMSNO.data <- list(y = y, 
                   days.scaled = days.scaled)

HMSNO.con <- list(nspecs = nspecs, parkS = parkS, parkE = parkE,
                  nsite = nsite, nreps = nreps, nstart = nstart, nend = nend)

save(HMSNO.data, file = "~/HMSNO/DataFormat/HMSNO.data.Rdata")

save(HMSNO.con, file = "~/HMSNO/DataFormat/HMSNO.con.Rdata")