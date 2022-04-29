#TODO: Remove all prediction data

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

print(colnames(UDZdf))

Data <- bind_rows(Data, UDZdf %>% select(Trap, Year, Species, ls(Temp, pattern = "R")) %>% select(-R7) %>%
  rename(`deployment ID` = Trap, species = Species) %>%
  mutate(species = recode(species, "Cephalophus harveyi" = "harveyi", "Cephalophus spadix" = "spadix", 
                          "Nesotragus moschatus" = "moschatus", "Tragelaphus scriptus" = "scriptus"),
         park = as.factor(5), density = 0)) # %>%
  #relocate(park, `deployment ID`))

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

#Add gps data

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

#Reorder park levels
# park.levels <- c("1" = "2", "2" = "3", "3" = "4", "4" = "1", "5" = "6", "6" = "5")
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


edge <- elevation <- array(NA, dim = c(nyrs, nsites, nparks))
days <- array(NA, dim = c(nreps, nyrs, nsites, nparks))

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
  edge[Cov$yearID[i], Cov$siteID[i], Cov$parkID[i]] <- Cov$edge[i]
  elevation[Cov$yearID[i], Cov$siteID[i], Cov$parkID[i]] <- Cov$elevation[i]
}

days <- (days - mean(days, na.rm = TRUE))/sd(days, na.rm = TRUE)
edge <- (edge - mean(edge, na.rm = TRUE))/sd(edge, na.rm = TRUE)
elevation <- (elevation - mean(elevation, na.rm = TRUE))/sd(elevation, na.rm = TRUE)
density <- as.numeric(Data %>% group_by(park) %>% distinct(density) %>% select(density) %>% .$density)

#----------#
#-Rainfall-#
#----------#

#Fix some overlap 3 months...
#Number of months sampled
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}


Cov <- Cov %>% mutate(
  n_mo = elapsed_months(end_date, start_date),
  str_mo = format(start_date, "%Y-%m"),
  mid_mo = ifelse(n_mo == 2, format(lubridate::floor_date(start_date, unit="month") + months(1), "%Y-%m"), NA),
  end_mo = ifelse(n_mo > 0, format(end_date, "%Y-%m"), NA)
)

Geomdf <- st_as_sf(Cov)

Geomdf <- Geomdf %>% mutate(
  str_precip = NA,
  mid_precip = NA,
  end_precip = NA,
  annual_precip = 0,
)

# dates <- unique(c(Geomdf$str_mo, Geomdf$mid_mo, Geomdf$end_mo))
# dates <- dates[!is.na(dates)]
# year <- substr(dates, 1, 4)
# month <- substr(dates, 6, 7)

year <- as.character(rep(2009:2020, each = 12))
month <- as.character(rep(sprintf('%0.2d', 1:12), length(2009:2020)))
dates <- paste0(year, "-", month)
year_ID <- as.numeric(as.factor(year))

# for(i in 1:length(year)){
#   url <- paste0("https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs/chirps-v2.0.", year[i], ".", month[i], ".tif.gz")
#   httr::GET(url = url,
#             httr::write_disk(path = paste0(getwd(), "/DataFormat/CHIRPS/", year[i], month[i], ".tif.gz"), overwrite = TRUE))
# }

filenames <- list.files("~/HMSNO/DataFormat/CHIRPS/", pattern = "chirps", full.names = TRUE)
for(i in 1:length(filenames)){
  rainfall <- raster::raster(filenames[i])
  raster::values(rainfall)[raster::values(rainfall) < 0] <- NA
  Geomdf <- Geomdf %>% mutate(str_precip = ifelse(str_mo == dates[i], raster::extract(rainfall, .), str_precip),
                            mid_precip = ifelse(mid_mo == dates[i], raster::extract(rainfall, .), mid_precip),
                            end_precip = ifelse(end_mo == dates[i], raster::extract(rainfall, .), end_precip),
                            annual_precip = ifelse(yearID == year_ID[i], raster::extract(rainfall, .) + annual_precip, annual_precip))
}

Geomdf <- Geomdf %>% mutate(precip = ifelse(n_mo == 1, (as.numeric(lubridate::ceiling_date(start_date, unit = "month") - start_date)*str_precip +
                                                          as.numeric(end_date - lubridate::floor_date(end_date, unit = "month"))*end_precip)/
                                              as.numeric(end_date - start_date),
                                            ifelse(n_mo == 2,
                                                   (as.numeric(lubridate::ceiling_date(start_date, unit = "month") - start_date)*str_precip +
                                                      (as.numeric(end_date - start_date) - (as.numeric(lubridate::ceiling_date(start_date, unit = "month") - start_date) +
                                                                                            as.numeric(end_date - lubridate::floor_date(end_date, unit = "month"))))*mid_precip +
                                                      as.numeric(end_date - lubridate::floor_date(end_date, unit = "month"))*end_precip)/
                                                     as.numeric(end_date - start_date), str_precip)))


# rainfall <- raster::stack(filenames)
# # rainfall_55 <- raster::raster("~/HMSNO/DataFormat/AfroClim/mbio_rcp45_2055/bio12_rcp45_2055_ensemble-mean-v3_wc30s.tif")
# # rainfall_85 <- raster::raster("~/HMSNO/DataFormat/AfroClim/mbio_rcp45_2085/bio12_rcp45_2085_ensemble-mean-v3_wc30s.tif")
# rainfall_55m <- raster::stack(list.files("~/HMSNO/DataFormat/AfroClim/rcp45_2055/", pattern = "pr", full.names = TRUE))
# rainfall_85m <- raster::stack(list.files("~/HMSNO/DataFormat/AfroClim/rcp45_2085/", pattern = "pr", full.names = TRUE))
# 
# temp <- c(raster::extract(rainfall, Geomdf[2,]), raster::extract(rainfall_55m, Geomdf[2,]), raster::extract(rainfall_85m, Geomdf[2,]))

precip.m <- precip.a <- array(NA, dim = c(22, nsites, nparks))

for(i in 1:dim(Geomdf)[1]){
  precip.m[Geomdf$yearID[i], Geomdf$siteID[i], Geomdf$parkID[i]] <- Geomdf$precip[i]
  precip.a[Geomdf$yearID[i], Geomdf$siteID[i], Geomdf$parkID[i]] <- Geomdf$annual_precip[i]
}

Futuredf <- Geomdf %>% group_by(parkID, siteID) %>% filter(yearID == max(yearID))

mean.AFRInames.rcp45.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/monthly/", pattern = "mean", full.names = TRUE)[c(1,5:12,2:4)]
max.AFRInames.rcp45.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/monthly/", pattern = "max", full.names = TRUE)[c(1,5:12,2:4)]
min.AFRInames.rcp45.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/monthly/", pattern = "min", full.names = TRUE)[c(1,5:12,2:4)]

mean.AFRInames.rcp85.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/monthly/", pattern = "mean", full.names = TRUE)[c(1,5:12,2:4)]
max.AFRInames.rcp85.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/monthly/", pattern = "max", full.names = TRUE)[c(1,5:12,2:4)]
min.AFRInames.rcp85.m <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/monthly/", pattern = "min", full.names = TRUE)[c(1,5:12,2:4)]

mean.AFRInames.rcp45.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/annual/", pattern = "mean", full.names = TRUE)
max.AFRInames.rcp45.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/annual/", pattern = "max", full.names = TRUE)
min.AFRInames.rcp45.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP45/annual/", pattern = "min", full.names = TRUE)

mean.AFRInames.rcp85.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/annual/", pattern = "mean", full.names = TRUE)
max.AFRInames.rcp85.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/annual/", pattern = "max", full.names = TRUE)
min.AFRInames.rcp85.a <- list.files("~/HMSNO/DataFormat/AFRICLIM/RCP85/annual/", pattern = "min", full.names = TRUE)

mean.precip.rcp85.a <- mean.precip.rcp45.a <- sd.precip.rcp85.a <- sd.precip.rcp45.a <-
  mean.precip.rcp85.m <- mean.precip.rcp45.m <- sd.precip.rcp85.m <- sd.precip.rcp45.m <- array(NA, dim = c(nsites, nparks))

for(i in 1:nparks){
  for(j in 1:nsite[i]){
    tmp <- Futuredf[Futuredf$parkID == i & Futuredf$siteID == j,]
    for(t in (nend[i]+1):12){
      add <- (t - tmp$yearID)*12
      rainfall <- raster::stack(filenames[which(year_ID %in% tmp$yearID)+add])
      precip.a[t,j,i] <- sum(raster::extract(rainfall, tmp))
      if(tmp$n_mo == 0){
        rainfall <- raster::raster(filenames[which(dates %in% tmp$str_mo)+add])
        precip.m[t,j,i] <- raster::extract(rainfall, tmp)
      }
      if(tmp$n_mo == 1){
        rainfall <- raster::raster(filenames[which(dates %in% tmp$str_mo)+add])
        str_precip <- raster::extract(rainfall, tmp)
        rainfall <- raster::raster(filenames[which(dates %in% tmp$end_mo)+add])
        end_precip <- raster::extract(rainfall, tmp)
        precip.m[t,j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str_precip +
          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end_precip)/
          as.numeric(tmp$end_date - tmp$start_date)
      }
      if(tmp$n_mo == 2){
        rainfall <- raster::raster(filenames[which(dates %in% tmp$str_mo)+add])
        str_precip <- raster::extract(rainfall, tmp)
        rainfall <- raster::raster(filenames[which(dates %in% tmp$end_mo)+add])
        end_precip <- raster::extract(rainfall, tmp)
        rainfall <- raster::raster(filenames[which(dates %in% tmp$mid_mo)+add])
        mid_precip <- raster::extract(rainfall, tmp)
        precip.m[t,j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str_precip +
          (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid_precip +
          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end_precip)/
          as.numeric(tmp$end_date - tmp$start_date)
      }#end if
    }#end t
    if(tmp$n_mo == 0){
      mo <- as.numeric(substr(tmp$str_mo, 6, 7))
      mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[mo])
      max.rainfall <- raster::raster(max.AFRInames.rcp45.m[mo])
      min.rainfall <- raster::raster(min.AFRInames.rcp45.m[mo])
      mean.precip.rcp45.m[j,i] <- raster::extract(mean.rainfall, tmp)
      sd.precip.rcp45.m[j,i] <- (raster::extract(max.rainfall, tmp) - raster::extract(min.rainfall, tmp))/4
      mean.rainfall <- raster::raster(AFRInames.rcp85.m[mo])
      max.rainfall <- raster::raster(max.AFRInames.rcp85.m[mo])
      min.rainfall <- raster::raster(min.AFRInames.rcp85.m[mo])
      mean.precip.rcp85.m[j,i] <- raster::extract(mean.rainfall, tmp)
      sd.precip.rcp45.m[j,i] <- (raster::extract(max.rainfall, tmp) - raster::extract(min.rainfall, tmp))/4
    }
    if(tmp$n_mo == 1){
      str_mo <- as.numeric(substr(tmp$str_mo, 6, 7))
      end_mo <- as.numeric(substr(tmp$end_mo, 6, 7))
      
      str.mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[str_mo])
      str.max.rainfall <- raster::raster(max.AFRInames.rcp45.m[str_mo])
      str.min.rainfall <- raster::raster(min.AFRInames.rcp45.m[str_mo])
      
      end.mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[end_mo])
      end.max.rainfall <- raster::raster(max.AFRInames.rcp45.m[end_mo])
      end.min.rainfall <- raster::raster(min.AFRInames.rcp45.m[end_mo])
      
      str.mean.precip <- raster::extract(str.mean.rainfall, tmp)
      str.max.precip <- raster::extract(str.max.rainfall, tmp)
      str.min.precip <- raster::extract(str.min.rainfall, tmp)
      
      end.mean.precip <- raster::extract(str.mean.rainfall, tmp)
      end.max.precip <- raster::extract(str.max.rainfall, tmp)
      end.min.precip <- raster::extract(str.min.rainfall, tmp)
      
      mean.precip.rcp45.m[j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.mean.precip +
                          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.mean.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      max.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.max.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.max.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      min.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.min.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.min.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      sd.precip.rcp45.m[j,i] <- (max.precip - min.precip)/4
      
      
      str.mean.rainfall <- raster::raster(mean.AFRInames.rcp85.m[str_mo])
      str.max.rainfall <- raster::raster(max.AFRInames.rcp85.m[str_mo])
      str.min.rainfall <- raster::raster(min.AFRInames.rcp85.m[str_mo])
      
      end.mean.rainfall <- raster::raster(mean.AFRInames.rcp85.m[end_mo])
      end.max.rainfall <- raster::raster(max.AFRInames.rcp85.m[end_mo])
      end.min.rainfall <- raster::raster(min.AFRInames.rcp85.m[end_mo])
      
      str.mean.precip <- raster::extract(str.mean.rainfall, tmp)
      str.max.precip <- raster::extract(str.max.rainfall, tmp)
      str.min.precip <- raster::extract(str.min.rainfall, tmp)
      
      end.mean.precip <- raster::extract(end.mean.rainfall, tmp)
      end.max.precip <- raster::extract(end.max.rainfall, tmp)
      end.min.precip <- raster::extract(end.min.rainfall, tmp)
      
      mean.precip.rcp85.m[j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.mean.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.mean.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      max.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.max.precip +
                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.max.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      min.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.min.precip +
                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.min.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      sd.precip.rcp85.m[j,i] <- (max.precip - min.precip)/4
      
    }
    if(tmp$n_mo == 2){

      str_mo <- as.numeric(substr(tmp$str_mo, 6, 7))
      end_mo <- as.numeric(substr(tmp$end_mo, 6, 7))
      mid_mo <- as.numeric(substr(tmp$mid_mo, 6, 7))
      
      str.mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[str_mo])
      str.max.rainfall <- raster::raster(max.AFRInames.rcp45.m[str_mo])
      str.min.rainfall <- raster::raster(min.AFRInames.rcp45.m[str_mo])
      
      end.mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[end_mo])
      end.max.rainfall <- raster::raster(max.AFRInames.rcp45.m[end_mo])
      end.min.rainfall <- raster::raster(min.AFRInames.rcp45.m[end_mo])
      
      mid.mean.rainfall <- raster::raster(mean.AFRInames.rcp45.m[min_mo])
      mid.max.rainfall <- raster::raster(max.AFRInames.rcp45.m[min_mo])
      mid.min.rainfall <- raster::raster(min.AFRInames.rcp45.m[min_mo])
      
      str.mean.precip <- raster::extract(str.mean.rainfall, tmp)
      str.max.precip <- raster::extract(str.max.rainfall, tmp)
      str.min.precip <- raster::extract(str.min.rainfall, tmp)
      
      end.mean.precip <- raster::extract(end.mean.rainfall, tmp)
      end.max.precip <- raster::extract(end.max.rainfall, tmp)
      end.min.precip <- raster::extract(end.min.rainfall, tmp)
      
      mid.mean.precip <- raster::extract(mid.mean.rainfall, tmp)
      mid.max.precip <- raster::extract(mid.max.rainfall, tmp)
      mid.min.precip <- raster::extract(mid.min.rainfall, tmp)
    
      mean.precip.rcp45.m[j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.mean.precip +
                          (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.mean.precip +
                          as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.mean.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      max.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.max.precip +
                                     (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.max.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.max.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      min.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.min.precip +
                                     (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.min.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.min.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      sd.precip.rcp45.m[j,i] <- (max.precip - min.precip)/4
      
      str.mean.rainfall <- raster::raster(mean.AFRInames.rcp85.m[str_mo])
      str.max.rainfall <- raster::raster(max.AFRInames.rcp85.m[str_mo])
      str.min.rainfall <- raster::raster(min.AFRInames.rcp85.m[str_mo])
      
      end.mean.rainfall <- raster::raster(mean.AFRInames.rcp85.m[end_mo])
      end.max.rainfall <- raster::raster(max.AFRInames.rcp85.m[end_mo])
      end.min.rainfall <- raster::raster(min.AFRInames.rcp85.m[end_mo])
      
      mid.mean.rainfall <- raster::raster(mean.AFRInames.rcp85.m[mid_mo])
      mid.max.rainfall <- raster::raster(max.AFRInames.rcp85.m[mid_mo])
      mid.min.rainfall <- raster::raster(min.AFRInames.rcp85.m[mid_mo])
      
      str.mean.precip <- raster::extract(str.mean.rainfall, tmp)
      str.max.precip <- raster::extract(str.max.rainfall, tmp)
      str.min.precip <- raster::extract(str.min.rainfall, tmp)
      
      end.mean.precip <- raster::extract(end.mean.rainfall, tmp)
      end.max.precip <- raster::extract(end.max.rainfall, tmp)
      end.min.precip <- raster::extract(end.min.rainfall, tmp)
      
      mid.mean.precip <- raster::extract(mid.mean.rainfall, tmp)
      mid.max.precip <- raster::extract(mid.max.rainfall, tmp)
      mid.min.precip <- raster::extract(mid.min.rainfall, tmp)
      
      mean.precip.rcp85.m[j,i] <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.mean.precip +
                                     (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.mean.precip +
                                     as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.mean.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      max.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.max.precip +
                       (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.max.precip +
                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.max.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      min.precip <- (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date)*str.min.precip +
                       (as.numeric(tmp$end_date - tmp$start_date) - (as.numeric(lubridate::ceiling_date(tmp$start_date, unit = "month") - tmp$start_date) +
                                                                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))))*mid.min.precip +
                       as.numeric(tmp$end_date - lubridate::floor_date(tmp$end_date, unit = "month"))*end.min.precip)/
        as.numeric(tmp$end_date - tmp$start_date)
      sd.precip.rcp85.m[j,i] <- (max.precip - min.precip)/4
    }#end if
    
    mean.rainfall <- raster::raster(mean.AFRInames.rcp45.a)
    max.rainfall <- raster::raster(max.AFRInames.rcp45.a)
    min.rainfall <- raster::raster(min.AFRInames.rcp45.a)
    
    mean.precip <- raster::extract(mean.rainfall, tmp)
    max.precip <- raster::extract(max.rainfall, tmp)
    min.precip <- raster::extract(min.rainfall, tmp)
    
    mean.precip.rcp45.a[j,i] <- mean.precip
    sd.precip.rcp45.a[j,i] <- (max.precip - min.precip)/4
    
    mean.rainfall <- raster::raster(mean.AFRInames.rcp85.a)
    max.rainfall <- raster::raster(max.AFRInames.rcp85.a)
    min.rainfall <- raster::raster(min.AFRInames.rcp85.a)
      
    mean.precip <- raster::extract(mean.rainfall, tmp)
    max.precip <- raster::extract(max.rainfall, tmp)
    min.precip <- raster::extract(min.rainfall, tmp)
    
    mean.precip.rcp85.a[j,i] <- mean.precip
    sd.precip.rcp85.a[j,i] <- (max.precip - min.precip)/4
  }#end j
  print(i)
}#end i

#--------------#
#-Compile data-#
#--------------#

HMSNO.data <- list(y = y, days = days, elevation = elevation, edge = edge, density = density,
                   precip.a = precip.a, precip.m = precip.m, 
                   mean.precip.rcp45.a = mean.precip.rcp45.a, sd.precip.rcp45.a = sd.precip.rcp45.a,
                   mean.precip.rcp45.m = mean.precip.rcp45.m, sd.precip.rcp45.m = sd.precip.rcp45.m,
                   mean.precip.rcp85.a = mean.precip.rcp85.a, sd.precip.rcp85.a = sd.precip.rcp85.a,
                   mean.precip.rcp85.m = mean.precip.rcp85.m, sd.precip.rcp85.m = sd.precip.rcp85.m)

HMSNO.con <- list(nspecs = nspecs, parkS = parkS, parkE = parkE,
                  nsite = nsite, nreps = nreps, nstart = nstart, nend = nend)

save(HMSNO.data, file = "~/HMSNO/DataFormat/HMSNO.Adata.Rdata")

save(HMSNO.con, file = "~/HMSNO/DataFormat/HMSNO.Acon.Rdata")

