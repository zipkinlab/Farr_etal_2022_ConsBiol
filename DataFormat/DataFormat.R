#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(tidyverse)
library(reshape2)

#-----------#
#-Load Data-#
#-----------#

Data <- NULL
files <- list.files(path = "~/HMSNO/DataFormat/EditedData", pattern = "Data", full.names = TRUE)
for(i in 1:length(files)){
  sheets <- excel_sheets(files[i])
  for(j in 1:length(sheets)){
    Temp <- read_excel(files[i], sheet = sheets[j])
    Temp <- Temp %>% mutate(species = as.character(sheets[j]), park = as.factor(i)) %>% 
      select(park, `deployment ID`, Year, species, ls(Temp, pattern = "R"), density)
    Data <- bind_rows(Data, Temp)
  }
}

SiteData <- NULL
file <- list.files(path = "~/HMSNO/DataFormat/RawData", pattern = "metadata", full.names = TRUE)
sheets <- excel_sheets(file)
for(j in 1:length(sheets)){
  Temp <- read_excel(file, sheet = sheets[j])
  Temp <- Temp %>% select(`deployment ID`, Year, ls(Temp, pattern = "R"), 
                          elevation, edge, days)
  SiteData <- bind_rows(SiteData, Temp)
}

SiteData <- bind_rows(SiteData, read_excel("~/HMSNO/DataFormat/EditedData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 1) %>%
                       select(`deployment ID`, days, elevation, edge, Year)) %>%
  bind_rows(., read_excel("~/HMSNO/DataFormat/EditedData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 2) %>%
              select(`deployment ID`, days, elevation, edge, Year)) %>%
  distinct() %>%
  select(`deployment ID`, days, elevation, edge, Year)

#-------------#
#-Format Data-#
#-------------#

Data[,5:10] <- apply(Data[,5:10], 2, function(y) as.numeric(gsub("-", NA, y)))

Data <- Data %>% mutate(`deployment ID` = ifelse(species == "sylvilocutor" & Year == "2010" & park == 1,
                                                 gsub("TEAM-001", "CT", `deployment ID`), `deployment ID`),
                        `deployment ID` = ifelse(species == "sylvilocutor" & Year == "2010" & park == 1,
                                                 gsub("-2010", "", `deployment ID`), `deployment ID`)) %>% 
                arrange(park, Year, `deployment ID`)


Data <- left_join(Data, SiteData, by = c("deployment ID", "Year"))

Data <- Data %>% drop_na(days) %>%
  mutate(siteID = as.numeric(factor(`deployment ID`, levels = unique(`deployment ID`))),
         parkID = as.numeric(park),
         yearID = as.numeric(as.factor(Year)),
         specID = as.numeric(as.factor(species))) %>%
  arrange(parkID, yearID, siteID)

#--------------#
#-Extract data-#
#--------------#

spec <- Data$specID
site <- Data$siteID
yr <- Data$yearID
park <- as.numeric(Data %>% group_by(siteID) %>%
        distinct(parkID) %>% select(parkID) %>% .$parkID)

nspec <-  max(spec)
npark <- max(Data$parkID)
nsites <- max(site)
nyrs <- max(yr)
nreps <- 6
nstart <- as.numeric(Data %>% group_by(parkID, siteID) %>%
          summarize(nstart = min(yearID) - min(Data$yearID) + 1) %>%
          select(nstart) %>% .$nstart)
nend <- as.numeric(Data %>% group_by(parkID, siteID) %>%
        summarize(nend = max(yearID) - min(Data$yearID) + 1) %>%
        select(nend) %>% .$nend)
nyrstr <- as.numeric(Data %>% group_by(parkID) %>%
          summarize(nyrstr = min(yearID) - min(Data$yearID) + 1) %>%
          select(nyrstr) %>% .$nyrstr)
nyrend <- as.numeric(Data %>% group_by(parkID) %>%
          summarize(nyrend = max(yearID) - min(Data$yearID) + 1) %>%
          select(nyrend) %>% .$nyrend)
nsitestr <- as.numeric(Data %>% group_by(parkID) %>%
            summarize(nsitestr = min(siteID)) %>%
            select(nsitestr) %>% .$nsitestr)
nsiteend <- as.numeric(Data %>% group_by(parkID) %>%
            summarize(nsiteend = max(siteID)) %>%
            select(nsiteend) %>% .$nsiteend)
#SKIP specific parks for a given species...

#SKIP specific years or reps for a given site...
# as.numeric(Data %>% filter(days > 1) %>% 
#              group_by(park, `deployment ID`) %>%
#              summarize(nstart = min(Year) - min(Data$Year) + 1) %>%
#              select(nstart) %>% .$nstart)



occ <- array(NA, dim = c(nreps, nyrs, nsites, nspec))

for(s in 1:nspec){
  for(j in 1:nsites){
    #for(t in ifelse(nstart[j] == nend[j], nstart[j], nstart[j]:nend[j])){
    for(t in nstart[j]:nend[j]){
      occ[1:nreps,t,j,s] <- rep(0,nreps)
    }
  }
}
for(i in 1:dim(Data)[1]){
  for(k in 5:(nreps+5)){
    if(is.na(Data[i,k])){
      occ[k-4,yr[i],site[i],1:nspec] <- rep(NA, nspec)
    }
  }
  occ[1:6,yr[i],site[i],spec[i]] <- as.numeric(Data[i,5:10])
}

Covariatedf <- Data %>% group_by(siteID, yearID) %>%
  select(days, edge, elevation) %>% distinct(siteID, yearID, .keep_all = TRUE)

days <- edge <- elevation <- array(NA, dim = c(nyrs, nsites))

for(i in 1:dim(Covariatedf)[1]){
  days[Covariatedf$yearID[i], Covariatedf$siteID[i]] <- Covariatedf$days[i]
  edge[Covariatedf$yearID[i], Covariatedf$siteID[i]] <- Covariatedf$edge[i]
  elevation[Covariatedf$yearID[i], Covariatedf$siteID[i]] <- Covariatedf$elevation[i]
}

density <- Data %>% group_by(siteID) %>%
  select(density) %>% distinct() %>% .$density

#--------------#
#-Compile data-#
#--------------#

HMSNO.data <- list(y = occ)

HMSNO.con <- list(nspec = nspec, nsites = nsites, nstart = nstart, nend = nend, nreps = nreps,
                  npark = npark, nyrstr = nyrstr, nyrend = nyrend, nsitestr = nsitestr, nsiteend = nsiteend,
                  park = park, elevation = elevation, edge = edge, density = density, days = days)

save(HMSNO.data, file = "~/HMSNO/DataFormat/HMSNO.data.Rdata")

save(HMSNO.con, file = "~/HMSNO/DataFormat/HMSNO.con.Rdata")

