#-------------------#
#-Working Directory-#
#-------------------#

setwd("C:/Users/farrm/Documents/GitHub/HMSNO/DataFormat")

#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(dplyr)
library(reshape2)
library(jagsUI)

#-----------#
#-Load Data-#
#-----------#

Data <- NULL
files <- list.files(path = "./EditedData", pattern = "Data", full.names = TRUE)
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
file <- list.files(path = "./RawData", pattern = "metadata", full.names = TRUE)
sheets <- excel_sheets(file)
for(j in 1:length(sheets)){
  Temp <- read_excel(file, sheet = sheets[j])
  Temp <- Temp %>% select(`deployment ID`, Year, ls(Temp, pattern = "R"), 
                          elevation, edge, days)
  SiteData <-bind_rows(SiteData, Temp)
}

Data1 <- left_join(Data, SiteData, by = c("deployment ID", "Year"))
