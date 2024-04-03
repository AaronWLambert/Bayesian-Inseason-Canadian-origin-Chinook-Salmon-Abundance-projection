#===================================================================================
# Project Name: Yukon River INSEASON FORECAST - PSS and Eagle Yearly Data-Updating
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Purpose: To import PSS and Eagle Sonar passage estimates from ADFG website &
#          preprocess it into the correct format for running in the models.
#
#===================================================================================
# NOTES:
# 
# 
# Pilot Station Sonar ##############################################################
# PSS passage can be downloaded from:
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# Select Project -> Pilot Station Sonar
# Select years 1995 - last year available
# Select Species -> Chinook
# Import as xlsx file (excel)
# Save in the <data/ADFG PSS Daily Reports> file 

# Eagle Sonar ######################################################################
# Eagle Sonar passage can be downloaded from:
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# Select Project -> Eagle Escapement
# Select years 2005 - last year available
# Select Species -> Chinook
# Import as xlsx file (excel)
# Save in the <data/ADFG Eagle Daily Reports> file 

# Packages
library(tidyverse)
library(lubridate)
library(readxl)

# set to working directory (Change this to your working directory)
wd <- "C:/Users/aaron/OneDrive/Desktop/ADFG Yukon Model/Yukon Chinook Bayesian Inseason Projection Model ADFG"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Year to be added
new.year <- 2023

# Read in the file as xlsx
# Pilot Station Sonar Daily Passage
PSS <- as.data.frame(read_xlsx(file.path(dir.data,"ADFG PSS Daily Reports/Yukon Escapement Daily Final 2023.xlsx"),skip = 3))

# Eagle Sonar Daily Passage
Eagle <- as.data.frame(read_xlsx(file.path(dir.data,"ADFG Eagle Daily Reports/Yukon Escapement Daily Eagle Final 2023.xlsx"),skip = 3))

# Pilot Station data.frame in long format for Bayesian models #############################
# Add a day column with the appropriate day ranges
PSS$Day <- c(148:252)

# Convert to long format 
PSS_hist <- PSS[,2:ncol(PSS)] %>% 
  pivot_longer(cols = -c(Day),
               names_to = "Year",
               values_to = "count",
               names_transform  = list(Year=as.numeric)) %>% 
  as.data.frame() %>% 
  replace_na(list(count=0)) %>% # replace NA with zero
  arrange(Year)  # Order years from least to greatest (important for model)

# make sure it worked
head(PSS_hist)
dim(PSS_hist)

# Uncoment to save PSS_hist  
# saveRDS(object = PSS_hist, file = file.path(dir.data, "PSS passage Final 2023.RDS"))


# Eagle #############################################################################
# Convert to long format 
eagle_hist <- Eagle[,2:ncol(Eagle)] %>% 
  pivot_longer(cols = -c(Day), 
               names_to = "Year", 
               values_to = "count",
               names_transform  = list(Year=as.numeric)) %>% 
  as.data.frame()%>% 
  replace_na(list(count=0)) %>% # replace NA with zero
  arrange(Year)


# make sure it worked
head(eagle_hist)
dim(eagle_hist)

# Uncomment to save updated dataframe
# saveRDS(object = eagle_hist, file = file.path(dir.data, "Eagle passage Final 2023.RDS"))
