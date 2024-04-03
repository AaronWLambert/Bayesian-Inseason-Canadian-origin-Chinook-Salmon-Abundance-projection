#**********************************************************************************
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Version: 1.0 to 3.4
# Purpose: To make yearly updates to data files:
#           GSI
#           EOS Can passage
#
#   1) Read in data
#   2) Preprocess data
#   3) Call to Stan model to generate inseason projection
#
#
#NOTES:
# This script is the general working script that is used to run single iterations 
#  of the model with a Stan file for all versions up to 3.4.
#   
#   This script is the rewritten code for leave-one-out retro testing
# Next steps: 
# 
#
# Packages #########################################################################
require(tidyverse)
require(ggthemes)
require(lubridate)
require(reshape2)

# Define Workflow Paths ============================================

# set to working directory
wd <- "C:/Users/aaron/OneDrive/Desktop/ADFG Yukon Model/Yukon Chinook Bayesian Inseason Projection Model ADFG"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"stan")
dir.data <- file.path(wd,"data")
dir.R <- file.path(wd,"R")


# Import Data ###################################################################################
# Historical Canadian EOS reconstructed run size
# This is the reconstructed data from Curry 
#*** MUST UPDATE TO LAST YEARS FILE every time you do this ***

# Reconstructed Canadian-origin total passage
CAN_hist <- readRDS(file.path(dir.data,"Canadian Passage RR 21Mar23.RDS"))

# Read in historic preseason forecasts (2013 - 2021)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_12June23_eagle+harvest.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 4Apr23.RDS"))


# Realized CanadiChinook Salmon Runsize  ################################################################
# Add new row of observations for historic observed EOS Can-Chinook
head(CAN_hist)

# Year to be added
Year <- 2023

# The total reconstucted run size estimate
abund <- 15820

# str(CAN_hist)
# str(abund)
# str(Year)
# Add to DF

# Append the new year and abundance to the data.frame
CAN_hist[nrow(CAN_hist)+1,1] <- Year
CAN_hist[nrow(CAN_hist),2] <- abund

# Make sure it worked...
view(CAN_hist)

# Uncomment next line to save the updated data.frame 
# saveRDS(CAN_hist, file = file.path(dir.data,"Canadian Chinook Passage RR 2Apr24.RDS"))


# GSI Proportions ################################################################
tail(GSI_by_year,20)
 
# Fill in the relevant year, stra tum, start and end days, and the observed prop and sd
#*** Note that day 152 = June 1***

GSI_add <- data.frame(year = rep(2023,3),
                      stratum = c(1:3),
                      startday = c(158, 178, 188),  # First day of each stratum
                      endday = c(177, 187, 204),    # Last day of each stratum
                      propCan = c(0.44, 0.54, 0.48),# Proportion from 1 to 3 stratum
                      sd = c(NA,NA,NA)) # NOTE! The updated GSI is pulled from JTC report
                                        # No SD supplied. However, not used in model

# Add new df to previous years obs
GSI_by_year_updated <- rbind(GSI_by_year,GSI_add)

# Uncomment next line to save the updated data.frame
# saveRDS(GSI_by_year_updated, file = file.path(dir.data,"GSI by year unadj 4Apr24.RDS"))

# Preseason forecast ###########################################################

# The following years JTC preseason forecast
newYearPF <- data.frame("Year" = 2024,   # Year to be added
                        "mean" = 24000) # The forecast point estimate

# Bind to PF dataframe
Pf <- rbind(pf_hist,newYearPF)

# Uncomment next line to save the updated data.frame
# saveRDS(object = Pf, file = file.path(dir.data, "pf_ver3.1_4Apr24_needupdatingtonewmethod.RDS"))

