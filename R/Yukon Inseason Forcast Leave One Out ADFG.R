#**********************************************************************************
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 3.26.24

# Version: ADF&G 
# Purpose: To generate inseason Yukon River Canadian origin Chinook 
#           salmon abundance projections
#
#   1) Read in data
#   2) Preprocess data
#   3) Call to Stan model to generate inseason projection
#   4) Assorted plots
#
#NOTES:
#   This script is the rewritten code for leave-one-out retro testing


# Next steps: 
# 
#
# Packages #########################################################################
require(rstan)
require(bayesplot)
require(tidyverse)
require(mgcv)
require(ggthemes)
require(viridis)
require(shinystan)
require(lubridate)
require(reshape2)
require(dplyr)
require(tidybayes)


# Parralize for optimum model run time (Speeds things up...)
rstan_options(auto_write = TRUE) 
mc.cores = parallel::detectCores()

# Define Workflow Paths ============================================

# set to working directory (Set to directory with all the necessary folders)
wd <- "C:/Users/aaron/OneDrive/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
# To use this, one must have folders with these names in the working directory
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")


######### Import Data ###############
# Total Canadian-origin reconstructed run size for each year
CAN_hist <- readRDS(file.path(dir.data,"Canadian Passage RR 21Mar23.RDS"))

# PSS historical for Jun 1st = 152, to Aug 31 = 243
# Data obtained from ADFG website
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# NOTE!!!! See use "Dataframe preprocess.R" to get  most recent version.
#  If not using Dataframe preprocess.R, uncomment out the following to load relevant PSS historical data
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2022.RDS"))

# Eagle passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))

# Read in historic preseason forecasts (2007 - 2022)
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_14Aprl23_eagle+harvest.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 21Mar23.RDS"))

# Read in PSS observation error estimates ( Not currently used in any models)
PSS_sd <- readRDS(file = file.path(dir.data,"PSS SD 1995_2021.RDS"))

# Midpoint from cum dist fitting for individual years
logistic.all <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))

# Logistic curve parameters from fitting to each year (Used as informative prior)
logistic.mean <- readRDS(file = file.path(dir.data,"fit logistic for midpoint 18Apr23.RDS"))

# Normal dist parameters from fitting to each year (Used as informative prior)
normal.all <- read.csv(file = file.path(dir.output,"normal curve parameters All Chinook 1995_2022.csv"))

# Control Section #############################################################################################
# This is where users can input dates, stan model, and stan model controls.

# Model Version
model.version <- "PSSprop_ESprop"

# Range of years 1995 to current year
myYear <- 2019

# Range of days 152 -243
myDay <- 153

# MCMC Parameters
n.chains <- 4   # Number of chains to run
n.iter <- 3000  # Number of iterations to run
n.thin <- 2     # How many iterations to thin

# Start Years for Predictors
startYearPF <- 2007    # Preseason forecast
startYearGSI <- 2005   # Genetic stock identification
startYearPSS <- 2005   # Pilot station sonar
startYearEagle <- 2005 # Eagle Sonar

# Wont typically change;
#  *Day 152 is June 1*
startDayPSS <- 148
endDayPSS <- 250

#********************************* Prepare all data for Stan model *******************
# Preseason Forecast Prio Dist ######################

# Current Preseason forecast Can Origin
pf <- log(pf_hist$mean[pf_hist$Year == myYear])

# Past preseason forecasts for to compute sd for prior in stan model
PF_vect <- pf_hist$mean[pf_hist$Year != myYear]

# Names for accounting
names(PF_vect) <-  pf_hist$Year[pf_hist$Year != myYear]

# EOS reconstructed runsize for past years that have a preseason forecast 
EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 & 
                           CAN_hist$Year != myYear]

names(EOS) <- CAN_hist$Year[CAN_hist$Year >= 2007 &
                              CAN_hist$Year != myYear]

# SD for preseason forecast prior
pf_sigma <- sd(log(PF_vect/EOS))

# Years in PF
yearPF <- unique(pf_hist$Year[pf_hist$Year != myYear])

n_yearsPF <- length(yearPF)

# Reconstructed Can-origin Run Size ############################################

# Vector of run sizes excluding the year of interest
totalEOS <- CAN_hist$can.mean[ CAN_hist$Year >= startYearPSS &
                                 CAN_hist$Year != 1996 &
                                 CAN_hist$Year != myYear]

# Name the elements for accounting purposes
names(totalEOS) <- CAN_hist$Year[ CAN_hist$Year >= startYearPSS &
                                    CAN_hist$Year != 1996 &
                                    CAN_hist$Year != myYear]

# Number of years included in EOS
n_totalEOS <- length(totalEOS)

# Pilot Station Sonar Data ##############################################################

# Vector of historic PSS years excluding myYear
yearPSS <- unique(PSS_hist$Year[PSS_hist$Year != myYear &
                                  PSS_hist$Year >= startYearPSS &
                                  PSS_hist$Year != 1996])

# Number of years used in model
n_yearPSS <- length(yearPSS)

# # PSS days included up to myDay
# #  I.e., June 1 = 1, June 2 = 2 ....
dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & 
                                PSS_hist$Day>= startDayPSS])

# All days used from beginning of season to end of season
dayPSS_all <- c(startDayPSS:endDayPSS)

# Pointer vector for myDay
loc.allDays.myDay <- which(dayPSS_all %in% myDay)

# Number of days used
n_dayPSS <- length(dayPSS)

# Number of all days in season
n_dayPSS_all <- length(startDayPSS:endDayPSS )

# Current years PSS passage up to myDay
curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear &
                            PSS_hist$Day <= myDay &
                            PSS_hist$Day >= startDayPSS])

# Name the days for accounting
names(curr_PSS) <- (PSS_hist$Day[PSS_hist$Year == myYear &
                                 PSS_hist$Day <= myDay &
                                 PSS_hist$Day>=startDayPSS])

# Number of days used
n_curr_PSS <- length(curr_PSS)

# Cumulative current PSS counts
cum_curr_PSS <- cumsum(curr_PSS)

# Create matrix with PSS daily Chinook passage from day 152 to myday 
# Create matrix with dimensions days x historic Year
PSS_mat <- matrix(nrow = n_dayPSS,
                  #start 1996 to account for missing year 1996
                  ncol = n_yearPSS) 

# Give names to matrix for accounting
colnames(PSS_mat) <- yearPSS

rownames(PSS_mat) <- dayPSS

# Create vector of counts from relevant years and days to input into matrix
(count_vect <- PSS_hist$count[PSS_hist$Year != myYear &
                                PSS_hist$Year >= startYearPSS &
                                PSS_hist$Year <= max(CAN_hist$Year) &
                                PSS_hist$Day <= (myDay) &
                                PSS_hist$Day >= (startDayPSS)])

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:n_yearPSS){
  for (d in 1:n_dayPSS) {
    
    PSS_mat[d,y] <- count_vect[counter]
    counter = counter+1
    
  }
}

# Cumulative PSS matrix
cum_PSS_mat <-  apply(PSS_mat, 2, cumsum)

# Create matrix of PSS for every day in the season
PSS_mat_all <- matrix(nrow = n_dayPSS_all,
                      #start 1996 to account for missing year 1996
                      ncol = n_yearPSS) 

# Give names to matrix for accounting
colnames(PSS_mat_all) <- yearPSS

rownames(PSS_mat_all) <- c(startDayPSS:endDayPSS)

# Create vector of counts from relevant years and days to input into matrix
(count_vect_all <- PSS_hist$count[PSS_hist$Year != myYear &
                                    PSS_hist$Year >= startYearPSS &
                                    PSS_hist$Day <= endDayPSS &
                                    PSS_hist$Day >= startDayPSS 
                                  
])

# The number of observations in count_vector for PSS historic counts
(n_hist_counts_all <- length(count_vect_all))

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:n_yearPSS){
  for (d in 1:n_dayPSS_all) {
    
    PSS_mat_all[d,y] <- count_vect_all[counter]
    counter = counter+1
    
  }
}

# Matrix of cumulative pss pasage by day and year
cum_PSS_mat_all <- apply(PSS_mat_all,MARGIN = 2, FUN = cumsum)

# Copy matrix
PSS_mat_prop_all <- PSS_mat_all

# Populate with NA's
PSS_mat_prop_all[,] <- NA

# Fill with observed proportions by day and year
for (y in 1:n_yearPSS){
  for (d in 1:n_dayPSS_all) {
    
    
    PSS_mat_prop_all[d,y] <- cum_PSS_mat_all[d,y]/ cum_PSS_mat_all[n_dayPSS_all,y]
    
    
  }
}

# Total PSS passage by year up to myDay
cumPSS <- apply(PSS_mat, 2, sum)

# Calculate the proportion of PSS passage up to myDay of the EOS total Canadian Chinook 
obs.prop.PSS <- cumPSS/totalEOS

# Mean proportion on myDay at PSS
mean_propPSS <- mean(obs.prop.PSS)

# Eagle Sonar #####################################################################

# Vector of historic Eagle Sonar passage excluding myYear
yearEagle <- unique(Eagle_hist$Year[Eagle_hist$Year != myYear])

# Number of Eagle sonar years used in model
n_yearEagle <- length(yearEagle)

# Eagle days included up to myDay
dayEagle <- unique(Eagle_hist$Day[Eagle_hist$Day <= myDay & 
                                    Eagle_hist$Day>= startDayPSS])

# Number of days used
n_dayEagle <- length(dayEagle)

# Eagle daily passage estimate for days up to myDay for the year of interest (myYear)
curr_Eagle_vect <- (Eagle_hist$count[Eagle_hist$Year == myYear &
                                       Eagle_hist$Day <= myDay & 
                                       Eagle_hist$Day>=startDayPSS])

curr_Eagle <- array(data = curr_Eagle_vect,
                    dim = length(curr_Eagle_vect))

# Number of days used
n_curr_Eagle <- length(curr_Eagle)

# Cumulative Eagle counts for myYear
cum_curr_Eagle <- cumsum(curr_Eagle)

# Location pointers to index years in common for variance calculations
# 
loc.eagle.years <- which(yearPSS %in% yearEagle)

loc.PF.years.eagle <- which(yearEagle %in% yearPF)

loc.PF.years.PSS <- which(yearPSS %in% yearPF)

# Make sure it works....
# length(loc.eagle.years)
# 
# yearPSS[loc.eagle.years]


# Create matrix with dimensions myday and myYear
Eagle_mat <- matrix(nrow = n_dayPSS,
                    #start 1996 to account for missing year 1996
                    ncol = n_yearEagle) 

# Give names to matrix for accounting
colnames(Eagle_mat) <- yearEagle

rownames(Eagle_mat) <- dayPSS

# Create vector of counts from relevant years and days to input into matrix
# (count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) & *******************************************
#                                 PSS_hist$Day <= (myDay)])
(count_vect_Eagle <- Eagle_hist$count[Eagle_hist$Year != myYear &
                                        Eagle_hist$Day <= (myDay) &
                                        Eagle_hist$Day >= (startDayPSS)])

# The number of observations in count_vector for PSS historic counts
n_hist_Eaglecounts <- length(count_vect_Eagle)

# Use loop to populate matrix with counts by days
# Set counter

# This is to index where to start filling the Eagle matrix 
SE <- n_dayPSS- n_dayEagle + 1

#Counter for loop
counter <-1

# Loop to fill matix by year and day 
# !!! Ignore sucscript out of bounds warning!!!!!
for (y in 1:n_yearEagle){
  for (d in SE:n_dayPSS) {
    
    
    Eagle_mat[d,y] <- count_vect_Eagle[counter]
    
    counter = counter+1
    
  }
}

# Fill in the rest of the matrix with zeros
Eagle_mat[is.na(Eagle_mat)]<-0

# Cumulative count matrix for historic years
Eagle_cum_hist_mat <- apply(X = Eagle_mat,MARGIN = 2,FUN = cumsum)

# Cumulative Counts by year for plotting
cumEagle <- apply(Eagle_mat,MARGIN = 2,sum)

names(cumEagle) <- yearEagle

# Calculate the proportion for day D for each year Y
prop_Eagle <- cumEagle/totalEOS[loc.eagle.years]

# Mean cumulative proportion at eagle on myday
mean_prop_Eagle <- mean(prop_Eagle)

# Vector containing avg GSI proportions for each day ACROSS ALL YEARS (PSS_GSI) ####
meanGSI_vect <- vector(length = length(startDayPSS:myDay))

names(meanGSI_vect) <- c(startDayPSS:myDay)

# Vector of GSI sd
sdGSI_vect <- vector(length = length(startDayPSS:myDay))

names(sdGSI_vect) <- c(startDayPSS:myDay)

counter <- 1

for (d in startDayPSS:myDay) {
  # d = 175
  meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d &
                                                     GSI_by_year$endday >=d &
                                                     GSI_by_year$year != myYear &
                                                     GSI_by_year$year != 2013])
  
  sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d &
                                                 GSI_by_year$endday >=d &
                                                 GSI_by_year$year != myYear&
                                                 GSI_by_year$year != 2013])
  
  counter <- counter+1
}

# No GSI data for day 148 and 149. Use day 150 for these values
meanGSI_vect[1:2] <- meanGSI_vect[3]

sdGSI_vect[1:2] <- sdGSI_vect[3]

# Mean GSI by mean strata dates. This not currently used #######################
meanStartDay <- GSI_by_year %>% 
  filter(year!= myYear & year != 2013) %>%
  group_by(stratum) %>% 
  summarise("meanStartDay" = mean(startday)) %>% 
  as.data.frame()

GSI_mean_by_strata <- GSI_by_year %>% 
  filter(year!= myYear & year != 2013) %>%
  summarize("stratumMean" = c(mean(propCan[stratum == 1]),
                              mean(propCan[stratum == 2]),
                              mean(propCan[stratum == 3 | stratum == 4])),
            "stratumSD" = c(sd(propCan[stratum == 1]),
                            sd(propCan[stratum == 2]),
                            sd(propCan[stratum == 3 | stratum == 4]))) %>%  as.data.frame()

GSI <- cbind(GSI_mean_by_strata,round(meanStartDay[1:3,]))

GSI_avg <-c(rep(GSI$stratumMean[GSI$stratum == 1], 
                times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
            rep(GSI$stratumMean[GSI$stratum == 2], 
                times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
            rep(GSI$stratumMean[GSI$stratum == 3],
                times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))

GSI_avg_vect <- GSI_avg[1:n_dayPSS]

GSI_sd <- c(rep(GSI$stratumSD[GSI$stratum == 1], 
                times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
            rep(GSI$stratumSD[GSI$stratum == 2], 
                times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
            rep(GSI$stratumSD[GSI$stratum == 3],
                times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))

GSI_sd_vect <- GSI_sd[1:n_dayPSS]

# Adjusted PSS by avg gsi for testing (Not currently used)#####################################
# Create matrix with dimensions myday and myYear
PSS_mat_adj <- matrix(nrow = n_dayPSS,
                      ncol = n_yearPSS) 

# Give names to matrix
colnames(PSS_mat_adj) <- yearPSS

rownames(PSS_mat_adj) <- dayPSS

for (y in 1:n_yearPSS) {
  for (d in 1:n_dayPSS) {
    # y = 1
    # d = 15
    PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
    
  }
  
}

# Cumulative GSI adjusted historic counts
cumPSS_adj <- apply(PSS_mat_adj, 2, sum)
adj_curr_PSS <- vector(length = n_dayPSS)

# GSI adjusted for this year
for (d in 1:n_dayPSS) {
  adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
}

# GSI Beta parameters
meanpropCAN <- meanGSI_vect

sd_meanpropCAN <- sdGSI_vect

paramA <- vector(length = n_dayPSS)
paramB <- vector(length = n_dayPSS)

for(x in 1:n_dayPSS){
  # x = 1
  paramA[x] <- ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;
  
  paramB[x] <- paramA[x] * ((1/meanpropCAN[x])-1);
  
}

# PSS var calculations for model incorporating uncertainty (Not currently used) ##########################
tempSD <- PSS_sd[PSS_sd$Day<= myDay &
                   PSS_sd$Year != myYear,] %>%
  group_by(Year) %>%
  summarise(sd = sqrt(sum(Var)))


full_tempSD <- data.frame("Year" = yearPSS, "sd" = 0)

full_tempSD$sd <- tempSD$sd[match(full_tempSD$Year,tempSD$Year)]

full_tempSD$sd[ is.na(full_tempSD$sd)] <- 0

PSS_year_sd <- full_tempSD$sd

names(PSS_year_sd) <- yearPSS

curr_PSS_year_sd <- sqrt(sum(PSS_sd$Var[PSS_sd$Day<= myDay &
                                          PSS_sd$Year == myYear] ))

# Normal and Logistic parameters for PSSnormal and PSSlogistic models
# Store each parameter as a vector for model input

# Normal dist arrival parameters
(ps_mu <- normal.all$mid[normal.all$year != myYear &
                           normal.all$year >= startYearPSS])

(ps_sd <- normal.all$sd[normal.all$year != myYear &
                          normal.all$year >= startYearPSS])

(ps_alpha_norm <- normal.all$alpha[normal.all$year != myYear&
                                     normal.all$year >= startYearPSS])

n_ps_mu <- length(ps_mu)
n_ps_sd <- length(ps_sd)
n_ps_alpha <- length(ps_alpha_norm)

# Logistic curve priors
ps_m <- logistic.all$mid[logistic.all$year != myYear &
                           logistic.all$year >= startYearPSS]

ps_s <- logistic.all$sd[logistic.all$year != myYear &
                          logistic.all$year >= startYearPSS]

ps_alpha.log <- logistic.all$alpha[logistic.all$year != myYear &
                                     logistic.all$year >= startYearPSS]

n_ps_m <- length(ps_m)
n_ps_s <- length(ps_s)
n_ps_alpha_log <- length(ps_alpha.log)

# Stan Model Call ######################################

# Inits list for PSSreg_GSI
# inits <- function(){
#   list(
#         sigma = runif(1,0,2),
#         alpha = runif(1,150000,200000),
#         beta = runif(1,0,0.5)
#         # ln_phi = runif(1,-1,1),
#         # propCAN_logit = runif(n_dayPSS,0,1),
#         # mid = runif(1,174,176),
#         # shape = runif(1, 4,6),
#         # scale = runif(1,140000,160000),
#         # alpha_sst = runif(1,1,2),
#         # beta_sst = runif(1,-2,-1),
#         # phi = runif(1,1,50)
# 
#   )
#   }

# Inits list for model PSSnorm_ESprop
# inits<-  function(){
#   list(
#     # "ps_alpha_curr" = runif(1,15e3,17e3),
#     "ps_alpha_curr" = runif(1,10e4,11e4),
#     "ps_mu_curr" = runif(1,172,175),
#     "ps_sd_curr" = runif(1,5,7),
#     "ps_alpha_hist" = runif(n_yearPSS,10e4,11e4),
#     "ps_mu_hist" = runif(n_yearPSS,172,175),
#     "ps_sd_hist" = runif(n_yearPSS,5,7),
#     "sigma" = runif(1,2,3),
#     "sigma_hist" = runif(n_yearPSS,2,3),
#     "beta" = runif(1,0,0.5),
#     "sigma_reg" = runif(1,0.2,0.5),
#     "alpha" = runif(1,10,11),
#     
#   )
# }

# Inits list for PSSlogistic_ESprop
# inits <-function(){
#   
#   list("ps_alpha_curr" = runif(1,1e5,11e4),
#        "ps_mid_curr" = runif(1,168,170),
#        "ps_shape_curr" = runif(1,2,3),
#        "ps_shape_hist" = runif(n_yearPSS,2,3),
#        "ps_mid_hist" = runif(n_yearPSS,168,170),
#        "ps_alpha_hist"=runif(n_yearPSS,1e5,11e4),
#        "sigma" = runif(1,0.1,0.5),
#        "sigma_hist" = runif(n_yearPSS,0,1),
#        "beta" = runif(1,0.3,0.5),
#        "sigma_reg" = runif(1,0.2,0.5),
#        "alpha" = runif(1,100,200),
#        "beta_sst" = runif(1,-2,0),
#        "alpha_t" = runif(1,170,180)
#        
#   )
# }

inits_ll <- list(inits(), inits(), inits(), inits())

# Run the stan model
fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ",
                                            model.version ,".stan", sep = "")), #model.version from above
            data = list("PSS_mat"=PSS_mat,
                        "PSS_mat_all" = PSS_mat_all,
                        "cum_PSS_mat_all" = cum_PSS_mat_all,
                        "n_totalEOS"=n_totalEOS,
                        "totalEOS"=totalEOS, 
                        "n_dayPSS"=n_dayPSS,
                        "dayPSS"=dayPSS, 
                        "n_yearPSS"=n_yearPSS,
                        "yearPSS"=yearPSS,
                        "n_curr_PSS"=n_curr_PSS,
                        "curr_PSS"=curr_PSS,
                        "Pf" = pf, 
                        "Pf_sigma" = pf_sigma,
                        "meanpropCAN" = meanGSI_vect,
                        "sd_meanpropCAN" = sdGSI_vect,
                        # "meanpropCAN" = GSI_avg_vect,
                        # "sd_meanpropCAN" = GSI_sd_vect,
                        # "mean_adj_PSS_mat" = PSS_mat_adj,
                        # "mean_adj_curr_PSS"= adj_curr_PSS,
                        # "PSS_cum_hist_mat" = PSS_cum_hist_mat,
                        # "cum_curr_PSS" = cum_curr_PSS,
                        # "meanpropCAN"= GSI_avg_vect,
                        # "sd_meanpropCAN" = GSI_sd_vect,
                        # "GSI_mean"= GSI_avg_vect,
                        # "GSI_sd" = GSI_sd_vect,
                        "Eagle_mat" = Eagle_mat,
                        "n_totalEOS" = n_totalEOS,
                        "totalEOS" = totalEOS, 
                        "n_dayEagle" = n_dayEagle,
                        "dayEagle" = dayEagle, 
                        "n_yearEagle" = n_yearEagle,
                        "yearEagle" = yearEagle,
                        "n_curr_Eagle" = n_curr_Eagle,
                        "curr_Eagle" = curr_Eagle,
                        "loc_eagle_years" = loc.eagle.years,
                        # "loc_pf_years" = loc.PF.years,
                        "n_yearsPF" = n_yearsPF,
                        "loc_pf_years_Eagle" = loc.PF.years.eagle,
                        "loc_pf_years_PSS" = loc.PF.years.PSS,
                        "paramA" = paramA,
                        "paramB" = paramB,
                        "curr_PSS_year_sd" = curr_PSS_year_sd,
                        "PSS_year_sd" = PSS_year_sd,
                        "propEagle" = prop_Eagle,
                        "mean_propEagle" = mean_prop_Eagle,
                        # "may_sst_hist" = may.sst.hist,
                        # "n_yearSST" = n_yearSST,
                        # "may_sst_curr" = may.sst.curr,
                        # "dev_hist" = deviations.mp,
                        # "hist_MP" = hist.midpoint,
                        "myDay" = myDay,
                        # "avg_mid" = log_mid,
                        # "shape" = log_s,
                        "n_dayPSS_all" = n_dayPSS_all,
                        "PSS_mat_prop_all" = PSS_mat_prop_all,
                        "dayPSS_all" = dayPSS_all,
                        "mean_propPSS" = mean_propPSS,
                        "ps_mu"= ps_mu,
                        "ps_sd"= ps_sd,
                        "ps_alpha_norm"= ps_alpha_norm,
                        "n_ps_mu" = n_ps_mu,
                        "n_ps_sd" = n_ps_sd,
                        "n_ps_alpha" = n_ps_alpha,
                        "loc_allDays_myDay" = loc.allDays.myDay,
                        # "n_sst" = n_sst,
                        # "loc_devMP_allYears" = loc.devMP.allYears,
                        # "PSS_mat_all" = PSS_mat_all,
                        "ps_alpha_log"= ps_alpha.log,
                        "ps_m"= ps_m,
                        "ps_s"= ps_s,
                        "n_ps_m" = n_ps_m,
                        "n_ps_s" = n_ps_s,
                        "n_ps_alpha_log" = n_ps_alpha_log,
                        "cum_PSS_mat" = cum_PSS_mat,
                        "cum_curr_PSS" = cum_curr_PSS 
                        
            ),
            # init = inits_ll, # Inits list. Uncomment to use
            chains = n.chains,
            iter = n.iter, 
            thin = n.thin, 
            cores = mc.cores,
            # control = list(max_treedepth = 25, adapt_delta = 0.99), # adjust this for treedepth warnings/div chains warnings
            verbose = F,
            save_warmup = T
            
)


# Figures Section #####################################################

# Trace plots to check for convergence
traceplot(object = fit, c(
  "alpha",
  "beta",
  "sigma"
  # "ps_alpha_hist"
  # "ps_mu_hist"
  # "ps_sd_hist"
  # "alpha_sst",
  # "sigma_predPSS",
  # "sigma_predEagle",
  # "RunSize"
))


# Launch shiny app (Eady way to look at diagnostic plots)
# shinystan::launch_shinystan(as.shinystan(fit)) # Uncomment to use

# Extract parameter estimates for plotting & analysis
pars <- rstan::extract(fit)

# Note that the following figures only work for some models

# Regression plot (PSSreg_ESprop) ##############################################

# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$predPSS, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

# Get the quantiles for predPSS from model output
# quant.predPSS2 <- apply(X = pars$cum_predicted_histPSS, 
#                         MARGIN = 2, 
#                         FUN = quantile, 
#                         probs=c(0.025, 0.25, 0.5, 0.75, 0.975))


# Get order of cumPSS for polygon below
ord <- order(cumPSS)

# Regression plot for non-GSI-adjusted PSS passage
plot(x = cumPSS/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = expression(Sigma ~ PSS[D] ~("1000's"~ Chinook ~Salmon)),
     ylab = "Total EOS Canadian Chinook Salmom (1000's)",
     xlim = c(0,(max(cumPSS))/1000),
     ylim = c(0, (max(totalEOS)/1000)+10),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(cumPSS[ord]/1000, rev(cumPSS[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS[ord]/1000, rev(cumPSS[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)
points(x = pars$cum_predicted_histPSS[1,]/1000, y = totalEOS/1000, cex = 3)

# Add years to plot
text(x = cumPSS/1000, y = (totalEOS/1000)+8, labels = yearPSS)
points(x = sum(PSS_hist$count[PSS_hist$Day <= myDay &
                                PSS_hist$Year == myYear])/1000, 
       y = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
       cex = 3,
       lwd = 5,
       col = "Blue",
       pch = 23,
       bg = "red"
)
text(x = sum(PSS_hist$count[PSS_hist$Day <= myDay &
                              PSS_hist$Year == myYear])/1000, 
     y = (CAN_hist$can.mean[CAN_hist$Year == myYear]/1000)+10,
     labels = myYear,
     col = "Blue",
     cex = 2)

# Regression plot for GSI adjusted passage (PSSreg_GSI) ########################
# Order for polygon
ord <- order(cumPSS_adj)

# Regression plot
plot(x = cumPSS_adj/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     # xlab = paste("\sigma PSS Passage (1000's Chinook Salmon), 1995-",
     #              myYear-1,")" ),
     xlab = expression(Sigma ~ PSS[D] ~("1000's"~ Chinook ~Salmon)),
     ylab = "Total EOS Canadian Chinook Salmom (1000's)",
     # main = paste("Predicted PSS Fit","1995 to",
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS_adj)/1000),
     ylim = c(min(quant.predPSS/1000, totalEOS/1000), max(quant.predPSS/1000)+10),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(cumPSS_adj[ord]/1000, rev(cumPSS_adj[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS_adj[ord]/1000, rev(cumPSS_adj[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS_adj[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)
text(x = cumPSS_adj/1000, y = (totalEOS/1000)+10, labels = yearPSS)
points(x = max(cum_curr_PSS)/1000* mean(GSI_avg), 
       y = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
       col = "green",
       pch = 25,
       bg = "blue",
       lwd = 2,
       cex = 2
)
text(x = max(cum_curr_PSS)/1000* mean(GSI_avg)+12, 
     y = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
     labels = myYear,
     cex = 2)

# Regression for EOS PSS vs Can EOS (PSSnormal_ESprop)
# Get order of cumPSS for polygon below
cumPSS_all <- apply(PSS_mat_all, 2, sum)

ord <- order(cumPSS_all)

# Regression plot for non-GSI-adjusted PSS passage
plot(x = cumPSS_all/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Total EOS PSS Chinook Salmon (1000's)" ),
     ylab = ("Total EOS Canadian Chinook Salmon (1000's)"),
     # ylab = paste("1000's Canadian Chinook Salmon (1995-",myYear-1,")"),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     xlim = c(min(cumPSS_all/1000)-20,(max(cumPSS_all/1000)+20)),
     ylim = c(min(quant.predPSS/1000)-30, max(quant.predPSS/1000)+30),
     # ylim = c(0,350),
     cex.axis = 1.2,
     cex.lab = 1.2)
polygon(x = c(cumPSS_all[ord]/1000, rev(cumPSS_all[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS_all[ord]/1000, rev(cumPSS_all[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS_all[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)

# Add years to plot
text(x = cumPSS_all/1000, y = (totalEOS/1000)+5, labels = yearPSS, cex = .75)
points(x = sum(PSS_hist$count[
  PSS_hist$Year == myYear])/1000, 
  y = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
  cex = 3,
  lwd = 5,
  col = "Blue",
  pch = 23,
  bg = "red"
)
text(x = sum(PSS_hist$count[
  PSS_hist$Year == myYear])/1000, 
  y = (CAN_hist$can.mean[CAN_hist$Year == myYear]/1000)+10,
  labels = myYear,
  col = "Blue",
  cex = 1)

# Eagle regression plot (PSSreg_ESreg)######################################################

# Get the quantiles for predPSS from model output
quant.predEagle <- apply(X = pars$predEagle, 
                         MARGIN = 2, 
                         FUN = quantile, 
                         probs=c(0.025, 0.25, 0.5, 0.75, 0.975))


# Get order of cumPSS for polygon below
ord <- order(cumEagle)

# Regression plot for Eagle Sonar
plot(x = cumEagle/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Eagle Cumulative Chinook Salmon (1000's; Day = ",myDay,")" ),
     ylab = paste("EOS Eagle + Harvest (1000's Chinook Salmon)"),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     # xlim = c(0,4),
     ylim = c(min(quant.predEagle/1000)-10, max(quant.predEagle/1000)+10),
     # ylim = c(0,150),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(cumEagle[ord]/1000, rev(cumEagle[ord]/1000)),
        y = c(quant.predEagle[1,ord]/1000,rev(quant.predEagle[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumEagle[ord]/1000, rev(cumEagle[ord]/1000)),
        y = c(quant.predEagle[2,ord]/1000,rev(quant.predEagle[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumEagle[ord]/1000, quant.predEagle[3,ord]/1000, col = "red", lw = 2)
text(x = cumEagle/1000, y = (totalEOS/1000) + 3, labels = yearEagle)
points(x = sum(Eagle_hist$count[Eagle_hist$Day <= myDay &
                                  Eagle_hist$Year == myYear])/1000, 
       y = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
       cex = 3,
       lwd = 5,
       col = "Blue",
       pch = 23,
       bg = "red"
)

text(x = sum(Eagle_hist$count[Eagle_hist$Day <= myDay &
                                Eagle_hist$Year == myYear])/1000, 
     y = median(pars$curr_predEagle)/1000 +13,
     labels = myYear,
     cex = 2,
     col = "blue")
abline(0,1, lty = 2, lwd = 5)

median(pars$curr_predEagle)

# Normal Curve Fit Plots PSSnormal_ESprop ##############################################

# Loop through years to calculate quantiles and put into data frame for plotting
for(y in 1:n_yearPSS) {
  
  # y = 1
  year <- yearPSS[y]
  
  quant.predPSS <- apply(X = pars$ps_pred_hist[,,y], 
                         MARGIN = 2, 
                         FUN = quantile, 
                         probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  
  # If statement to plot without retrospective testing
  #  When Retrospective = FALSE, then the realized Canadian-origin line is omitted
  
  
  pred.plot.df <- data.frame(dayPSS_all,
                             PSS_hist$count[PSS_hist$Year == year & 
                                              PSS_hist$Day <= endDayPSS],
                             t(quant.predPSS))
  names(pred.plot.df) <- c("Days","count","low95","low50","median",
                           "up50","up95")
  
  pred.plot.df$Obs <- ifelse(pred.plot.df$Days <= myDay, "Obs", "Realized")
  
  pred.plot.df$Year <- year
  
  ifelse(y==1,masterDF <- pred.plot.df, masterDF <-rbind(masterDF,pred.plot.df))
}

# Quantiles for the current year myYear 
# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$ps_pred_curr, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))



pred.plot.df <- data.frame(dayPSS_all,
                           PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= endDayPSS],
                           t(quant.predPSS))
names(pred.plot.df) <- c("Days","count","low95","low50","median",
                         "up50","up95")

pred.plot.df$Obs <- ifelse(pred.plot.df$Days <= myDay, "Obs", "Realized")

pred.plot.df$Year <- myYear


# Bind into a data.frame
masterDF <- rbind(masterDF, pred.plot.df)


# Plot results
ggplot(masterDF, aes(x = Days,  y = count/1000))+
  geom_point(aes(color = Obs)) +
  # geom_col(alpha = .3)+
  geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.3) +
  geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.3) +
  geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
  coord_cartesian(x = c(148,215))+
  labs(x = "",
       y = "")+
  # theme_linedraw(base_line_size = 0)+
  # theme(text = element_text(size = 20, family = "serif"),
  #       panel.grid.major =element_line(),    #strip major gridlines
  #       panel.grid.minor = element_blank(),    #strip minor gridlines
  #       axis.ticks = element_blank(),
  #       panel.background = element_blank(),
  #       axis.line = element_line(color = "black"))+
  
  scale_fill_colorblind(name = "")+
  scale_color_colorblind(name = "")+
  facet_wrap(~Year, nrow = 6, scales = "free_y")


# Logistic Curve  model fit PSSlogistic_ESprop ##########################

# Loop through years to calculate quantiles and put into data frame for plotting
for(y in 1:n_yearPSS) {
  
  # y = 1
  year <- yearPSS[y]
  
  quant.predPSS <- apply(X = pars$ps_pred_hist[,,y], 
                         MARGIN = 2, 
                         FUN = quantile, 
                         probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  
  # If statement to plot without retrospective testing
  #  When Retrospective = FALSE, then the realized Canadian-origin line is omitted
  
  pred.plot.df <- data.frame(dayPSS_all,
                             cumsum(PSS_hist$count[PSS_hist$Year == year & PSS_hist$Day <= endDayPSS]),
                             t(quant.predPSS))
  names(pred.plot.df) <- c("Days","cumcount","low95","low50","median",
                           "up50","up95")
  
  pred.plot.df$Obs <- ifelse(pred.plot.df$Days <= myDay, "Obs", "Realized")
  
  pred.plot.df$Year <- year
  
  ifelse(y==1,masterDF <- pred.plot.df, masterDF <-rbind(masterDF,pred.plot.df))
}

# Quantiles for the current year myYear 
# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$ps_pred_curr, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

pred.plot.df <- data.frame(dayPSS_all,
                           cumsum(PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= endDayPSS]),
                           t(quant.predPSS))
names(pred.plot.df) <- c("Days","cumcount","low95","low50","median",
                         "up50","up95")

pred.plot.df$Obs <- ifelse(pred.plot.df$Days <= myDay, "Obs", "Realized")

pred.plot.df$Year <- myYear

# Bind into a data.frame
masterDF <- rbind(masterDF, pred.plot.df)

# Plot results
ggplot(masterDF, aes(x = Days,  y = cumcount/1000))+
  geom_point(aes(color = Obs)) +
  # geom_col(alpha = .3)+
  geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.3) +
  geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.3) +
  geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
  coord_cartesian(x = c(148,215))+
  labs(x = "",
       y = "")+
  theme_linedraw(base_line_size = 0)+
  # theme(text = element_text(size = 20, family = "serif"),
  #       panel.grid.major =element_line(),    #strip major gridlines
  #       panel.grid.minor = element_blank(),    #strip minor gridlines
  #       axis.ticks = element_blank(),
  #       panel.background = element_blank(),
  #       axis.line = element_line(color = "black"))+
  
  scale_fill_colorblind(name = "")+
  scale_color_colorblind(name = "")+
  facet_wrap(~Year, nrow = 6)

# PSS only dens plots (Works for all models #################################################################

# Posterior EOS Can-orig run size estimate
df.run <- data.frame("par" = "Runsize",
                     "value" = pars$RunSize)

# PF (prior)
df.prior <- data.frame("par" = "Prior", 
                       "value" = pars$prior_pf)

# PSS prediction density
df.postPredPss <- data.frame("par" = "PSS Prediction", 
                             "value" = pars$post_curr_predPSS)

# Bind into one data frame with 2 columns for plotting
df.comb <- rbind(df.prior,
                 df.postPredPss,
                 df.run)
str(df.comb)
df.comb$par <- as.factor(df.comb$par)

df.comb$par <- relevel(df.comb$par,"Prior")

# Density plots comparing pf, linear prediction, and posterior estimate
ggplot(df.comb, aes(x = value/1000, fill = par))+
  geom_density(alpha = .5)+
  ggtitle(paste("Density plot for",myYear,", day",myDay,"\n Version", model.version))+
  theme_classic()+
  # geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
  #                color = "End of Season Runsize"),
  #            linetype = 2,
  #            size = 1)+
  geom_vline(aes(xintercept  = mean(pars$RunSize)/1000))+
  
  # xlim(c(0,300000))+
  coord_cartesian(xlim = c(0,200))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", 
       fill = "Parameters", 
       y = "Probability Density")+
  # scale_fill_discrete(name = "Parameters", 
  # labels = c( "Preseason Forecast (Prior)", 
  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  scale_fill_colorblind(name = "", 
                        labels = c( "Preseason Forecast (Prior)",
                                    "PSS Prediction",
                                    "Runsize"))+
  # theme_tidybayes(element_text(size = 20))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "top",
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")


ggplot(df.comb, aes(x = value/1000, y =fct_rev( par)))+
  stat_eye(aes(fill = after_stat(level)), .width = c(.8,.95,1), alpha = .5)+
  geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
                 color = "End of Season Runsize"),
             linetype = 2,
             size = 1)+
  coord_cartesian(xlim = c(0,200))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", 
       fill = "Parameters", 
       y = "Parameter")+
  # scale_fill_discrete(name = "Parameters", 
  # labels = c( "Preseason Forecast (Prior)", 
  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  # scale_fill_colorblind(name = "", 
  #                       labels = c( "Preseason Forecast (Prior)",
  #                                   "PSS Prediction",
  #                                   "Runsize"))+
  theme_tidybayes()+
  theme(text = element_text(size = 16, family = "serif"),
        legend.position = "top",
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")+
  guides(fill = guide_legend(override.aes = list(size = 5)))

# PSS & Eagle density plots #################################################
# Posterior EOS Can-orig run size estimate
df.run <- data.frame("par" = "Runsize",
                     "value" = pars$RunSize)

# PF (prior)
df.prior <- data.frame("par" = "Prior", 
                       "value" = pars$prior_pf)

# PSS prediction density
df.postPredPss <- data.frame("par" = "PSS Prediction", 
                             "value" = pars$post_curr_predPSS)

# Eagle prediction density
df.postPredEagle <- data.frame("par" = "Eagle Prediction", 
                               "value" = pars$post_curr_predEagle)

# Bind into one data frame with 2 columns for plotting
df.comb <- rbind(df.prior,
                 df.postPredPss,
                 df.postPredEagle,
                 df.run)
str(df.comb)
df.comb$par <- as.factor(df.comb$par)

df.comb$par <- relevel(df.comb$par,"Prior")

# Density plots comparing pf, linear prediction, and posterior estimate
ggplot(df.comb, aes(x = value/1000, fill = par))+
  geom_density(alpha = .5)+
  ggtitle(paste("Density plot for",myYear,", day",myDay,"\n Version", model.version))+
  theme_classic()+
  geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
                 color = "End of Season Runsize"),
             linetype = 2,
             size = 1)+
  
  # xlim(c(0,300000))+
  coord_cartesian(xlim = c(0,250))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", 
       fill = "Parameters", 
       y = "Probability Density")+
  # scale_fill_discrete(name = "Parameters", 
  # labels = c( "Preseason Forecast (Prior)", 
  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  scale_fill_colorblind(name = "",
                        labels = c( "Preseason Forecast (Prior)",
                                    "Eagle Prediction",
                                    "PSS Prediction",
                                    "Runsize"))+
  # theme_tidybayes(element_text(size = 20))+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")


# Logistic curve fit ###################################################

df.plot <- as.data.frame(cum_PSS_mat_all)
df.plot$day <- dayPSS_all
plot.df <- pivot_longer(df.plot, cols = -day, names_to = "Year", values_to = "count")

plot.df <- plot.df %>% group_by(Year) %>% mutate(obs.prop = count/max(count))

# 
# 
# test <- pars$ps_cumm_pred
# 
# test[,] <- NA
# 
# for (d in 1:n_dayPSS_all) {
#   
#   test[,d] <- pars$ps_cumm_pred[,d]/ pars$scale
# }
# test <- pars
pars$ps_prop_pred
# # Limits
x.lim <- c(145,225)
# y.lim <- c(0,max(apply(pars$ps_pred_curr, 2, quantile, probs=0.975), max(cumsum(PSS_hist$count[PSS_hist$Year == myYear]))))
y.lim <- c(0,325000)

# Plot Observed
plot(x=plot.df$day,
     y=plot.df$obs.prop,
     type='p', pch=21, 
     col='gray',
     xlim=c(148,225),
     # ylim=y.lim,
     xlab='Day',
     ylab='Proportion of EOS Total PSS Passage',
     # main="PSS: Current Year", 
     xaxt='n',
     cex.axis = 1)

axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
axis(side=1, at=seq(from=161, to=201, by=10),
     labels=c("June-10","June-20","June-30","July-10","July-20"),
     cex.axis = 1)
# segments(x0=dayPSS, y0=rep(0, length(dayPSS)), x1=dayPSS, y1=cum_curr_PSS, col='black', lwd=2)
# points(x=dayPSS, y=cum_curr_PSS, type='p', pch=21, bg='gray')
# points(x = c((myDay+1):252), y = cumsum(PSS_hist$count[PSS_hist$Year == myYear])[((myDay+2)-148):(252-147)],
#        type = "p",pch = 23, bg = "green")
# Plot Predicted
polygon(x=c(dayPSS_all,rev(dayPSS_all)), y=c(apply(pars$ps_prop_pred, 2, quantile, probs=0.025),
                                             rev(apply(pars$ps_prop_pred, 2, quantile, probs=0.975))),
        border=FALSE, col=rgb(1,0,0,0.2))

polygon(x=c(dayPSS_all,rev(dayPSS_all)), y=c(apply(pars$ps_prop_pred, 2, quantile, probs=0.25),
                                             rev(apply(pars$ps_prop_pred, 2, quantile, probs=0.75))), border=FALSE, col=rgb(1,0,0,0.2))
lines(x=dayPSS_all, y=apply(pars$ps_prop_pred, 2, median), col=rgb(1,0,0, alpha=0.2), lwd=3)
grid()

