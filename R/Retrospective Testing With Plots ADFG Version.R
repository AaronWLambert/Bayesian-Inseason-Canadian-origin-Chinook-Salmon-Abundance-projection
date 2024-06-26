#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.27.24
# Version: ADFG Retrospective testing
# Purpose: Compare model perfomance using leave-one-out retro testing.
# 
#
#   1) Read in data
#   2) Call the function
#   3) Generate retrospective outputs such as figures and values.
#
#
#=================================================================================
# NOTES:
# 
#
# 
# Next steps: 

# 
#
#=================================================================================
library(rstan)
library(bayesplot)
library(tidyverse)
library(ggthemes)
library(viridis)
library(shinystan)
library(lubridate)
library(ggpubr)
library(gridExtra)
library(tidybayes)
library(wesanderson)
library(grid)
library(bbmle)

# Parralize for optimum model run time
rstan_options(auto_write = TRUE)
#
mc.cores = parallel::detectCores()
# mc.cores <-1

# Define Workflow Paths =========================================================================

# Change to the working directory on your computer
wd <- "C:/Users/aaron/OneDrive/Desktop/ADFG Yukon Model/Yukon Chinook Bayesian Inseason Projection Model ADFG/"

setwd(wd)

# Folders that must be named in the working directory
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Functions #################################################################################

# Functions to run model and get day of year using the date
source(file = file.path(dir.R,"model_run_function LOO_ADFG.R"))

# Function to get retrospective stats for MAPE, RMSE, and PE
source(file = file.path(dir.R,"Retro Function ADFG.R"))

######### Import Data #######################################################################
# Note that these files must be updated to reflect the current methods

# JTC preseason forecasts 
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_4Apr24_needupdatingtonewmethod.RDS"))

# EOS Can-orig reconstructed counts (NEW METHOD; Connor et al)
CAN_hist <- readRDS(file.path(dir.data,"Canadian Chinook Passage RR 2Apr24.RDS"))

# PSS Daily Passage
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2023.RDS"))

# Eagle Sonar daily passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage Final 2023.RDS"))

# Mean GSI by strata (naive estimator)
# GSI_mean <- readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# GSI data by year and strata  (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 4Apr24.RDS"))

# Read in PSS observation error estimates
PSS_sd <- readRDS(file = file.path(dir.data,"PSS SD 1995_2021.RDS"))

# Shape parameters for logistic arrival distribution
logistic.all <- read.csv(file = file.path(dir.data,
                                          "logistic curve parameters All Chinook 1995_2022.csv"))

# Shape parameters for normal arrival distribution
# normal.all <- read.csv(file = file.path(dir.data,
#                                         "normal curve parameters All Chinook 1995_2022.csv"))
normal.all <- readRDS(file = file.path(dir.data,"normal curve parameters All Chinook 1995_2023.RDS"))



# Control Section #################################################################################

# Year to conduct retro testing through
end.year = 2022

# Input model for retro testing
model.version <- "PSSnormal_ESprop"

# MCMC Parameters
n.chains <- 4   # Number of chains to run
n.iter <- 1000  # Number of iterations per chain
n.thin <- 2     # Thinning rate 

# Days to use in retrospective testing runs ##############
# Test days used in full season run (every 5 days starting June 2)
testDays <- seq(from = 153, to = 243, by = 5)

# Years included in retrospective testing
testYears <- c(2007:end.year)

# Retrospective testing loop ######################

# # List to store outputs
outputList<-list()

# Increse memory limit (May need this for PSSnormal and PSSlogistic)
# memory.limit(size = 60000)

# Display warnings while looping over days and years
options(warn = 1)
for(y in c(testYears)){
  for(d in c(testDays)){
    
    outputList[[paste("",y,"_",d, sep = "")]]<-InSeasonProjection(model.version = model.version,
                                                                  myYear = y,
                                                                  myDay = d,
                                                                  n.chains = n.chains,
                                                                  CAN_hist = CAN_hist,
                                                                  pf_hist = pf_hist,
                                                                  PSS_hist = PSS_hist,
                                                                  PSS_sd = PSS_sd,
                                                                  n.thin = n.thin,
                                                                  n.iter = n.iter,
                                                                  GSI_by_year = GSI_by_year,
                                                                  Eagle_hist = Eagle_hist,
                                                                  normal = TRUE,
                                                                  logistic = FALSE,
                                                                  prior.df.log = logistic.all,
                                                                  prior.df.norm = normal.all,
                                                                  startDayPSS = 148,
                                                                  startYearPSS = 2005,
                                                                  endYear = end.year
                                                                  )
    print(paste("Day =",d,"Year =",y))
  } #dloop
  
  print(paste("Finally done with year",y))
} #yloop


# Save or read in outputlist #######################################

# Save the resulting model output here. this will save to the output folder.
saveRDS(object = outputList, file = file.path(dir.output,
                                              "VerPSSreg 27Mar24.RDS"))

# Read in outputs from retro testing 

# **Note** Increase memory limit if using PSSnorm or PSSlogistic
# memory.limit(size = 60000)

# Read in the model output
# outputlist_verPSSreg <- readRDS(file = file.path(dir.output,"VerPSSreg 27Mar24.RDS")) 

##### Calculations for retrospective testing #############################################

# Uses retrospective.function
#  This function calculates RMSE, MAPE, PE, and precision (not currently used)

# PSSreg 
RetroList_verPSSreg <- retrospective.function(outputList = outputlist_verPSSreg, # The model results from above
                                         testYears = testYears,                  # The years tested 
                                         testDays = testDays,                    # The days across the season tested
                                         CAN_hist = CAN_hist,                    # Reconstructed run size
                                         pf_hist = pf_hist,                      # PF (only needed for PF retro)
                                         startYearRetro = 2007,                  # First year of retro testing (default 2007)
                                         endYearRetro = 2022,                    # End year retro testing (default 2022)
                                         pf = FALSE)                             # Is this the PF? (no)

# Preseason forecast
RetroList_pf <- retrospective.function(outputList = outputlist_verPSSreg, # Use any output for the PF retro
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_new,
                                           pf_hist = pf_hist,
                                           startYearRetro = 2007,
                                           endYearRetro = 2022,
                                           pf = TRUE)


# RMSE ######################################################################################

# Extract RMSE into dataframe for each model tested

# PSSreg
rmseDF_verPSSreg <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_verPSSreg$RMSE_by_day_vect,
                              "Version" = "PSSreg")


# Preseason forecast
rmseDF_pf <- data.frame("Day" = testDays,
                            "RMSE" = RetroList_pf$RMSE_by_day_vect,
                            "Version" = "PF")

# rmseDF_pf_old <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_pf_old$RMSE_by_day_vect,
#                             "Version" = "PF Eagle+Harv")

# Combine data frames into one data.frame for plotting
full_rmseDF <- rbind(
  rmseDF_verPSSreg,
  rmseDF_pf_new # Uncomment if you want this in the excel-like table below
)

# Make version a factor
full_rmseDF$Version <- as.factor(full_rmseDF$Version)

# Day names for better looking plots
day.names <- c(
  "June 2",
  "June 7",
  "June 12",
  "June 17",
  "June 22",
  "June 27",
  "July 2",
  "July 7",
  "July 12",
  "July 17",
  "July 22",
  "July 27",
  "Aug 1", #,
  "Aug 6",
  "Aug 11",
  "Aug 16",
  "Aug 21",
  "Aug 26",
  "Aug 31"
)

# Plot to look at RMSE as a color coded excel-like table
ggplot(full_rmseDF, aes(x = factor(Day), y = factor(Version), fill = RMSE))+
  geom_tile()+
  geom_text(aes(label = round(RMSE)))+
  scale_fill_distiller(palette = "Spectral")+
  labs(x = "Day of Year",
       y = "Model Versions")+
  scale_x_discrete(breaks = c(testDays), labels = day.names)+
  theme(text = element_text(size = 14)
  )

# Add PF for plotting
full_rmseDF$PF <- rmseDF_pf_new$RMSE

# Save plot as png
# png(filename = file.path(dir.figs,"Chap 1RMSE long Plot 28Feb24.png"), width = 1000, height = 800)

# RMSE plot
ggplot(full_rmseDF, aes(x = Day, y = RMSE/1000
                        # col = Version,
                        # shape= Version, 
                        # alpha = 0.7
                        ))+
  geom_point()+
  geom_line()+
  labs(y = "RMSE (1000's)")+
  # geom_col(data = short_rmseDF, position = "dodge", width = 3)+
  coord_cartesian(ylim = c(
    min(full_rmseDF$RMSE/1000)-.5,
    max(full_rmseDF$RMSE/1000)+.5
    # 10000,15000
    
  ))+
  guides(alpha = "none")+
  scale_x_continuous(breaks = c(testDays), labels = day.names)+
  geom_hline(aes(yintercept = PF/1000, # PF for a reference
                 color = factor(PF)),
             size = 1,linetype = 2,
             show.legend = T)+
  facet_wrap(~Version, ncol = 2)+ # facet if compairing at more than a few models at a time
  scale_color_discrete(labels = "Preseason Forecast",
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 0.5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))

# Uncomment to turn off the save function
# dev.off()

# MAPE ####################################################################

# Create MAPE data frame
mapeDF_verPSSreg <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_verPSSreg$MAPE_vect,
                          "Version" = "Ver PSSreg")

mapeDF_PF_new <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_new$MAPE_vect,
                            "Version" = "PF")

# Combine into one data-frame
full_MAPE.df <- rbind(
  mapeDF_verPSSreg
  # mapeDF_PF_new # Uncomment for use in excel-like table
)



# MAPE Table for papers
# MAPE_table<-as.data.frame(pivot_wider(full_MAPE.df,names_from = Version, values_from = MAPE))
# MAPE_table_round <- format(MAPE_table, digits = 3)

# write.table(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 15Jul22"),
#             row.names = F)

# Add PF for plotting
full_MAPE.df$PF <- mapeDF_PF_new$MAPE

ggplot(full_MAPE.df, aes(x = factor(Day), y = factor(Version), fill = MAPE*100))+
  geom_tile()+
  geom_text(aes(label = round(MAPE*100,1)))+
  scale_fill_distiller(palette = "Spectral", name = "MAPE")+
  scale_x_discrete(breaks = c(testDays), labels = day.names)+
  labs(x = "Day of Year",
       y = "Model Versions")+
  theme(text = element_text(size = 14))

# Save the MAPE plot
# png(file = file.path(dir.figs,"MAPE 27Mar24.png"),width = 1000, height = 800)

# MAPE plot
ggplot(full_MAPE.df, aes(x = Day,
                         y = MAPE*100, 
                         # col = Version,
                         # shape = Version
))+
  # geom_col(position = "dodge",width = 3) + 
  
  geom_point(size = 1.5)+
  geom_line(aes(group = Version))+
  geom_hline(aes(yintercept = PF*100, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  coord_cartesian(ylim = c(min(full_MAPE.df$MAPE*100)-2,max(full_MAPE.df$MAPE*100)+2))+
  ylab('MAPE (%)')+
  scale_x_continuous(breaks = c(testDays), labels = day.names)+
  facet_wrap(~Version, ncol = 2)+
  scale_color_discrete(labels = "Preseason Forecast",
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))

# dev.off()


# Percent Error Box Plots ###########################################################

PE_verPSSreg <- as.data.frame(RetroList_verPSSreg$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_verPSSreg$version <- "Ver PSSreg"

PE_pf <- as.data.frame(RetroList_pf$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf$version <- "PF"

peDF_total <- rbind(
  PE_verPSSreg,
  PE_pf
)

# For plotting PF point shape 
newDF <- left_join(x = pf_hist, y = CAN_hist)

newDF$Year <- as.factor(newDF$Year)

newDF$PF_old <- (newDF$mean-newDF$can.mean)/newDF$can.mean

# Uncomment to save PE plot as png
# png(file = file.path(dir.figs,"Chap 1 PE Plot across years yaxis Cropped 28Feb24.png"),
# width = 1000, height = 800)

# PE plot for each day across years
# ggplot(peDF_total[!peDF_total$Year %in% c(2021,2022),], 
ggplot(peDF_total, 
       aes(x = version,
           y = PE*100, 
           fill =version))+
  geom_boxplot(
    # outlier.color = NA
    # position = position_dodge(width = 1.01)
    
  )+ 
  # geom_point()+
  # scale_fill_colorblind()+
  scale_fill_manual(values = cust.col.pe)+
  # scale_color_colorblind(name = "")+
  geom_hline(yintercept  = 0, linetype = 2, color = "orange", size = 1)+
  guides(alpha = "none")+
  labs(fill = "Model Version",
       y = "Percent Error",
       x = "",
       color = "")+
  coord_cartesian(ylim = c(-40,60))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text.x =   element_blank(),
        axis.ticks.x = element_blank())+
  facet_wrap(~Day, labeller = as_labeller(label_names)
  )

# dev.off()
# ggtitle("New Can_hist ver 2")


# Turn off PDF save
# dev.off()




