#**********************************************************************************
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Version: 1.0 to 3.4
# Purpose: To make yearly updates to data files including,
#           GSI
#           Total Yearly Reconstructed Canadian Chinook passage
#           Normal curve priors used in PSSnormal_ESprop
#  
#
#
#NOTES:
# 
#   
#    
# 
#
# Packages #########################################################################
require(tidyverse)
require(ggthemes)
require(lubridate)
require(reshape2)
library(bbmle)

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
CAN_hist <- readRDS(file.path(dir.data,"Canadian Chinook Passage RR 2Apr24.RDS"))

# Read in historic preseason forecasts (2013 - 2021)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_12June23_eagle+harvest.RDS"))

# PSS passage
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2023.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 4Apr23.RDS"))


# Update Yearly Reconstructed Canadian Chinook Salmon Run Size  ################################################################
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


# Update normal dist curve priors using the new years PSS Chinook passage #######################

#  Function for fitting normal dist curve to daily PSS passage 
extract.data <- function(df, year) {
  ### TESTING ###
  # df <- PSS_hist
  # year <- 1995
  ###############
  temp.data <- df[df$Year == year,]
  min.day <- min(temp.data$Day)
  max.day <- max(temp.data$Day)
  days <- min.day:max.day
  n.days <- length(days)
  rets <- df$count[df$Year==year]
  cum.rets <- cumsum(rets)
  
  #OUTPUT SECTION
  output <- NULL
  output$min.day <- min.day
  output$max.day <- max.day
  output$days <- days
  output$n.days <- n.days
  output$rets <- rets
  output$cum.rets <- cum.rets
  return(output)
}

# For testing that the function worked....
# out <- extract.data(df = PSS_hist, year = 1995)

# FUNCTION for fitting normal distribution using BBMLE2
like.norm <- function(days, rets, mu, sd, alpha, sigma, plot) {
  ### TESTING ###
  # rets <- out$rets
  # days <- out$days
  # sigma <- -2.3
  # mu <- 5.25
  # sd <- 12.26
  # alpha <- 20.2
  # plot <- TRUE
  # dist <- 'norm'
  ###############
  
  mu <- exp(mu)
  sd <- exp(sd)
  alpha <- exp(alpha)
  sigma <- exp(sigma)
  
  #days <- out$days
  #rets <- out$rets
  n.days <- length(days)
  
  pred <- vector(length=n.days)
  logLike <- vector(length=n.days)
  
  # d <- 5
  # day <- 152
  for(d in 1:n.days) {
    day <- days[d]
    
    # Predicted daily passage
    pred[d] <- alpha*dnorm(day,mu,sd)
    
    # Calculate negative log likelihood for minimization function
    logLike[d] <- -1*(rets[d] - pred[d])^2
    
    
  }#next d
  
  # Plotting the fit to daily PSS arrivals
  if(plot==TRUE) {
    y.lim <- c(0,max(rets, pred, na.rm=TRUE))
    #plot(rets ~ days, type='l', col='gray', xlab='Day', ylab='Catch + Esc', lty=3, ylim=y.lim)
    plot(rets ~ days,
         type='h',
         col='black',
         xlab='Day',
         ylab='PSS Chinook Salmon',
         lty=1,
         ylim=y.lim,
         xlim = c(148,220),
         lwd=5,
         cex.axis = 1.5,
         cex.lab = 1.5,
         xaxt='n')
    #points(x=days, y=rets, pch=21, bg='gray')
    #Model data
    lines(x=days, y=pred, lwd=2, col='red')
    # axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
    axis(side=1, at=seq(from=152, to=202, by=10),
         labels=c("June-1","June-11","June-21","July-1","July-11","July-21"),
         cex = 1.5)
  }
  
  #Negative Log Likelihood
  NLL <- -1*sum(logLike)
  return(NLL)
}
# # End of functions

# Fitting to each year in the PSS dataframe
newYear <- 2023

years <- unique(PSS_hist$Year[PSS_hist$Year <= newYear])
n.years <- length(years)

norm.mu <- vector(length=n.years)
norm.sd <- vector(length=n.years)
norm.alpha <- vector(length=n.years)
# pdf(file = file.path(dir.output,"4.0 Normal MLE Fit 15Jul22"), height=8, width=8)
# par(mfrow = c(3,3))
# y <- 1
for(y in 1:n.years) {
  
  # y = 2
  
  year <- years[y]
  print(year)
  # Extract data NOTE! Change GSI if needed for model
  out <- extract.data(df = PSS_hist, year=year)
  #Fit normal data
  fit.norm <-  mle2(like.norm,
                    start=list(mu=log(170), sd=log(20), alpha=log(6000)),
                    fixed=list(plot=FALSE, sigma=log(20)),
                    data=list(days=out$days, rets=out$rets),
                    method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
  
  # Extract param values
  norm.mu[y] <- exp(coef(fit.norm)[1])
  norm.sd[y] <- exp(coef(fit.norm)[2])
  norm.alpha[y] <- exp(coef(fit.norm)[3])
  
  # Plot the fit
  like.norm(days=out$days,
            rets=out$rets,
            coef(fit.norm)[1],
            coef(fit.norm)[2],
            coef(fit.norm)[3],
            coef(fit.norm)[4],
            plot=TRUE)
  mtext(year, side=3, line=2, font=2)
}#next y

# dev.off()
# Data frame of parameter estimates
prior.df <- data.frame(norm.mu,norm.sd,norm.alpha,years)

# Name the colums
names(prior.df) <- c('mid','sd','alpha',"year")

# Save as RDS
# saveRDS(object = prior.df, file = file.path(dir.data,"normal curve parameters All Chinook 1995_2023.RDS"))

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

