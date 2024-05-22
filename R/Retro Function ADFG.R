#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 05.11.24
# Version: ADFG handoff
# Purpose: Functions that will calculate retrospective testing metrics
# 
#
#
#
#
#=================================================================================
#NOTES:
# 
# 
# Next steps: 
#
#=================================================================================

# Function for MAPE, PE, RMSE, and other things listed in the function description

#' Title
#'
#' @param outputList Model output list from retrospective testing
#' @param testYears Years included in the retrospective testing
#' @param testDays Days included in retrospective testing
#' @param CAN_hist A data.frame of Yukon River Canadian-origin Chinook salmon reconstructed counts
#' @param pf_hist A data.frame of preseason forecast point estimates for years 2007-2022
#' @param param The parameter of interest. Default is RunSize, the posterior integrated run size estimate
#' @param startYearRetro Year retrospective testing begins. Default 2007
#' @param endYearRetro Year retrospective testing ends. Default is 2022
#' @param pf True or False that the retrospective testing is for the preseason forecast.
#'
#' @return Returns a list of consisting of: a matrix of the median of the parameter of interest. 
#' A matrix of observed run sizes. A matrix of percent error. A matrix of absolute perecent error.
#' A vector of MAPE calculated across years for each day tested.
#' A vector of RMSE calculated across years for each day tested.
#' 3 Matrices indicating if the true runsize was observed in the 95, 80, and 50 CI.
#' A vector of testdays used.
#' A vector of years included in retrospective testing.
#' @export
#'
#' @examples
retrospective.function<- function(outputList, 
                                  testYears, 
                                  testDays, 
                                  CAN_hist,
                                  pf_hist,
                                  param = "RunSize",
                                  startYearRetro = 2007,
                                  endYearRetro = 2022,
                                  pf = FALSE){
  # OutputList <- a list of outputs from the Stan model runs for days and years
  # testYears <- vector of days included in OutputList
  # testDays <- vector of years included in OutputList
  # CAN_hist <- dataframe of historic EOS reconstructed run sizes
  # Returns a list including:
  #     Matrix of medians from posterior RunSize projections
  #     Matrix of historic EOS Reconstructed Can-orig Runsize
  #     Percent error matrix
  #     APE matrix
  #     MAPE for days across years vector
  #     RMSE for days across years vector
  #     RMSE for years across days vector
  #     Matrices for how often the "true" EOS can-orig runsize falls in the 
  #      50, 80, and 95% intervals.
  #
  # Calculations for retrospecitve testing #############################################
  if(pf == FALSE){
    # Create matrix of observed and median outputs
    # 
    # Create matrix with dimensions of days by years that are in the outputList from model run
    median_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(median_mat) <- c(testYears)
    rownames(median_mat) <- c(testDays)
    
    # Make sure it works...
    dim(median_mat)
    str(median_mat)
    
    # Extract medians from posterior RunSize prediction and place in vector for matrix pop.
    # populate matrix with medians from posterior RunSize prediction
    # counter <- 1
    # for (y in 1:length(testYears)) {
    #   for (d in 1:length(testDays)) {
    #     
    #     median_mat[d,y]<- median(outputList[[counter]]$pars$RunSize)
    #     
    #     counter <- counter+1
    #   }
    # }
    
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        y1<-testYears[y]
        d1<-testDays[d]
        
        median_mat[d,y]<- median(outputList[[paste(y1,"_",d1,sep = "")]]$pars[[param]])
        
        # print(d)
        # print(y)
      }
    }
    median_mat
    
    # Matrix with realized reconstructed counts by year and day 
    
    realized_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(realized_mat) <- c(testYears)
    rownames(realized_mat) <- c(testDays)
    
    # Make sure it works...
    dim(realized_mat)
    
    # Extract values for realized runsize
    realized_vect <- CAN_hist$can.mean[CAN_hist$Year >= startYearRetro] 
    
    counter <- 1
    # Populate the realized matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        realized_mat[d,y]<- realized_vect[counter]
        
      }
      counter <- counter+1
    }
    realized_mat
    # Percent Error (PE) matrix
    PE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(PE_mat) <- c(testYears)
    rownames(PE_mat) <- c(testDays)
    
    # Populate the PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        PE_mat[d,y]<- (median_mat[d,y]- realized_mat[d,y])/realized_mat[d,y]
        
      }
    }
    # Percent Error (APE) matrix
    APE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(APE_mat) <- c(testYears)
    rownames(APE_mat) <- c(testDays)
    
    # Populate the absolute PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        APE_mat[d,y]<- abs(realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
        
      }
      
    }
    
    
    # Calculate MAPE by Days
    MAPE_day_vect <- vector(length = length(PE_mat[,1]))
    names(MAPE_day_vect) <- testDays
    
    for (d in 1:length(MAPE_day_vect)) {
      MAPE_day_vect[d] <- mean(APE_mat[d,])
    }
    
    # Root Mean Square error for days across years
    # Square Error (SE) matrix
    SE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(SE_mat) <- c(testYears)
    rownames(SE_mat) <- c(testDays)
    
    # Populate the SE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        SE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])^2
        
      }
      
    }
    
    # Calculate RMSE by days
    RMSE_day_vect <- vector(length = length(SE_mat[,1]))
    names(RMSE_day_vect) <- testDays
    
    for (d in 1:length(RMSE_day_vect)) {
      RMSE_day_vect[d] <- sqrt(mean(SE_mat[d,]))
    }
    
    # calculate RMSE by year
    RMSE_year_vect <- vector(length = length(SE_mat[1,]))
    names(RMSE_year_vect) <- as.character(testYears)
    
    for (y in 1:length(RMSE_year_vect)) {
      RMSE_year_vect[y] <- sqrt(mean(SE_mat[,y]))
    }
    
    # How often does the interval capture the true return size? ######
    # 95% Interval
    true95Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true80Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true50Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    # Give names to matrix
    colnames(true95Mat) <- c(testYears)
    rownames(true95Mat) <- c(testDays)
    
    colnames(true80Mat) <- c(testYears)
    rownames(true80Mat) <- c(testDays)
    
    colnames(true50Mat) <- c(testYears)
    rownames(true50Mat) <- c(testDays)
    
    count <- 1
    
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        # 95%
        test95<-quantile(outputList[[count]]$pars$RunSize, probs = c(.025,.975))
        
        true95Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test95[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test95[1],
                                 1,
                                 0)
        # 80%
        test80<-quantile(outputList[[count]]$pars$RunSize, probs = c(.1,.9))
        
        true80Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test80[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test80[1],
                                 1,
                                 0)
        # 50%
        test50<-quantile(outputList[[count]]$pars$RunSize, probs = c(.25,.75))
        
        true50Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test50[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test50[1],
                                 1,
                                 0)
        count <- count+1
      }
    }
    
  }else{
    # Create matrix with dimensions of days by years that are in the outputList from model run
    median_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(median_mat) <- c(testYears)
    rownames(median_mat) <- c(testDays)
    
    # Make sure it works...
    dim(median_mat)
    str(median_mat)
    # Extract values for realized runsize
    pf_vect <- pf_hist$mean[pf_hist$Year>=startYearRetro &
                              pf_hist$Year <= endYearRetro]
    
    counter <- 1
    # Populate the realized matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        median_mat[d,y]<- pf_vect[counter]
        
      }
      counter <- counter+1
    }
    # Extract medians from posterior RunSize prediction and place in vector for matrix pop.
    # populate matrix with medians from posterior RunSize prediction
    # counter <- 1
    # for (y in 1:length(testYears)) {
    #   for (d in 1:length(testDays)) {
    #     
    #     median_mat[d,y]<- median(outputList[[counter]]$pars$prior_pf )
    #     
    #     counter <- counter+1
    #   }
    # }
    # median_mat
    
    # Matrix with realized reconstructed counts by year and day 
    
    realized_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(realized_mat) <- c(testYears)
    rownames(realized_mat) <- c(testDays)
    
    # Make sure it works...
    dim(realized_mat)
    
    # Extract values for realized runsize
    realized_vect <- CAN_hist$can.mean[CAN_hist$Year >= startYearRetro &
                                         CAN_hist$Year <= endYearRetro ]
    
    counter <- 1
    # Populate the realized matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        realized_mat[d,y]<- realized_vect[counter]
        
      }
      counter <- counter+1
    }
    realized_mat
    # Percent Error (PE) matrix
    PE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(PE_mat) <- c(testYears)
    rownames(PE_mat) <- c(testDays)
    
    # Populate the PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        PE_mat[d,y]<- (median_mat[d,y]- realized_mat[d,y])/realized_mat[d,y]
        
      }
    }
    # Percent Error (APE) matrix
    APE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(APE_mat) <- c(testYears)
    rownames(APE_mat) <- c(testDays)
    
    # Populate the absolute PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        APE_mat[d,y]<- abs(realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
        
      }
      
    }
    
    
    # Calculate MAPE by Days
    MAPE_day_vect <- vector(length = length(PE_mat[,1]))
    names(MAPE_day_vect) <- testDays
    
    for (d in 1:length(MAPE_day_vect)) {
      MAPE_day_vect[d] <- mean(APE_mat[d,])
    }
    
    # Root Mean Square error for days across years
    # Square Error (SE) matrix
    SE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(SE_mat) <- c(testYears)
    rownames(SE_mat) <- c(testDays)
    
    # Populate the SE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        SE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])^2
        
      }
      
    }
    
    # Calculate RMSE by days
    RMSE_day_vect <- vector(length = length(SE_mat[,1]))
    names(RMSE_day_vect) <- testDays
    
    for (d in 1:length(RMSE_day_vect)) {
      RMSE_day_vect[d] <- sqrt(mean(SE_mat[d,]))
    }
    
    # calculate RMSE by year
    RMSE_year_vect <- vector(length = length(SE_mat[1,]))
    names(RMSE_year_vect) <- as.character(testYears)
    
    for (y in 1:length(RMSE_year_vect)) {
      RMSE_year_vect[y] <- sqrt(mean(SE_mat[,y]))
    }
    
    # How often does the interval capture the true return size? ######
    # 95% Interval
    true95Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true80Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true50Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    # Give names to matrix
    colnames(true95Mat) <- c(testYears)
    rownames(true95Mat) <- c(testDays)
    
    colnames(true80Mat) <- c(testYears)
    rownames(true80Mat) <- c(testDays)
    
    colnames(true50Mat) <- c(testYears)
    rownames(true50Mat) <- c(testDays)
    
    count <- 1
    
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        # 95%
        test95<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.025,.975))
        
        true95Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test95[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test95[1],
                                 1,
                                 0)
        # 80%
        test80<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.1,.9))
        
        true80Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test80[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test80[1],
                                 1,
                                 0)
        # 50%
        test50<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.25,.75))
        
        true50Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test50[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test50[1],
                                 1,
                                 0)
        count <- count+1
      }
    }
  } 
  
  RetroList <- list("Med_mat" = median_mat, 
                    "Realized_mat" = realized_mat,
                    "PE_mat" = PE_mat,
                    "APE_mat" = APE_mat,
                    "MAPE_vect" = MAPE_day_vect,
                    "RMSE_by_day_vect" = RMSE_day_vect,
                    "RMSE_by_year_vect" = RMSE_year_vect,
                    "True95" = true95Mat,
                    "True80" = true80Mat,
                    "True50" = true50Mat,
                    "Version" = outputList$`2013_153`$version,
                    "TestDaysIncluded" <- testDays,
                    "TestYearsIncluded" <- testYears)
  
  return(RetroList)
}

# Retro Plots function ###########################################################3

# Functions for extracting parameter values ################################################

runsize_trend_func <- function(mod, testYears, testDays, CAN_hist, PF){
  
  # Matrix of Runsize estimates
  RunSize_vect <- vector(length=length(mod))
  
  for (p in 1:length(RunSize_vect)) {
    RunSize_vect[p]<- median(mod[[p]]$pars$RunSize)
  }
  
  # RunSize_mat <-matrix(RunSize_vect, 
  #                      nrow = length(testDays), 
  #                      ncol = length(testYears),
  #                      byrow = FALSE) 
  RunSize_mat <-matrix(RunSize_vect, 
                       nrow = length(testDays), 
                       ncol = length(testYears),
                       byrow = FALSE)
  colnames(RunSize_mat)<- c(testYears)
  RunSize_DF <- as.data.frame((RunSize_mat))
  # RunSize_DF$day <- testDays
  RunSize_DF$day <- testDays
  RunSize_DF <- as.data.frame(pivot_longer(RunSize_DF,cols = -day))
  names(RunSize_DF)<- c("Day", "Year", "RunSize")
  RunSize_DF$Year<- as.double(RunSize_DF$Year)
  RunSize_DF <- left_join(RunSize_DF,CAN_hist)
  RunSize_DF <- left_join(RunSize_DF, PF)
  
  # Matrix of 2007 Runsize estimates
  PSSpred_vect <- vector(length=length(mod))
  
  for (p in 1:length(RunSize_vect)) {
    PSSpred_vect[p]<- median(mod[[p]]$pars$post_curr_predPSS)
  }
  
  # Matrix of PSS prediction
  PSSpred_mat <-matrix(PSSpred_vect,
                       nrow = length(testDays),
                       ncol = length(testYears),
                       byrow = FALSE)
  # PSSpred_mat <-matrix(PSSpred_vect, 
  #                      nrow = length(testDays_short), 
  #                      ncol = length(testYears),
  #                      byrow = FALSE) 
  colnames(PSSpred_mat)<- c(testYears)
  PSSpred_DF <- as.data.frame((PSSpred_mat))
  # PSSpred_DF$day <- testDays
  PSSpred_DF$day <- testDays
  PSSpred_DF <- as.data.frame(pivot_longer(PSSpred_DF,cols = -day))
  names(PSSpred_DF)<- c("Day", "Year", "PSSpred")
  PSSpred_DF$Year<- as.double(PSSpred_DF$Year)
  RunSize_DF <- left_join(PSSpred_DF,RunSize_DF)
  RunSize_DF$Version <- as.factor(mod[[1]]$version)
  
  return(RunSize_DF)
} # End of function