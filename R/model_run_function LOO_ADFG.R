#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.26.24
# Version: LOO-ADFG *This version is cleaned up for ADFG*

# Purpose: Function for looping through days and years using leave-one-out 
#          retrospective testing. 
# 
#
#
#
#
#=================================================================================
#NOTES:
# This function will be called from another script.
# 
# Next steps: 
#
#=================================================================================



InSeasonProjection <- function(model.version,      # Model version to run. (character string)
                               myYear,             # Year the projection is made
                               myDay ,             # Day the projection is made
                               n.chains,           # Number of chains to run in Stan
                               n.iter,             # Number of iterations
                               n.thin ,            # Thinning rate
                               CAN_hist,           # EOS total Canadian-origin Chinook est (DF)
                               PSS_hist,           # Daily PSS passage (DF)
                               GSI_by_year,        # Genetic apportionment of Can-origing Chinook
                               pf_hist,            # PF for 2007-present
                               Eagle_hist,         # Eagle sonar passage
                               PSS_sd,             # Variance estimates for PSS
                               norton.sst,         # SST
                               Emmonak_Air_Temp,   # Air Temp
                               prior.df.log,       # priors for logistic params
                               prior.df.norm,      # priors for normal params
                               savefit = FALSE,    # Whether to save the fit object
                               logistic = FALSE,   # Running version 5?
                               normal = FALSE,     # Running version 4?
                               multiplier = 1,     # Remove me!!!!!
                               startDayPSS = 148,  # First day of PSS calculation
                               startYearPSS = 1995,# First year to include in PSS
                               endYear             # Last Complete Year
){
  
  ## For testing
  # model.version = model.version
  # myYear = 2011
  # myDay = 183
  # n.chains = 4
  # n.iter = 10000
  # n.thin = n.thin
  # CAN_hist = CAN_hist_old
  # PSS_hist = PSS_hist
  # Eagle_hist = Eagle_hist
  # # GSI_mean = GSI_mean,
  # pf_hist = pf_hist
  # GSI_by_year = GSI_by_year
  # normal = FALSE
  # logistic = FALSE
  # startDayPSS = 148
  # startYearPSS = 1995
  # prior.df.log = logistic.all
  # prior.df.norm = normal.all
  # PSS_sd = PSS_sd
  # 
  
  
  
  # Start Years for Predictors
  startYearPF <- 2007
  startYearGSI <- 2005
  # startYearPSS <-1997
  
  # Wont typicaly change;
  # day 152 = June 1
  # startDayPSS <-148
  endDayPSS <- 250
  
  # Preseason Forcast######################
  
  # Current Preseason forecast CAN-orig Chinook
  pf <- log(pf_hist$mean[pf_hist$Year == myYear])
  
  # Centering PF draws around true run size
  # pf_milo <- log(CAN_hist$can.mean[CAN_hist == myYear])
  
  # Vector of historical preseason forecasts for to compute sd for prior in stan model
  PF_vect <- pf_hist$mean[
    pf_hist$Year <= endYear &
    pf_hist$Year != myYear]
  
  names(PF_vect) <-  pf_hist$Year[
    pf_hist$Year <= endYear &
    pf_hist$Year != myYear]
  
  
  # EOS reconstructed runsize for historic years that have a preseason forecast 
  EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 & 
                             CAN_hist$Year != myYear]
  
  names(EOS) <- CAN_hist$Year[CAN_hist$Year >= 2007 & 
                                CAN_hist$Year != myYear
  ]
  
  # SD for preseason forecast prior
  pf_sigma <- sd(log(PF_vect/EOS))
  
  # Years in PF
  yearPF <- unique(pf_hist$Year[pf_hist$Year != myYear 
                                &  pf_hist$Year <= endYear
  ])
  
  n_yearsPF <- length(yearPF)
  
  
  # Historic EOS Reconcstructed Can-origin Run Size ############################################
  
  # Vector of run sizes excluding the year of interest
  totalEOS <- CAN_hist$can.mean[CAN_hist$Year != 1996 &   # No PSS 1996
                                  CAN_hist$Year != myYear &
                                  CAN_hist$Year >= startYearPSS ]
  
  # Name the elements for accounting purposes
  names(totalEOS) <- CAN_hist$Year[CAN_hist$Year != 1996 &   # No PSS 1996
                                     CAN_hist$Year != myYear &
                                     CAN_hist$Year >= startYearPSS ]
  
  # Number of years included in EOS
  n_totalEOS <- length(totalEOS)
  
  # Pilot Station Sonar Data ##############################################################
  
  # Vector of historic PSS years excluding myYear
  yearPSS <- unique(PSS_hist$Year[PSS_hist$Year != myYear &
                                    PSS_hist$Year != 1996 &
                                    PSS_hist$Year >= startYearPSS
  ])
  
  # Number of years used in model
  n_yearPSS <- length(yearPSS)
  
  # # PSS days included up to myDay
  # #  I.e., June 1 = 1, June 2 = 2 ....
  # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])-147 *******
  # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & PSS_hist$Day>= startDayPSS])-147
  dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & 
                                  PSS_hist$Day>= startDayPSS])
  
  # Number of days used
  n_dayPSS <- length(dayPSS)
  
  # Number of total days in the season
  dayPSS_all <- startDayPSS:endDayPSS
  
  n_dayPSS_all <- length(dayPSS_all)
  
  loc.allDays.myDay <- which(dayPSS_all %in% myDay)
  
  # PSS daily passage estimate for days up to myDay for the year of interest (myYear)
  curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear &
                                PSS_hist$Day <= myDay &
                                PSS_hist$Day>=startDayPSS])
  
  names(curr_PSS) <- (PSS_hist$Day[PSS_hist$Year == myYear &
                                     PSS_hist$Day <= myDay &
                                     PSS_hist$Day>=startDayPSS])
  
  # Number of days used
  n_curr_PSS <- length(curr_PSS)
  
  # Cumulative current PSS counts
  cum_curr_PSS <- cumsum(curr_PSS)
  
  # Create matrix with dimensions days x historic Year
  PSS_mat <- matrix(nrow = n_dayPSS,
                    #start 1996 to account for missing year 1996
                    ncol = n_yearPSS) 
  
  # Give names to matrix for accounting
  colnames(PSS_mat) <- yearPSS
  
  rownames(PSS_mat) <- dayPSS
  
  # Create vector of counts from relevant years and days to input into matrix
  
  (count_vect <- PSS_hist$count[PSS_hist$Year != myYear &
                                  PSS_hist$Year <= max(CAN_hist$Year) &
                                  PSS_hist$Day <= (myDay) &
                                  PSS_hist$Day >= (startDayPSS) &
                                  PSS_hist$Year >= startYearPSS])
  
  # The number of observations in count_vector for PSS historic counts
  # n_hist_counts <- length(count_vect)
  
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
  
  cumPSS <- apply(PSS_mat,2,sum)
  
  # PSS martix of all the passage for historic/retro years
  
  # Create matrix with dimensions myday and myYear
  PSS_mat_all <- matrix(nrow = n_dayPSS_all,
                        #start 1996 to account for missing year 1996
                        ncol = n_yearPSS) 
  
  # Give names to matrix for accounting
  colnames(PSS_mat_all) <- yearPSS
  
  rownames(PSS_mat_all) <- dayPSS_all
  
  # Create vector of counts from relevant years and days to input into matrix
  (count_vect_all <- PSS_hist$count[PSS_hist$Year != myYear &
                                      PSS_hist$Year >= startYearPSS &
                                      PSS_hist$Day <= endDayPSS &
                                      PSS_hist$Day >= startDayPSS ])
  
  
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
  
  cumPSS_all <- apply(PSS_mat_all, 2, sum)
  
  # Copy matrix
  PSS_mat_prop_all <- PSS_mat_all
  
  cum_PSS_mat_all <- apply(PSS_mat_all,2,cumsum)
  
  # Populate with NA's
  PSS_mat_prop_all[,] <- NA
  
  # Fill with observed proportions by day and year
  for (y in 1:n_yearPSS){
    for (d in 1:n_dayPSS_all) {
      
      
      PSS_mat_prop_all[d,y] <- cum_PSS_mat_all[d,y]/ cum_PSS_mat_all[n_dayPSS_all,y]
      
      
    }
  }
  # #### Double first entry that is non-zero to make up for non-true zero days #############
  # for(y in 1:n_yearPSS){
  # 
  #   PSS_mat[,y][PSS_mat[,y]>0][1] <- PSS_mat[,y][PSS_mat[,y]>0][1]*2
  # }
  ###########################################################################
  
  
  # Eagle Sonar preprocessing #####################################################################
  
  # Vector of historic Eagle Sonar passage excluding myYear
  yearEagle <- unique(Eagle_hist$Year[Eagle_hist$Year != myYear 
                                      # & Eagle_hist$Year <= 2021
  ])
  
  # Number of Eagle sonar years used in model
  n_yearEagle <- length(yearEagle)
  
  # # Eagle days included up to myDay
  #     151 is subtracted for use in logistic model
  #      #  I.e., June 1 = 1, June 2 = 2 ....
  # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])-147 *******
  
  dayEagle <- unique(Eagle_hist$Day[Eagle_hist$Day <=(myDay) & 
                                      Eagle_hist$Day>= startDayPSS])
  
  # Number of days used
  n_dayEagle <- length(dayEagle)
  
  # Eagle daily passage estimate for days up to myDay for the year of interest (myYear)
  # curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay ]) **********
  curr_Eagle_vect <- (Eagle_hist$count[Eagle_hist$Year == myYear &
                                         Eagle_hist$Day <= myDay & 
                                         Eagle_hist$Day>=startDayPSS])
  
  curr_Eagle <- array(data = curr_Eagle_vect,
                      dim = length(curr_Eagle_vect))
  
  # Number of days used
  # Need this so stan will accept a 0 value for curr_pss
  # if(length(curr_Eagle)==1 & curr_Eagle[1]<1){
  #   curr_Eagle <- numeric(0)
  # }
  
  # Number of days used
  n_curr_Eagle <- length(curr_Eagle)
  
  # Cumulative Eagle counts for myYear
  cum_curr_Eagle <- cumsum(curr_Eagle)
  
  # Location pointers to index years in common for variance calculations
  loc.eagle.years <- which(yearPSS %in% yearEagle)
  
  loc.PF.years.eagle <- which(yearEagle %in% yearPF)
  
  loc.PF.years.PSS <- which(yearPSS %in% yearPF)
  
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
                                          # Eagle_hist$Year <= 2021 &
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
  
  if(myDay >= 178){
    # Loop to fill matix by year and day
    for (y in 1:n_yearEagle){
      for (d in SE:n_dayPSS) {
        
        
        Eagle_mat[d,y] <- count_vect_Eagle[counter]
        
        counter = counter+1
        
      }
    }
  }
  # Fill in the rest of the matrix with zeros
  Eagle_mat[is.na(Eagle_mat)]<-0
  
  # Cumulative count matrix for historic years
  Eagle_cum_hist_mat <- apply(X = Eagle_mat,MARGIN = 2,FUN = cumsum)
  
  # Cumulative Counts by year for plotting
  cumEagle <- apply(Eagle_mat,MARGIN = 2,sum)
  
  # Calculate the proportion for day D for each year Y
  prop_Eagle <- cumEagle/totalEOS[loc.eagle.years]
  
  # Mean proportion
  mean_prop_Eagle <- mean(prop_Eagle)
  
  # Vector containing avg GSI proportions for each day ACROSS ALL YEARS (for versions 3.0,3.1,3.2,3.3/ PSSreg_GSI) ####
  meanGSI_vect <- vector(length = n_dayPSS)
  
  names(meanGSI_vect) <- dayPSS
  
  # Vector of GSI sd
  sdGSI_vect <- vector(length = n_dayPSS)
  
  names(sdGSI_vect) <- dayPSS
  
  counter <- 1
  for (d in dayPSS) {
    # d = 175
    meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d &
                                                       GSI_by_year$endday >=d &
                                                       GSI_by_year$year != myYear &
                                                       GSI_by_year$year != 2013])
    
    sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d &
                                                   GSI_by_year$endday >=d &
                                                   GSI_by_year$year != myYear &
                                                   GSI_by_year$year != 2013])
    
    counter <- counter+1
  }
  
  # No GSI data for day 148 and 149. Use day 150 for these values
  meanGSI_vect[1:2] <- meanGSI_vect[3]
  
  sdGSI_vect[1:2] <- sdGSI_vect[3]
  
  # Mean across days for all season
  meanGSI_vect_all <- vector(length = n_dayPSS_all)
  
  names(meanGSI_vect_all) <- dayPSS_all
  
  sdGSI_vect_all <- vector(length = n_dayPSS_all)
  
  names(sdGSI_vect_all) <- dayPSS_all
  
  counter <- 1
  for (d in dayPSS_all) {
    # d = 175
    meanGSI_vect_all[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
    
    sdGSI_vect_all[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
    
    counter <- counter+1
  }
  
  # No GSI data for day 148 and 149. Use day 150 for these values
  meanGSI_vect_all[1:2] <- meanGSI_vect_all[3]
  
  sdGSI_vect_all[1:2] <- sdGSI_vect_all[3]
  
  # No GSI data for days 249:252
  #  add the last day available GSI information to these days
  meanGSI_vect_all[is.na(meanGSI_vect_all)] <- meanGSI_vect_all["248"]
  
  sdGSI_vect_all[is.na(sdGSI_vect_all)] <- sdGSI_vect_all["248"]
  
  # Mean GSI by strata & dates by mean start and end of stata dates (Not used)#################
  meanStartDay <-GSI_by_year %>%
    filter(year != 2013 & year != myYear) %>% 
    group_by(stratum) %>% 
    summarise("meanStartDay" = mean(startday)) %>% 
    as.data.frame()
  
  GSI_mean_by_strata <-GSI_by_year %>% 
    filter(year != 2013 & year != myYear) %>% 
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
  
  GSI_avg_vect_all <- GSI_avg[1:n_dayPSS_all]
  
  GSI_sd <- c(rep(GSI$stratumSD[GSI$stratum == 1], 
                  times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
              rep(GSI$stratumSD[GSI$stratum == 2], 
                  times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
              rep(GSI$stratumSD[GSI$stratum == 3],
                  times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))
  
  GSI_sd_vect <- GSI_sd[1:n_dayPSS]
  
  GSI_sd_vect_all <- GSI_sd[1:n_dayPSS_all]
  
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
  # Create matrix with dimensions myday and myYear
  PSS_mat_all_adj <- matrix(nrow = (length(startDayPSS:endDayPSS )),
                            #start 1996 to account for missing year 1996
                            ncol = n_yearPSS) 
  
  # Give names to matrix for accounting
  colnames(PSS_mat_all_adj) <- yearPSS
  
  rownames(PSS_mat_all_adj) <- c(startDayPSS:endDayPSS)
  
  # Use loop to populate matrix with counts by days
  # Set counter
  counter <-1
  
  for (y in 1:n_yearPSS){
    for (d in 1:(length(startDayPSS:endDayPSS))) {
      
      
      PSS_mat_all_adj[d,y] <- PSS_mat_all[d,y] * meanGSI_vect_all[d]
      
    }
  }
  
  # Cumulative GSI adjusted historic counts
  cumPSS_adj <- apply(PSS_mat_adj, 2, sum)
  
  
  adj_curr_PSS <- vector(length = n_dayPSS)
  for (d in 1:n_dayPSS) {
    
    
    adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
  }
  
  # GSI Beta parameters
  meanpropCAN <- meanGSI_vect
  
  sd_meanpropCAN <- sdGSI_vect
  
  # paramA <- vector(length = n_dayPSS)
  # paramB <- vector(length = n_dayPSS)
  # 
  # for(x in 1:n_dayPSS){
  #   # x = 1
  #   paramA[x] <- ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;
  #   
  #   paramB[x] <- paramA[x] * ((1/meanpropCAN[x])-1);
  #   
  # }
  
  # if(logistic == TRUE){
  
  # Logistic parameters for informative prior (PSSlogistic_ESprop)
  ps_m <- prior.df.log$mid[prior.df.log$year != myYear &
                             prior.df.log$year >= startYearPSS &
                             prior.df.log$year <= 2022]
  
  ps_s <- prior.df.log$sd[prior.df.log$year != myYear&
                            prior.df.log$year >= startYearPSS &
                            prior.df.log$year <= 2022]
  
  ps_alpha_log <- prior.df.log$alpha[prior.df.log$year != myYear&
                                       prior.df.log$year >= startYearPSS &
                                       prior.df.log$year <= 2022]
  
  n_ps_m <- length(ps_m)
  n_ps_s <- length(ps_s)
  n_ps_alpha_log <- length(ps_alpha_log)
  
  # Normal distribution parameters for informative priors (PSSnormal_ESprop)
  ps_mu <- prior.df.norm$mid[prior.df.norm$year != myYear &
                               prior.df.norm$year >=startYearPSS &
                               prior.df.norm$year <= 2022]
  
  ps_sd <- prior.df.norm$sd[prior.df.norm$year != myYear&
                              prior.df.norm$year >=startYearPSS &
                              prior.df.norm$year <= 2022]
  
  ps_alpha_norm <- prior.df.norm$alpha[prior.df.norm$year != myYear&
                                         prior.df.norm$year >=startYearPSS &
                                         prior.df.norm$year <= 2022]
  
  n_ps_mu <- length(ps_mu)
  n_ps_sd <- length(ps_sd)
  n_ps_alpha_norm <- length(ps_alpha_norm)
  
  
  # PSS SD (Uncertainty)
  tempSD <- PSS_sd[PSS_sd$Day<= myDay &
                     PSS_sd$Year != myYear &
                     PSS_sd$Year >= startYearPSS,] %>% 
    group_by(Year) %>% 
    summarise(sd = sqrt(sum(Var)))
  
  
  full_tempSD <- data.frame("Year" = yearPSS, "sd" = 0)
  
  full_tempSD$sd <- tempSD$sd[match(full_tempSD$Year,tempSD$Year)]
  
  full_tempSD$sd[ is.na(full_tempSD$sd)] <- 0
  
  PSS_year_sd <- full_tempSD$sd
  
  names(PSS_year_sd) <- yearPSS
  
  curr_PSS_year_sd <- sqrt(sum(PSS_sd$Var[PSS_sd$Day<= myDay &
                                            PSS_sd$Year == myYear] ))
  
  # # SST Info
  # may.sst.hist <- norton.sst$month.mean[norton.sst$month == 5 &
  #                                         norton.sst$year != 1996 &
  #                                         norton.sst$year >= 2000 &
  #                                         # norton.sst$year >= startYearPSS &
  #                                         norton.sst$year != myYear]
  # 
  # yearSST <- norton.sst$year[norton.sst$month == 5 &
  #                              norton.sst$year != 1996 &
  #                              norton.sst$year >= 2000 &
  #                              # norton.sst$year >= startYearPSS &
  #                              norton.sst$year != myYear]
  # 
  # n_yearSST <- length(yearSST)
  # 
  # may.sst.curr <-norton.sst$month.mean[norton.sst$month == 5 &
  #                                        norton.sst$year == myYear]
  # Point SST 2000 to 2022 for may
  # may.sst.hist <- norton.sst$sst[norton.sst$year != myYear]
  # 
  # yearSST <- norton.sst$year[norton.sst$year != myYear]
  # 
  # names(may.sst.hist) <- yearSST
  # 
  # n_yearSST <- length(yearSST)
  # 
  # # Current year sst
  # may.sst.curr <-norton.sst$sst[norton.sst$year == myYear]
  # 
  # n_sst <- length(may.sst.hist)
  # 
  # # Emmonak Air
  # # Air Temp Data
  # # NOTE that 2020 is interpolated with sst airtemp relationship
  # april.air.hist <- Emmonak_Air_Temp$air.mean[Emmonak_Air_Temp$year >= 2000 &
  #                                               Emmonak_Air_Temp$year != myYear]
  # 
  # yearAir <- Emmonak_Air_Temp$year[Emmonak_Air_Temp$year >= 2000 &
  #                                    Emmonak_Air_Temp$year != myYear ]
  # 
  # names(april.air.hist) <- yearAir
  # 
  # n.year.air <- length(yearAir)
  # 
  # april.air.curr <- Emmonak_Air_Temp$air.mean[Emmonak_Air_Temp$year == myYear ]
  
  # PSS timing deviations###################################################
  
  # Get vector of historic midpoints for all years excluding myYear
  # hist.midpoint <- logistic.all$mid[logistic.all$year != myYear &
  #                                     logistic.all$year >= 2000]
  # 
  # # Name each entry with the year
  # names(hist.midpoint) <-logistic.all$year[logistic.all$year != myYear &
  #                                            logistic.all$year >= 2000]
  # 
  # # Current years midpoint
  # curr.midpoint <- logistic.all$mid[logistic.all$year == myYear]
  # 
  # # Get historic mean mp
  # avg_midpoint <- mean(hist.midpoint)
  # 
  # # Calculate historic deviations from the mean mp
  # deviations.mp <- hist.midpoint - avg_midpoint
  # 
  # deviations.mp.curr <- curr.midpoint - avg_midpoint
  # 
  # # Pointer for years used in retrospective testing
  # loc.devMP.allYears <- which(yearSST %in% yearPSS)
  
  
  
  # Stan Model Call ######################################
  if(logistic==TRUE){
    
    # Logistic model intits list
    init_fn <-function(){
      
      list("ps_alpha_curr" = runif(1,1e5,11e4),
           "ps_mid_curr" = runif(1,168,170),
           "ps_shape_curr" = runif(1,2,3),
           "ps_shape_hist" = runif(n_yearPSS,2,3),
           "ps_mid_hist" = runif(n_yearPSS,168,170),
           "ps_alpha_hist"=runif(n_yearPSS,1e5,11e4),
           "sigma" = runif(1,0.1,0.5),
           "sigma_hist" = runif(n_yearPSS,0,1),
           "beta" = runif(1,0.3,0.5),
           "sigma_reg" = runif(1,0,1),
           "alpha" = runif(1,500,1500)
           
      )}
    # 
    
    # 
    inits_ll <- list(init_fn() , init_fn(), init_fn(), init_fn())
  }
  
  if(normal == TRUE){
    
    # Normal Distribution Inits List
    init_fn <-  function(){
      list(    "ps_alpha_curr" = runif(1,10e4,11e4),
               "ps_mu_curr" = runif(1,174,176),
               "ps_sd_curr" = runif(1,5,7),
               "ps_alpha_hist" = runif(n_yearPSS,10e4,11e4),
               "ps_mu_hist" = runif(n_yearPSS,174,176),
               "ps_sd_hist" = runif(n_yearPSS,5,7),
               "sigma" = runif(1,2,3),
               "sigma_hist" = runif(n_yearPSS,2,3),
               "beta_sst" = runif(1,-2,-1),
               "sigma_env" = runif(1,2.5,3),
               "alpha_env" = runif(1,177,178)
      )
    }
    
    
    # List of list for initial values in stan model
    inits_ll <- list(init_fn(), init_fn(), init_fn(),init_fn())
  }
  
  if(normal == FALSE & logistic == FALSE){
    #Used for all regression based models except below
    # inits <- function(){ 
    #   list( alpha = runif(1,150000,200000),
    #         beta = runif(1,0,0.5),
    #         sigma = runif(1,0,2),
    #         ln_phi = runif(1,-1,1),
    #         propCAN_logit = runif(n_dayPSS,0,1)
    #   )
    # }
    # MMSST model with beta combination of gamma and tau
    # inits <-  function(){
    #   list(    "tau" = runif(1,0.35,0.5),
    #            "gamma" = runif(1,0,.05),
    #            "sigma" = runif(1,0.2,0.5),
    #            "alpha" = runif(1,10000,12000)
    #   )
    # }
    
    inits_ll <- list(inits(), inits(), inits(), inits())
    
  }
  
  fit <- stan(file = file.path(dir.stan,paste("ADFG Yukon Inseason Forecast ", model.version ,".stan", sep = "")),
              data = list("PSS_mat"=PSS_mat,
                          "n_totalEOS"=n_totalEOS,
                          "totalEOS"=totalEOS, 
                          "n_dayPSS"=n_dayPSS,
                          "dayPSS"=dayPSS, 
                          "n_yearPSS"=n_yearPSS,
                          "yearPSS"=yearPSS,
                          "n_curr_PSS"=n_curr_PSS,
                          "curr_PSS"=curr_PSS,
                          "Pf" = pf, 
                          # "Pf_milo" = pf_milo,
                          "Pf_sigma" = pf_sigma,
                          # "meanpropCAN_all" = meanGSI_vect_all,
                          # "sd_meanpropCAN_all" = sdGSI_vect_all,
                          # "meanpropCAN" = meanGSI_vect,
                          # "sd_meanpropCAN" = sdGSI_vect,
                          "mean_adj_PSS_mat"=PSS_mat_adj,
                          "mean_adj_curr_PSS"= adj_curr_PSS,
                          "meanpropCAN_all"= GSI_avg_vect_all,
                          "sd_meanpropCAN_all" = GSI_sd_vect_all,
                          "meanpropCAN" = GSI_avg_vect,
                          "sd_meanpropCAN" = GSI_sd_vect,
                          "ps_m"= ps_m,
                          "ps_s"= ps_s,
                          "ps_alpha_log" = ps_alpha_log,
                          "n_ps_m"= n_ps_m,
                          "n_ps_s"= n_ps_s,
                          "n_ps_alpha_log" = n_ps_alpha_log,
                          "ps_mu"= ps_mu,
                          "ps_sd"= ps_sd,
                          "ps_alpha_norm"= ps_alpha_norm,
                          "n_ps_mu" = n_ps_mu,
                          "n_ps_sd" = n_ps_sd,
                          "n_ps_alpha" = n_ps_alpha_norm,
                          "pss_days_all" = dayPSS_all,
                          "dayPSS_all" = dayPSS_all,
                          "n_pss_days_all" = n_dayPSS_all,
                          "n_dayPSS_all" = n_dayPSS_all,
                          # "PSS_mat_adj" = PSS_mat_adj,
                          "myDay" = myDay,
                          "PSS_mat_all" = PSS_mat_all,
                          "cum_PSS_mat" = cum_PSS_mat,
                          "cum_curr_PSS" = cum_curr_PSS,
                          "Eagle_mat" = Eagle_mat, 
                          "n_dayEagle" = n_dayEagle,
                          "dayEagle" = dayEagle, 
                          "n_yearEagle" = n_yearEagle,
                          "yearEagle" = yearEagle,
                          "n_curr_Eagle" = n_curr_Eagle,
                          "curr_Eagle" = curr_Eagle,
                          "loc_pf_years_PSS" = loc.PF.years.PSS,
                          "loc_pf_years_Eagle" = loc.PF.years.eagle,
                          "loc_eagle_years" = loc.eagle.years,
                          "loc_allDays_myDay" = loc.allDays.myDay,
                          "n_yearsPF" = n_yearsPF,
                          "PSS_year_sd" = PSS_year_sd,
                          "curr_PSS_year_sd" = curr_PSS_year_sd,
                          "mean_propEagle" = mean_prop_Eagle,
                          # "may_sst_hist" = may.sst.hist,
                          # "may_sst_curr" = may.sst.curr,
                          # "n_sst" = n_sst,
                          # "n_yearSST" = n_yearSST,
                          # "loc_devMP_allYears" = loc.devMP.allYears,
                          "PSS_mat_prop_all" = PSS_mat_prop_all
                          # "dev_hist" = deviations.mp,
                          # "dev_curr" = deviations.mp.curr,
                          # "hist_MP" = hist.midpoint,
                          # "curr_MP" = curr.midpoint,
                          # "n_yearAir" = n.year.air,
                          # "april_air_hist" = april.air.hist,
                          # "april_air_curr" = april.air.curr
                          
              ),
              init = inits_ll,
              chains = n.chains,
              # chains = 1,
              # verbose = T,
              iter = n.iter, 
              thin = n.thin, 
              cores = n.chains,
              control = list(max_treedepth = 15, #Adjust for stan errors
                             adapt_delta = 0.85
              )
  )
  
  
  # Record any divergent transistions
  divergent <- get_divergent_iterations(fit)
  n <- sum(divergent)
  N <- length(divergent)
  
  div.trans <- n/N
  
  # Extract parameter estimates
  pars <- rstan::extract(fit)
  
  # Extract fit summary for saving
  fit.summary <- summary(fit)
  if(savefit == FALSE){
    outs <- list(
      "cumPSS" = cumPSS,
      "cum_curr_PSS" = cum_curr_PSS,
      "dayPSS" = dayPSS,
      "curr_PSS" = curr_PSS,
      "cumPSS_adj"= cumPSS_adj,
      "cumEagle" = cumEagle,
      "cumPSS_all" = cumPSS_all,
      "pars" = pars, 
      "totalEOS" = totalEOS,
      "summary" = fit.summary,
      "myYear" = myYear,
      "version" = model.version,
      "myDay" = myDay,
      "pss_days_all" = dayPSS_all,
      "cum_PSS_mat" = cum_PSS_mat,
      "PSS_mat" = PSS_mat,
      "PSS_mat_all_adj" = PSS_mat_all_adj,
      "PSS_mat_all" = PSS_mat_all,
      "yearPSS" = yearPSS,
      "myDay" = myDay,
      "myYear" = myYear,
      "divergent_trans" = div.trans
    )
  }else{
    outs <- list("cumPSS" = cumPSS,
                 "cum_curr_PSS" = cum_curr_PSS,
                 "dayPSS" = dayPSS,
                 "curr_PSS" = curr_PSS,
                 "cumPSS_adj"= cumPSS_adj,
                 "cumEagle" = cumEagle,
                 "cumPSS_all" = cumPSS_all,
                 "pars" = pars, 
                 "totalEOS" = totalEOS,
                 "summary" = fit.summary,
                 "myYear" = myYear,
                 "version" = model.version,
                 "myDay" = myDay,
                 "pss_days_all" = dayPSS_all,
                 "cum_PSS_mat" = cum_PSS_mat,
                 "PSS_mat" = PSS_mat,
                 "PSS_mat_all_adj" = PSS_mat_all_adj,
                 "PSS_mat_all" = PSS_mat_all,
                 "yearPSS" = yearPSS,
                 "myDay" = myDay,
                 "myYear" = myYear,
                 "divergent_trans" = div.trans,
                 "fit" = fit,
                 "dayPSS_all" = dayPSS_all)
  }#Savefit for false
}



# Funciton for plotting retro outputs ######################################################
# 
# outPlots <- function(outputList,  # Model object 
#                      CAN_hist,    # EOS run size df
#                      GSI = FALSE, # True for plotting with GSI
#                      Retrospective = FALSE, # True for plotting retro stats
#                      eagle = FALSE){        # True if Eagle is incorporated in model
#   outputPlots <- list()
#   
#   
#   # i =1
#   # Use loop output list from above to call plots and save them to a list ############
#   
#   # Extract relevant list objects for plotting
#   pars <- outputList$pars
#   fit <- outputList$fit
#   cumPSS <- outputList$cumPSS
#   cumPSS_adj <- outputList$cumPSS_adj
#   cumEagle <- outputList$cumEagle
#   totalEOS <- outputList$totalEOS
#   trueCan <- as.double(CAN_hist$can.mean[CAN_hist$Year == outputList$myYear]/1000)
#   myD <- outputList$myDay
#   myY <- outputList$myYear
#   Ver <- outputList$version
#   
#   # Get the quantiles for predPSS from model output
#   quant.predPSS <- apply(X = pars$predPSS, 
#                          MARGIN = 2, 
#                          FUN = quantile, 
#                          probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
#   
#   
#   # If statement to plot without retrospective testing
#   #  When Retrospective = FALSE, then the realized Canadian-origin line is omitted
#   if(Retrospective == TRUE){
#     # If statement to plot with GSI adjusted if TRUE
#     if(GSI == FALSE){
#       
#       pred.plot.df <- data.frame(cumPSS,
#                                  totalEOS,
#                                  t(quant.predPSS))
#       names(pred.plot.df) <- c("cumPSS", "totalEOS","low95","low50","median",
#                                "up50","up95")
#       
#       
#       
#       zz <- ggplot(pred.plot.df, aes(x = cumPSS/1000, y = totalEOS/1000 ))+
#         geom_point() +
#         geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
#         geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
#         geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
#         labs(x = "Cummulative PSS Chinook Passage (1000's)",
#              y = "Total Reconstructed EOS Chinook Passage (1000's)")+
#         theme(text = element_text(size = 20, family = "serif"),
#               panel.grid.major =element_line(),    #strip major gridlines
#               panel.grid.minor = element_blank(),    #strip minor gridlines
#               axis.ticks = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(color = "black"))+
#         scale_fill_colorblind(name = "")+
#         scale_color_colorblind(name = "")
#     } else{
#       pred.plot.df <- data.frame(cumPSS_adj,
#                                  totalEOS,
#                                  t(quant.predPSS))
#       names(pred.plot.df) <- c("cumPSS_adj", "totalEOS","low95","low50","median",
#                                "up50","up95")
#       
#       
#       
#       zz <- ggplot(pred.plot.df, aes(x = cumPSS_adj/1000, y = totalEOS/1000 ))+
#         geom_point() +
#         geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
#         geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
#         geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
#         labs(x = "Cummulative Canadian-adj PSS Chinook Passage (1000's)",
#              y = "Total Reconstructed EOS Chinook Passage (1000's)")+
#         theme(text = element_text(size = 20, family = "serif"),
#               panel.grid.major =element_line(),    #strip major gridlines
#               panel.grid.minor = element_blank(),    #strip minor gridlines
#               axis.ticks = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(color = "black"))+
#         scale_fill_colorblind(name = "")+
#         scale_color_colorblind(name = "")
#       
#     }
#     
#     
#     
#     # Extract estimates for prior, likelihood, and posterior 
#     df.run <- data.frame("par" = "Runsize", "value" = pars$RunSize)
#     df.prior <- data.frame("par" = "Prior", "value" = pars$prior_pf)
#     df.pssPred <- data.frame("par" = "PSS_Pred", "value" = pars$curr_predPSS)
#     df.postPredPss <- data.frame("par" = "PSS_post_pred", "value" = pars$post_curr_predPSS)
#     
#     # Bind into one dataframe with 2 columns
#     df.comb <- rbind(df.prior,df.postPredPss,df.run)
#     # df.comb$par <- levels(df.comb$par, levels = c("Prior","PSS_post_pred","Runsize"))
#     # Density plots comparing pf, linear prediction, and posterior estimate
#     postDense <-ggplot(df.comb, aes(x = value/1000, fill = par))+
#       geom_density(alpha = .65)+
#       geom_vline( aes(xintercept = trueCan,
#                       col = "True Runsize"),
#                   linetype = 2,
#                   size = 1.5)+
#       # ylab("Relative Probability")+
#       ylab("")+
#       xlab("")+
#       # ggtitle(paste("Version = ",Ver,"\n Day = ",myD,"\n Year = ",myY, sep = ""))+
#       # xlab("1000's of Chinook Salmon")+
#       coord_cartesian(xlim = c(0,200),
#                       ylim = c(0,.035))+
#       scale_fill_colorblind(name = "", 
#                             labels = c( "Preseason Forecast (Prior)", 
#                                         "PSS Prediction","Runsize Projection"))+
#       theme_tidybayes()+
#       theme(text = element_text(size = 20, family = "serif"),plot.title = element_text(size = 12))+
#       scale_color_discrete(name = "")
#     
#   }else{ # Else statement for Retrospective argument
#     
#     # If statement to plot with GSI adjusted if TRUE
#     if(GSI == FALSE){
#       
#       pred.plot.df <- data.frame(cumPSS,
#                                  totalEOS,
#                                  t(quant.predPSS))
#       names(pred.plot.df) <- c("cumPSS", "totalEOS","low95","low50","median",
#                                "up50","up95")
#       
#       
#       
#       zz <- ggplot(pred.plot.df, aes(x = cumPSS/1000, y = totalEOS/1000 ))+
#         geom_point() +
#         geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
#         geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
#         geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
#         labs(x = "Cummulative PSS Chinook Passage (1000's)",
#              y = "Total Reconstructed EOS Chinook Passage (1000's)")+
#         theme(text = element_text(size = 20, family = "serif"),
#               panel.grid.major =element_line(),    #strip major gridlines
#               panel.grid.minor = element_blank(),    #strip minor gridlines
#               axis.ticks = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(color = "black"))+
#         scale_fill_colorblind(name = "")+
#         scale_color_colorblind(name = "")
#     } else{
#       pred.plot.df <- data.frame(cumPSS_adj,
#                                  totalEOS,
#                                  t(quant.predPSS))
#       names(pred.plot.df) <- c("cumPSS_adj", "totalEOS","low95","low50","median",
#                                "up50","up95")
#       
#       
#       
#       zz <- ggplot(pred.plot.df, aes(x = cumPSS_adj/1000, y = totalEOS/1000 ))+
#         geom_point() +
#         geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
#         geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
#         geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
#         labs(x = "Cummulative Canadian-adj PSS Chinook Passage (1000's)",
#              y = "Total Reconstructed EOS Chinook Passage (1000's)")+
#         theme(text = element_text(size = 20, family = "serif"),
#               panel.grid.major =element_line(),    #strip major gridlines
#               panel.grid.minor = element_blank(),    #strip minor gridlines
#               axis.ticks = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(color = "black"))+
#         scale_fill_colorblind(name = "")+
#         scale_color_colorblind(name = "")
#       
#     }
#     
#     
#     
#     # Extract estimates for prior, likelihood, and posterior 
#     df.run <- data.frame("par" = "Runsize",
#                          "value" = pars$RunSize)
#     
#     # PF (prior)
#     df.prior <- data.frame("par" = "Prior", 
#                            "value" = pars$prior_pf)
#     
#     # PSS prediction density
#     df.postPredPss <- data.frame("par" = "PSS Prediction", 
#                                  "value" = pars$post_curr_predPSS)
#     
#     # Eagle prediction density
#     df.postPredEagle <- data.frame("par" = "Eagle Prediction", 
#                                    "value" = pars$post_curr_predEagle)
#     
#     # Bind into one data frame with 2 columns for plotting
#     df.comb <- rbind(df.prior,
#                      df.postPredPss,
#                      df.postPredEagle,
#                      df.run)
#     # df.comb$par <- levels(df.comb$par, levels = c("Prior","PSS_post_pred","Runsize"))
#     # Density plots comparing pf, linear prediction, and posterior estimate
#     postDense <-ggplot(df.comb, aes(x = value/1000, fill = par))+
#       geom_density(alpha = .65)+
#       geom_vline( aes(xintercept = median(pars$RunSize)/1000,
#                       col = "Median Projection"),
#                   linetype = 2,
#                   size = 1.5)+
#       # ylab("Relative Probability")+
#       ylab("")+
#       xlab("")+
#       # ggtitle(paste("Version = ",Ver,"\n Day = ",myD,"\n Year = ",myY, sep = ""))+
#       # xlab("1000's of Chinook Salmon")+
#       coord_cartesian(xlim = c(0,200))+
#       scale_fill_colorblind(name = "", 
#                             labels = c( "Preseason Forecast (Prior)", 
#                                         "PSS Prediction","Runsize Projection"))+
#       theme_tidybayes()+
#       theme(text = element_text(size = 20, family = "serif"),plot.title = element_text(size = 12))+
#       scale_color_discrete(name = "")
#     
#     eyeDense <- ggplot(df.comb, aes(x = value/1000, y =fct_rev( par)))+
#       stat_eye(aes(fill = after_stat(level)), .width = c(.8,.95,1), alpha = .5)+
#       geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
#                      color = "End of Season Runsize"),
#                  linetype = 2,
#                  size = 1)+
#       coord_cartesian(xlim = c(20,300))+
#       # ylim(c(0,500))
#       labs(x = "1000's of Chinook Salmon", 
#            fill = "Parameters", 
#            y = "Parameter")+
#       # scale_fill_discrete(name = "Parameters", 
#       # labels = c( "Preseason Forecast (Prior)", 
#       # "PSS Prediction","Runsize"))+
#       scale_x_continuous()+
#       # scale_fill_colorblind(name = "", 
#       #                       labels = c( "Preseason Forecast (Prior)",
#       #                                   "PSS Prediction",
#       #                                   "Runsize"))+
#       theme_tidybayes()+
#       theme(text = element_text(size = 16, family = "serif"),
#             legend.position = "top",
#             panel.grid.major =element_line(),    #strip major gridlines
#             panel.grid.minor = element_blank(),    #strip minor gridlines
#             axis.ticks = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(color = "black"))+
#       scale_color_discrete(name = "")+
#       guides(fill = guide_legend(override.aes = list(size = 5)))
#   }
#   
#   
#   
#   # Eagle regression plot
#   if(eagle == TRUE){
#     
#     quant.predEagle <- apply(X = pars$predEagle, 
#                              MARGIN = 2, 
#                              FUN = quantile, 
#                              probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
#     
#     eagle.pred.plot.df <- data.frame(cumEagle,
#                                      totalEOS[10:max(length(totalEOS))],
#                                      t(quant.predEagle))
#     names(eagle.pred.plot.df) <- c("cumEagle", "totalEOS","low95","low50","median",
#                                    "up50","up95")
#     
#     
#     
#     zz.Eagle <- ggplot(eagle.pred.plot.df, aes(x = cumEagle/1000, y = totalEOS/1000 ))+
#       geom_point() +
#       geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
#       geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
#       geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
#       labs(x = "",
#            y = "")+
#       theme(text = element_text(size = 20, family = "serif"),
#             panel.grid.major =element_line(),    #strip major gridlines
#             panel.grid.minor = element_blank(),    #strip minor gridlines
#             axis.ticks = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(color = "black"))+
#       scale_fill_colorblind(name = "")+
#       scale_color_colorblind(name = "")
#     
#     
#     
#     
#     # Posterior EOS Can-orig run size estimate
#     df.run <- data.frame("par" = "Runsize",
#                          "value" = pars$RunSize)
#     
#     # PF (prior)
#     df.prior <- data.frame("par" = "Prior", 
#                            "value" = pars$prior_pf)
#     
#     # PSS prediction density
#     df.postPredPss <- data.frame("par" = "PSS_post_pred", 
#                                  "value" = pars$post_curr_predPSS)
#     
#     # Eagle prediction density
#     df.postPredEagle <- data.frame("par" = "Eagle_post_pred", 
#                                    "value" = pars$post_curr_predEagle)
#     
#     # Bind into one data frame with 2 columns for plotting
#     df.comb <- rbind(df.prior,
#                      df.postPredPss,
#                      df.postPredEagle,
#                      df.run)
#     str(df.comb)
#     df.comb$par <- as.factor(df.comb$par)
#     
#     df.comb$par <- relevel(df.comb$par,"Prior")
#     
#     # Density plots comparing pf, linear prediction, and posterior estimate
#     Eagle.Dense<- ggplot(df.comb, aes(x = value/1000, fill = par))+
#       geom_density(alpha = .5)+
#       # ggtitle(paste("Density plot for",myY,", day",myD,"\n Version", model.version))+
#       theme_classic()+
#       geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myY]/1000,
#                      color = "End of Season Runsize"),
#                  linetype = 2,
#                  size = 1)+
#       
#       # xlim(c(0,300000))+
#       coord_cartesian(xlim = c(0,200),
#                       ylim = c(0,.05))+
#       # ylim(c(0,500))
#       labs(x = "", 
#            fill = "Parameters", 
#            y = "")+
#       # scale_fill_discrete(name = "Parameters", 
#       # labels = c( "Preseason Forecast (Prior)", 
#       # "PSS Prediction","Runsize"))+
#       scale_x_continuous()+
#       scale_fill_colorblind(name = "", 
#                             labels = c( "Preseason Forecast (Prior)", 
#                                         "Eagle Prediction",
#                                         "PSS Prediction",
#                                         "Runsize"))+
#       # theme_tidybayes(element_text(size = 20))+
#       theme(text = element_text(size = 20, family = "serif"),
#             panel.grid.major =element_line(),    #strip major gridlines
#             panel.grid.minor = element_blank(),    #strip minor gridlines
#             axis.ticks = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(color = "black"))+
#       scale_color_discrete(name = "")
#     
#     
#   }
#   
#   
#   
#   if(eagle == FALSE){
#     outputPlots <- list("DensPlot" = postDense, "PredPlot"= zz, "eyePlot" = eyeDense)
#   }else{
#     outputPlots <- list("DensPlot" = Eagle.Dense, "PredPlot"= zz, "EaglePlot" = zz.Eagle, )
#   }
#   
#   
#   
#   return(outputPlots)
# }


# Function for easily getting day of year for this project
myDay_func <- function(Month, Day){
  # Function to get doy for the Inseason Bayesian Projection Model
  # Enter Month as a number
  #Testing
  # Month <- 6 ; Day <- 25
  
  # Required package for function
  library(lubridate)
  
  # Df for storing dates
  date.df <- data.frame("dayofyear" = c(148:260), 
                        "date" = seq(from = as.Date("2022-05-28"), 
                                     to = as.Date("2022-9-17"), by=1) )
  # Extract month and day
  date.df$month <-  month(date.df$date)
  date.df$day <- day(date.df$date)
  
  MYDAY <- as.double(date.df$dayofyear[date.df$month == Month & date.df$day == Day])
  
  return(MYDAY)
}
