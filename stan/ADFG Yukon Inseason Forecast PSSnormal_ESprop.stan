// Yukon King Inseason Forecast Version 4.0
// Notes: fitting a normal curve to total PSS passage to predict EOS Can-origin Chinook abundance


// Components:
//  1) Pre-season forecast  (Prior)
//  2) PSS Normal Curvve fitting to arival dist
//  3) EOS PSS regression
//  

// Next to address:
//  Sigma pred needs to be reflective of curve fitting, not regression with actual counts



data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  int <lower=0> n_yearsPF;

  // EOS Canadian counts by year up to myyear -1
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  vector <lower=0> [n_dayPSS]dayPSS;
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // All the days in PSS season
  int <lower=0> n_dayPSS_all;
  vector <lower=0> [n_dayPSS_all]dayPSS_all;
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // Estimated Mean, sd, and alphas for normal curve
  int  <lower=0> n_ps_mu;
  int  <lower=0> n_ps_sd;
  int  <lower=0> n_ps_alpha;
  
  real <lower=0> ps_mu[n_ps_mu]; 
  real <lower=0> ps_sd[n_ps_sd];
  real <lower=0> ps_alpha_norm[n_ps_alpha];
  
  // Historic PSS years excluding myYear
  int<lower=0> n_yearEagle;
  
  // Matrix of historic counts by days & years
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Location vector of integers
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_allDays_myDay; 
  
  // Propotion of observed fish at Eagle
  real mean_propEagle;
  
  int <lower=0> myDay;
}

// 
parameters {
   
  // Normal curve paramaters for current year
  real <lower=10e3,upper=4e5> ps_alpha_curr;
  // real <lower=0> ps_alpha_curr;
  real <lower=165,upper=181> ps_mu_curr;
  // real <lower=0> ps_mu_curr;
  real <lower=1,upper = 14> ps_sd_curr;
  // real <lower=0> ps_sd_curr;
  real <lower=0> sigma;
  
  // Normal curve parameters for histric years
  real <lower=10e3,upper = 4e5> ps_alpha_hist[n_yearPSS];
  // real <lower=0> ps_alpha_hist[n_yearPSS];
  real <lower=165,upper=181> ps_mu_hist[n_yearPSS];
  // real <lower=0> ps_mu_hist[n_yearPSS];
  real <lower=1,upper=14> ps_sd_hist[n_yearPSS];
  // real <lower=0> ps_sd_hist[n_yearPSS];
  real <lower=0> sigma_hist[n_yearPSS];
  
  // Regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma_reg;
  
  
  // 
  real ln_RunSize;
  
}

transformed parameters{
  
  // Vector for holding this years daily predicted passage
  vector [n_dayPSS_all]ps_pred_curr;
  
  // Vector for fitting to all observed days
  vector [n_dayPSS]ps_pred_curr_obs;
  
  // // Matrix for fitting curve to complete historic passage/days
  matrix [n_dayPSS_all,n_yearPSS]ps_pred_hist;
  // 
  // // Matrix for fitting curve to historic years up to myDay
  matrix [n_dayPSS, n_yearPSS]ps_pred_hist_obs;
  
  // Vector to hold cumulative sums for historic passage up to myDay 
   //Empty vector for cum sum of counts for each year to myyear -1
  vector [n_yearPSS]cumHistPSS;
  
  // Historic cumulative curve predictions to add to cumHistPSS
  
  vector [n_yearPSS]cum_predicted_histPSS;
  
    // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]predPSS;
  
  // Vector of cumulative complete PSS passage up to myYear-1
  vector [n_yearPSS]cumHistPSS_all;
  
  // PSS prediction from current year
  real curr_predPSS;
  
  // PSS historic prediction for emiricap sd for update
  vector [n_yearPSS]predPSS_hist;
  
  // Empirical sd for update
  real sigma_predPSS;
  
  real sigma_predEagle;
  
  real cum_current_Eagle;
  
  real curr_predEagle;
  
  //Empty vector for cum sum of counts for each year excluding myYear
  vector [n_yearEagle]cumHistEagle;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearEagle]predEagle;
  
  // Curve fit to current season/year of interest
  ps_pred_curr[1:n_dayPSS_all] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((dayPSS_all - ps_mu_curr)/ps_sd_curr));
  
  
  // Curve fit to days up to myDay for the current year of interest
  ps_pred_curr_obs[1:n_dayPSS] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((dayPSS - ps_mu_curr)/ps_sd_curr));
  
  // Curves fit to complete and partial historic passage by year
  for(y in 1:n_yearPSS){
    for(d in 1:n_dayPSS_all)
    
    ps_pred_hist[d,y] = ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((dayPSS_all[d] - ps_mu_hist[y])/ps_sd_hist[y]));

    for(d in 1:n_dayPSS){
      
        ps_pred_hist_obs[d,y] = ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((dayPSS[d] - ps_mu_hist[y])/ps_sd_hist[y]));
    }

  }
  
   
   // Sum up prediction from Normal curve fitting for days myDay +1 to 252
   // Used to calculate empirical sd in update
   
    for( i in 1:n_yearPSS){
      
      cumHistPSS[i] = sum(PSS_mat[,i]);
    }
    
    cum_predicted_histPSS[1:n_yearPSS] = cumHistPSS + sum(ps_pred_hist[(loc_allDays_myDay + 1):n_dayPSS_all,]);
   
    
    // Sum the observed and current year estimated passage to end of season
    cum_current_PSS = sum(curr_PSS)+sum(ps_pred_curr[(loc_allDays_myDay + 1):n_dayPSS_all]);
 
    // Sum up complete PSS histoic days for use in regression
     for( i in 1:n_yearPSS){
      cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
    }
    
    //  Calculate predPSS from aplha, beta and cumHistPSS for model section
    // for (i in 1:n_yearPSS){
    //   predPSS[i] = alpha + beta * cumHistPSS_all[i];}
    predPSS = alpha +beta * cumHistPSS_all;
    
    // Historic prediction
    predPSS_hist = alpha + beta * cum_predicted_histPSS;
    

    // Empirical SD for update  
    // sigma_pred = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS_hist[loc_pf_years_PSS]));
    sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - log(predPSS_hist[loc_pf_years_PSS]));

    // Get current years pss prediction for update
    curr_predPSS = alpha + beta * cum_current_PSS;

    // Summ current year Eagle passage
    cum_current_Eagle = sum(curr_Eagle);
       
    // Eagle Component using observed propotions
    if(mean_propEagle > 0 && cum_current_Eagle >0){
        
        // Loop to get cumulative Eagle counts
    for (i in 1:n_yearEagle){
      cumHistEagle[i] = sum(Eagle_mat[,i]);
        }
 
    // // Summ current year Eagle passage
    //    cum_current_Eagle = sum(curr_Eagle);
   
       curr_predEagle = (cum_current_Eagle + 1e-3)/mean_propEagle;
  
       predEagle =  (cumHistEagle + 1e-3)/mean_propEagle;
 

  // Empirical SD for update  
      sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));

      }// end if statement
  
  
  
 // Print statements used to debug script
 // print("ps_mu=",ps_mu);
 // print("ps_sd=",ps_sd);
 // print("alpha=",ps_alpha);
 // print("ps_pred_curr=",ps_pred_curr);
 // print("ps_pred_curr_obs=",ps_pred_curr_obs);
 // print("ps_pred_hist=",ps_pred_hist);
 // print("ps_pred_curr_alpha=",ps_pred_curr_alpha);
  // print("cum_predicted_histPSS=",cum_predicted_histPSS)
  // print("cum_current_PSS:", cum_current_PSS);
  // print("ln_predPSS_hist:",ln_predPSS_hist);
  // print("sigma_pred:",sigma_pred);
  //   print("ln_predPSS:",ln_predPSS);

}

model {
  
  // Arival Dist Priors
  ps_alpha_curr ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(45e3,5e5);
  ps_mu_curr ~ normal(mean(ps_mu), sd(ps_mu));
  // ps_mu_curr ~ uniform(172,181);
  ps_sd_curr ~ normal(mean(ps_sd), sd(ps_sd));
  sigma ~ normal(0,5);
  // sigma ~ uniform(0,3);
  
  ps_alpha_hist ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(45e3,5e5);
  ps_mu_hist ~ normal(mean(ps_mu), sd(ps_mu));
  // ps_mu_hist ~ uniform(172,181);
  ps_sd_hist ~ normal(mean(ps_sd), sd(ps_sd));
  sigma_hist ~ normal(0,5);
  // sigma_hist ~ uniform(0,3);
  
  // Regression priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma_reg ~ normal(0,10);
  
  // Preseason Forcast likelihood
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  
  // Arrival dist curve fitting likelihood
    for(d in 1:n_dayPSS){
    if(curr_PSS[d] > 0){
      
      log(curr_PSS[d] ) ~ normal(log(ps_pred_curr_obs[d] ), sigma);}
  
    for(y in 1:n_yearPSS){
      if(PSS_mat[d,y] > 0){
        
      log(PSS_mat[d,y]) ~ normal(log(ps_pred_hist_obs[d,y]), sigma_hist[y]);}
    }

  }
  
  // Regression likelihood
  log(totalEOS[1:n_yearPSS]) ~ normal(log(predPSS[1:n_yearPSS]), sigma_reg);

  // Update run size projection with PSS projection
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), sigma_predPSS);
  
  // Update likelihood for the posterior prediction with the eagle projection
  if(mean_propEagle > 0 && cum_current_Eagle >0){
  target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}

}

generated quantities{
  //
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  // Projection
  real RunSize;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf, Pf_sigma);
  
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(mean_propEagle > 0 && cum_current_Eagle >0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}

  // Exponentiate the projection
  RunSize = exp(ln_RunSize);

}



