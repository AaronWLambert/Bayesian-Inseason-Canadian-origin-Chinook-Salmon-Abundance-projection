// Yukon King Inseason Forecast Version PSSlogistic_ESprop

//Notes: 
// fitting a logistic curve to each years PSS passage to project total PSS passage


//Components:
//  1) Pre-season forecast  (Prior)
//  2) PSS logistic curve and Season total PSS regression
//  3) Eagle proportion estimator

data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  int <lower=0> n_yearsPF;

  // EOS Canadian abundacne by year
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  vector <lower=0> [n_dayPSS]dayPSS;
  
  // Historic PSS years
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // All the days in PSS season (June - August)
  int <lower=0> n_dayPSS_all;
  vector <lower=0> [n_dayPSS_all]dayPSS_all;
  
  // Matrix of historic PSS passage by days & years
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Matrix of cumulative historic PSS passage by days and years
  matrix <lower=0> [n_dayPSS, n_yearPSS]cum_PSS_mat;
  
  // Matrix of total historic PSS passage
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_all;
  
  // Current year passage up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // Cumulative current year counts up to myDay
  real<lower=0> cum_curr_PSS[n_dayPSS];
  
  // Estimated Mean, sd, and alpha for informative priors on fitting logistic curve
  int  <lower=0> n_ps_m;
  int  <lower=0> n_ps_s;
  int  <lower=0> n_ps_alpha_log;
  
  real <lower=0> ps_m[n_ps_m]; 
  real <lower=0> ps_s[n_ps_s];
  real <lower=0> ps_alpha_log[n_ps_alpha_log];
  
  // Matrix of historic Eagle passage by days & years
  int<lower=0> n_yearEagle;
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current Eagle passage up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Location vectors
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_allDays_myDay; 
  
  // Expected proportion of cumulative Eagle passage and run size
  real mean_propEagle;
  
  // Day of the projection
  int <lower=0> myDay;
}

// 
parameters {
   
  // logistic curve paramaters for current year
  // real <lower=10000, upper=1e6> ps_alpha_curr;
  real <lower=0> ps_alpha_curr;
  real <lower=0> ps_mid_curr;
  real <lower=0> ps_shape_curr;
  real <lower=0> sigma;
  
  // Logistic curve parameters for histric years
  // real <lower=10000, upper=1e6> ps_alpha_hist[n_yearPSS];
  vector <lower=0> [n_yearPSS]ps_alpha_hist;
  vector <lower=0> [n_yearPSS]ps_mid_hist;
  vector <lower=0> [n_yearPSS]ps_shape_hist;
  vector <lower=0> [n_yearPSS]sigma_hist;
  
  // PSS regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma_reg;

  // log run size
  real ln_RunSize;
  
}

transformed parameters{
  
  // This years daily predictions from logistic fitting
  vector [n_dayPSS_all]ps_pred_curr;
  
  // Fitting to just the observed days
  vector [n_dayPSS]ps_pred_curr_obs;
  
  // Predictions historically from logistic curve
  matrix [n_dayPSS_all,n_yearPSS]ps_pred_hist;
  // 
  matrix [n_dayPSS, n_yearPSS]ps_pred_hist_obs;
  
  // PSS predictions for fitting regression
  vector [n_yearPSS]predPSS;
  
  // Vector of cumulative complete PSS passage
  vector [n_yearPSS]cumHistPSS_all;
  
  // PSS prediction from current year
  real curr_predPSS;
  
  // PSS historic prediction for emirical sd for update (uses logistic estimated PSS passage)
  vector [n_yearPSS]predPSS_hist;
  
  // PSS empirical sd for update
  real sigma_predPSS;
  
  // Eagle empirical sd for update
  real sigma_predEagle;
  
  // Cumualtive eagele passage
  real cum_current_Eagle;
  
  // This years Eagle prediction
  real curr_predEagle;
  
  // Historic cumulative Eagle passage
  vector [n_yearEagle]cumHistEagle;
  
  // Historic Eagle predictions
  vector [n_yearEagle]predEagle;
  
  // Logistic curve fit and prediction for current year
  for(d in 1:n_dayPSS_all){
   ps_pred_curr[d] = ps_alpha_curr *(1/(1 + exp((-(dayPSS_all[d] - ps_mid_curr)) / ps_shape_curr)));}
   
  // Curve fit to days up to myDay for the current year of interest
  for(d in 1:n_dayPSS){
  ps_pred_curr_obs[d] = ps_alpha_curr * (1/(1 + exp((-(dayPSS[d] - ps_mid_curr))/ps_shape_curr)));}
  
  // Logistic curve fit and prediction for historic years
  for(y in 1:n_yearPSS){
    for(d in 1:n_dayPSS_all)
    ps_pred_hist[d,y] = ps_alpha_hist[y] * (1/(1 + exp((-(dayPSS_all[d] - ps_mid_hist[y]))/ps_shape_hist[y])));

    for(d in 1:n_dayPSS){
      ps_pred_hist_obs[d,y] = ps_alpha_hist[y] * (1/(1 + exp((-(dayPSS[d] - ps_mid_hist[y]))/ps_shape_hist[y])));
    }
  }

    // Sum up complete PSS passage for season total regression
     for( i in 1:n_yearPSS){
      cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
    }
    
    // PSS regression
    predPSS = alpha + beta * cumHistPSS_all;
    
    // Use alpha from fitting logistic curve to past years for empirical sd
    predPSS_hist = alpha + beta * ps_alpha_hist;
    
    // Empirical SD for update  
    sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS_hist[loc_pf_years_PSS]));
      
    // Get current years PSS prediction for update
    curr_predPSS = alpha + beta * ps_alpha_curr;
    
    // Summ current year Eagle passage
    cum_current_Eagle = sum(curr_Eagle);
    
    // Eagle Component using observed propotions
    if(mean_propEagle > 0 && cum_current_Eagle > 0){ // Only do when Eagle passage is observed
        
    // Loop to get cumulative Eagle counts
    for (i in 1:n_yearEagle){
      cumHistEagle[i] = sum(Eagle_mat[,i]);
        }
 
    // Summ current year Eagle passage
       cum_current_Eagle = sum(curr_Eagle);
       
    // Eagle predictions
       curr_predEagle = (cum_current_Eagle+1)/mean_propEagle;
  
       predEagle =  (cumHistEagle+1)/mean_propEagle;
 

    // Eagle empirical SD for update  
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
  // print("predPSS:",predPSS);
  // print("sigma_pred:",sigma_pred);
  
}

  model{
    
  // Logistic curve fitting priors
  ps_alpha_curr ~ normal(mean(ps_alpha_log), sd(ps_alpha_log));
  // ps_alpha_curr ~ uniform(30000,400000);
  ps_mid_curr ~ normal(mean(ps_m), sd(ps_m));
  // ps_mid_curr ~uniform(167,182);
  ps_shape_curr ~ normal(mean(ps_s), sd(ps_s));
  sigma ~ normal(0,5);
  
  ps_alpha_hist ~ normal(mean(ps_alpha_log), sd(ps_alpha_log));
  // ps_alpha_hist ~ uniform(30000,400000);
  ps_mid_hist ~ normal(mean(ps_m), sd(ps_m));
  // ps_mid_hist ~uniform(167,182);
  ps_shape_hist ~ normal(mean(ps_s), sd(ps_s));
  sigma_hist ~ normal(0,5);

  // Regression priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma_reg ~ normal(0,10);
  
  // log preseason Forcast run size
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // Fitting the logisitic curve
  for(d in 1:n_dayPSS){ 
    if(cum_curr_PSS[d] > 0){
      
  log(cum_curr_PSS[d] ) ~ normal(log(ps_pred_curr_obs[d]), sigma);}

    for(y in 1:n_yearPSS){
    if(cum_PSS_mat[d,y] > 0 ){
      
      log(cum_PSS_mat[d,y]) ~normal(log(ps_pred_hist_obs[d,y]), sigma_hist[y]);}
    }

  }
  
  // Fitting regression
  log(totalEOS[1:n_yearPSS])~normal(log(predPSS[1:n_yearPSS]), sigma_reg);

  // Update run size with PSS
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), sigma_predPSS);
  
  // Update likelihood for the posterior prediction with the eagle projection
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}

}

generated quantities{
  
  real RunSize;
  
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  // Exponentiate the projection
  RunSize = exp(ln_RunSize);
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);

  // Exponential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);

  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);

  post_curr_predPSS = exp(ln_post_curr_predPSS);

  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
}









