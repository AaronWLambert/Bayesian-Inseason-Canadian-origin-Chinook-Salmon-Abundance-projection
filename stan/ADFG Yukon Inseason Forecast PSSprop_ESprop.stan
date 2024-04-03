// Aaron Lambert
// 2/27/24

// Yukon Chinook Inseason Forecast Version 6.0.2

// Notes: 
//  Uses logistic fit to all years to generate the expecte proportion
//  on day D

//Components:
//  1) Pre-season forecast (Prior)
//  2) PSS logistic fit to proporitons  
//  3) Eagle simple proportion estimator
//  



data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  int <lower=0> n_yearsPF;

  // EOS Canadian abundance by year
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  vector <lower=0> [n_dayPSS]dayPSS;
  int <lower=0> n_dayPSS_all;
  vector <lower=0> [n_dayPSS_all]dayPSS_all;
  
  // Historic PSS years 
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // Matrix of PSS passage by day and year
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_all;
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_prop_all;
      // Matrix of cumulative historic PSS passage by days and years
  // matrix <lower=0> [n_dayPSS_all, n_yearPSS]cum_PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // Day the projection is made on
  real myDay;

  // Matrix of historic Eagel passage by days & years
  int<lower=0> n_yearEagle;
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year Eagle passage up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Pointer vectors for correct years
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  // int <lower=0> loc_devMP_allYears[n_yearPSS];
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  
  // mean eagle prpoporiton of total run size on myDay
  real mean_propEagle;
  
}

transformed data{
  
  // Vector to hold sum of counts for each historic year up to myDay
  vector [n_yearPSS]cumHistPSS;
  
  // Vector for the sum of each historic year to the EOS
  real cumHistPSS_all[n_yearPSS];
  
  // Sum of the current year up to myDay
  real cum_current_PSS;
  
  // matrix [n_dayPSS_all, n_yearPSS]PSS_mat_prop_logit;
  
  // Loop to get cumulative counts
    for (i in 1:n_yearPSS){
        cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
      }
  
  // Loop to get cumulative historic counts
    for (i in 1:n_yearPSS){
        cumHistPSS[i] = sum(PSS_mat[,i]);
      }
  // Sum up the current year PSS passage to myDay
  cum_current_PSS = sum(curr_PSS);
  
  
}

// 
parameters {
   
  // PSS regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // logistic curve paramaters
  real  <lower=166, upper = 182> mid;
  real  <lower=1,upper = 12> shape;
  real  <lower=0> sigma_logistic;
  
  // Log run size
  real ln_RunSize;
  
}

transformed parameters{
  
  // Historic PSS prediction
  vector [n_yearPSS]predPSS;
  
  // Current years PSS prediction
  real curr_predPSS;
  
  // Empirical sd for PSS for updating
  real sigma_predPSS;
  
  // Expected proportion for myDay (day of projection)
  real expectedProp;
  
  // Historic estimates of total PSS passage from expected proportion
  vector [n_yearPSS]est_predPSS;
  
  // This years (myYear) estimate for total PSS passage from expected proportion
  real curr_est_predPSS;
  
  // Predicted PSS historically for just years used for updating comparison
  // used to calc empircal sd
  vector [n_yearPSS]predPSS_hist;
  
  // PSS predicted proportions for fitting logistic
  vector [n_dayPSS_all]ps_prop_pred;
  
  // Cumulative Eagle counts up to myDay
  vector [n_yearEagle]cumHistEagle;
  
  // Predicted Eagle passage
  vector [n_yearEagle]predEagle;
  
  // This years Eagle passsage up to myDay
  real cum_current_Eagle;
  
  // This years predicted run size from Eagle
  real curr_predEagle;
  
  // Empirical sd for updating with Eagle est
  real sigma_predEagle;
  
  // Fit logistic to obs PSS proportions
  for(d in 1:n_dayPSS_all){
    ps_prop_pred[d] =  (1/(1 + exp((-(dayPSS_all[d] - mid)) / shape)));}

  // Regress cumulative total PSS against the EOS Can run size
  for (i in 1:n_yearPSS){
    predPSS[i] = alpha + beta * cumHistPSS_all[i];}
  
  // Calculate expected proportion for day of projection
  expectedProp = 1/(1+exp((-(myDay - mid))/ shape));

  // Estimated historic years by proportions
  for(y in 1:n_yearPSS){
  est_predPSS[y] = cumHistPSS[y] /expectedProp;}
  
  // Using the historic est total PSS passage to calc empirical variance
  predPSS_hist = alpha + beta * est_predPSS;
  
  // Empirical SD for update
  sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - log(predPSS_hist[loc_pf_years_PSS]));
  
  // This years est total PSS passage
  curr_est_predPSS = cum_current_PSS/expectedProp;
 
  // This years PSS est of run size
  curr_predPSS = alpha + beta * curr_est_predPSS;
  
  // Loop to get cumulative Eagle counts
  for (i in 1:n_yearEagle){
     cumHistEagle[i] = sum(Eagle_mat[,i]);
    }
 
  // Summ current year Eagle passage
   cum_current_Eagle = sum(curr_Eagle);
  
  // Eagle Component using observed propotions 
  // Only used when fish are observed at eagle in the projection year
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
    
  // This years Eagle est of run size
  curr_predEagle = (cum_current_Eagle+1)/mean_propEagle;
  
  // Historic est of run size
  predEagle =  (cumHistEagle+1)/mean_propEagle;
 

  // Empirical SD for Eagle update to run size
 sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));}
  
  
  
  
  
 
  // Print statements to diagnose issues
  // Comment out when resolved :
  
  // print("ln_predPSS",ln_predPSS);
  
}

model {
  
  // Prior on SD for logistic fitting
  sigma_logistic ~ normal(0,100); // 
  
  // Regression priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma ~ normal(0,5);
  
  // Run size drawn form preseason forcast prior dist
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // Logistic curve fitting 
  for(d in 1:n_dayPSS_all){
    for(y in 1:n_yearPSS){

    PSS_mat_prop_all[d,y]  ~ normal(ps_prop_pred[d], sigma_logistic);}}
    
  // Total PSS vs Run size regression fitting
  for(i in 1:n_yearPSS){
    log(totalEOS[i]) ~ normal(log(predPSS[i]),sigma);
  }

  // Update PF with PSS est
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), sigma_predPSS);
  
   // Update PF with Eagle est
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}

}

generated quantities{
  
  // Generated prediction intervals for analysis
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  real RunSize;
  
  ln_prior_pf = normal_rng(Pf, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
  // Bring it back to reality
  RunSize = exp(ln_RunSize);
  
  
 }



