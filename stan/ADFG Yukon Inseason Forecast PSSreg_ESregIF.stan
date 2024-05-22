// Yukon King Inseason Forecast version PSSreg_ESregIF

//Notes: Regression in normal space with eagle added into likelihood only 
//       after fish first appear at Eagle in the current year


//Components:
//  1) Pre-season forecast  (Prior)
//  2) EOS-Can ~ PSS in normal space
//  3) Eagle Sonar regression added to likelihood
//  



data {
  
  // This years Preseason forecast
  real <lower=0>Pf;
  real <lower=0>Pf_sigma;
  int <lower=0> n_yearsPF;
  

  // EOS Canadian counts by year
  int <lower=0>n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int<lower=0> n_dayPSS ;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years
  int<lower=0> n_yearPSS;
  int <lower=0>yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  matrix <lower=0>[n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year counts up to myDay
  int <lower=0> n_curr_PSS;
  real <lower=0> curr_PSS[n_curr_PSS];
  
  // Number of Eagle sonar years used
  int<lower=0> n_yearEagle;
  
  // Matrix of historic Eagle Sonar counts by days & years
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year Eagle counts up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Location vector of integers
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  
}

// 
parameters {
   
  // PSS regression priors
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // Eagle Sonar regression priors
  real <lower=0> alpha_eagle;
  real <lower=0> beta_eagle;
  real <lower=0> sigma_eagle;
  
  // log posterior run size
  real ln_RunSize;
  
}

transformed parameters{
  
  // Empty vector for cum sum of PSS counts for each year
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts
  vector [n_yearPSS]predPSS;
  
  // Empty vector for cum sum of Eagle counts for each year
  real cumHistEagle[n_yearEagle];
  
  // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // Current years PSS prediction
  real curr_predPSS;
  
  // PSS empirical SD
  real sigma_predPSS;
  
  // Empty vector for predicted Eagle counts
  vector [n_yearEagle]predEagle;
  
  // Current years Eagle counts
  real cum_current_Eagle;

  // Current years predicted Eagle passage
  real curr_predEagle;
  
  // Empirical Eagle SD
  real sigma_predEagle;
  
  // Loop to get cumulative PSS counts
  for (i in 1:n_yearPSS){
     cumHistPSS[i] = sum(PSS_mat[,i]);
    }
 
  // Loop to get cumulative Eagle counts
  for (i in 1:n_yearEagle){
     cumHistEagle[i] = sum(Eagle_mat[,i]);
    }
 
  // Summ current year Eagle passage
   cum_current_Eagle = sum(curr_Eagle);
   
  //  Calculate predPSS from aplha, beta and cumHistPSS for model section
  for (i in 1:n_yearPSS){
    predPSS[i] = alpha + beta * cumHistPSS[i];}
 
  // Sigma calculation for update
  //  Only years used in PF comparison are used here
    sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS[loc_pf_years_PSS]));

 
  // Sum curr_PSS
    cum_current_PSS = sum(curr_PSS);
 
  // Current years PSS prediction
    curr_predPSS = alpha + beta * cum_current_PSS;
 
 if(cum_current_Eagle > 0){ // Only use Eagle if Current years passage is observed
 
  // Eagle regression 
  for (i in 1:n_yearEagle){
 
 predEagle[i] = alpha_eagle + beta_eagle * cumHistEagle[i];}

  // Empirical SD for update  
 sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));
 
 //Eagle prediction 
 curr_predEagle = alpha_eagle + beta_eagle * cum_current_Eagle ;}
 
}

model {
  
  // PSS regression priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma ~ normal(0,100);
  
  // Eagle regression priors
  alpha_eagle ~ normal(0,1e15);
  beta_eagle ~ normal(0,1e15);
  sigma_eagle ~ normal(0,100);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // PSS regression likelihood 
  for(i in 1:n_yearPSS){
    
    log(totalEOS[i])~normal(log(predPSS[i]),sigma);
    
  }
  
  if(cum_current_Eagle > 0) // Only update likelihood when Eagle passage is observed in current year
  // Eagle regression likelihood
  for(i in 1:n_yearEagle){
    
    log(totalEOS[loc_eagle_years][i]) ~ normal( log(predEagle[i]), sigma_eagle);
    
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS),sigma_predPSS);
  
  
  if(cum_current_Eagle > 0){ // Only update likelihood when Eagle passage is observed in current year
  // Update likelihood for the posterior prediction
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
  
  // Variable for posterior abundance
  real RunSize;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(cum_current_Eagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
  //
  RunSize = exp(ln_RunSize);
  
}






