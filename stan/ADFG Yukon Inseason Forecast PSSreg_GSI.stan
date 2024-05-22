// Yukon King Inseason Forecast version PSSreg_GSI

//Notes: 
//       

//Components:
//  1) Pre-season forecast  (Prior)
//  2) EOS-Can ~ PSS in normal space
//  3) GSI incorporated (Beta likelihood)
//  



data {
  
  // This years Preseason forecast
  real <lower=0>Pf;
  real <lower=0>Pf_sigma;
  int <lower=0> n_yearsPF;
  

  // EOS Canadian counts excluding myYear
  int <lower=0>n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int<lower=0> n_dayPSS;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years
  int<lower=0> n_yearPSS;
  int <lower=0>yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  matrix <lower=0>[n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year PSS passage up to myDay
  int <lower=0> n_curr_PSS;
  real <lower=0> curr_PSS[n_curr_PSS];
  
  // Mean historic GSI proportions for GSI Likelihood comparison
  real <lower=0> meanpropCAN[n_dayPSS];
  real <lower=0> sd_meanpropCAN[n_dayPSS];
  
  // Pointer vector for PF years in PSS
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  
}

transformed data{
  
    // Parameters for GSI beta likelihood
  real <lower=0> paramA[n_dayPSS];
  // 
  real <lower=0> paramB[n_dayPSS];
  
  
  
    // Calculate the A an B parameter for beta dist from observed mean & sd GSI
  for(x in 1:n_dayPSS){

      paramA[x] = ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;

      paramB[x] = paramA[x] * ((1/meanpropCAN[x])-1);}
  
  
  
}

// 
parameters {
   
  // Regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // GSI parameters
  real ln_phi;
  
  // Logit canadian proportion draws
  real propCAN_logit[n_dayPSS];
  
  // Log run size (Integrated run size estimate)
  real ln_RunSize;
  
}

transformed parameters{
  
  // Empty vector for cum sum of counts for each year
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS
  vector [n_yearPSS]predPSS;
  
  // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // This years prediction
  real curr_predPSS;
  
  // Empirical sd
  real sigma_predPSS;
  
  real  phi;
  
  // First days proportion
  vector<lower=0, upper=1> [n_dayPSS]propCAN;
   
  phi = exp(ln_phi);
 
  // Bring propCAN_logit into beta space
   for(t in 1:n_dayPSS){
     
      propCAN[t] = exp(propCAN_logit[t])/(1+exp(propCAN_logit[t]));
      }

 // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   
   cumHistPSS[i] = 0;
   
   for(j in 1:n_dayPSS){
     
   cumHistPSS[i] += (PSS_mat[j,i] * propCAN[j]);
 }
 }
 
 // Calculate predPSS from aplha, beta and cumHistPSS for model section
  for (i in 1:n_yearPSS){
    predPSS[i] = alpha + beta * cumHistPSS[i];}
 
 // Sigma calculation for update
 //  Only years with PF are used here
 sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS[loc_pf_years_PSS]));
 
 
 // Sum curr_PSS after ajusting by propCAN
 cum_current_PSS = 0;
 
 for(j in 1:n_dayPSS){
   
 cum_current_PSS += (curr_PSS[j]*propCAN[j]);

 }
 
 // Current years prediction
 curr_predPSS = alpha + beta * (cum_current_PSS);
 
}

model {
  
  // Regresion priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma ~ normal(0,10);
  
  // GSI random walk prior
  ln_phi ~ normal(-1,1);
  
  // First day GSI
  propCAN_logit[1] ~ normal(0,1);
  
  
  // GSI random walk 
  propCAN_logit[2:n_dayPSS] ~normal(propCAN_logit[1:n_dayPSS-1], phi);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  //Likelihood Function
  for(i in 1:n_yearPSS){
    
    log(totalEOS[i])~normal(log(predPSS[i]),sigma);
    
  }
  
  // Canadian Proportion likelihood
  propCAN ~ beta(paramA, paramB);
  
  // Update for the posterior integrated run size
  target += normal_lpdf(ln_RunSize | log(curr_predPSS),sigma_predPSS);

}

generated quantities{
  
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real RunSize;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  //
  RunSize = exp(ln_RunSize);
}




