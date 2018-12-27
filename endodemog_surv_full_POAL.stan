
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> nplot;                     // number of plots (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> plot[N];                     // plot id
    
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    int<lower=0, upper = 1> endo[N];            // endophyte status 
    int<lower=0, upper=1> origin[N];            // origin status
    }
    
    parameters {
    real alpha;                  // intercept
    vector[K] beta;              // predictor parameters
    vector[nyear] beta_year;      // random year effect
    vector[nplot] beta_plot;      // random plot effect
    
    real<lower=0> sigma_y;        //year variance intercept
    real sigma_e;                //effect of endo on variance
    }
    
    transformed parameters {
    real<lower=0> sigma_0;       // year variance
    for(n in 1:N){
    sigma_0 = sigma_y + sigma_e*endo[n];
    }
    }
    
    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] = alpha + beta[1]*logsize_t[n] + beta[1]*endo[n] + beta[3]*origin[n] 
       + beta[4]*logsize_t[n]*endo[n]
       + beta_year[year_t[n]] + beta_plot[plot[n]];
    }
    // Priors
    alpha ~ normal(0,1e6);      // prior for fixed intercept
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    beta_year ~ normal(0,sigma_0);   // prior for year random effects
    beta_plot ~ normal(0, 1e6); // prior for plot random effects
 
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
  
      
