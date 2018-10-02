
    data { 
    int<lower=1> N;                   // number of observations
    int<lower=1> Y;                   // number of years
    int<lower=0> year_t[N];            // year of recruitment
    int <lower=0,upper=1> surv_t1[N];  // plant survival at time t+1 and target variable (response)
    real logsize_t[N];             // log of plant size at time t (predictor)
    }
    
    parameters {
    real alpha0;              // fixed intercept
    real beta0;               // fixed slope
    vector[Y] beta_year;      // random year intercept
    }
    
    model {
    vector[N] mu;
    // Priors
    beta_year ~ normal(0,100); // year random effects
    // Linear Predictor
    for(n in 1:N)
       mu[n] = alpha0 + beta_year[year_t[n]] + beta0 * logsize_t[n];
    // Likelihood
    surv_t1 ~ bernoulli_logit(mu);
    }
      
