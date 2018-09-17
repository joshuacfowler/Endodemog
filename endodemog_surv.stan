
    data { 
    int<lower=0> N;                     // number of observations
    int<lower=0,upper=1> surv_t1[N]; // plant survival at time t+1 and target variable
    vector [N] logsize_t;            // log of plant size at time t
    }
    
    parameters {
    real alpha; // intercept
    real beta;  // slope
    }

    model {
    vector[N] mu;
    for(n in 1:N)
    mu = alpha + beta*logsize_t;  // linear predictor
  //Priors
    alpha ~ normal(0,100);
    beta ~ normal(0,100);
    
  //Likelihood
    surv_t1 ~ bernoulli_logit(mu);
    }
    
