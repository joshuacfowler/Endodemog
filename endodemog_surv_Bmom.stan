
    data { 
    int<lower=0> N; // number of observations
    int<lower=0,upper=1> surv_t1[N]; // plant survival at time t+1
    }
    
    parameters {
    real<lower=0, upper=1> mu;
    }
    
    model {
//Priors
    mu ~ beta(1,1);
//Likelihood
    for(n in 1:N)
      surv_t1[n] ~ bernoulli(mu);
  }
    
