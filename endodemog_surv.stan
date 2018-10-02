
    data { 
    int<lower=0> N;                   // number of observations
    int<lower=0> Y;                   // number of years
    int<lower=0,upper=1> surv_t1[N];  // plant survival at time t+1 and target variable (response)
    vector [N] logsize_t;             // log of plant size at time t (predictor (fixed effect))
    }
    
    parameters {
    real alpha;       // intercept
    real beta_size[Y];   // size parameter
    real beta_year[Y];   // year effects parameter
    }
    
    transformed parameters {
    real beta0_y[Y];
    for(y in 1:Y)
        beta0_y[y] = beta_size[y] + beta_year[y];
    }
    }
    model {

    // Linear Predictor
    vector[N] mu;
    for(n in 1:N)
       mu[n] = alpha + beta0_y[n]*logsize_t[n];  // linear predictor

    // Priors
    alpha ~ normal(0,100);
    beta_size ~ normal(0,100);
    beta_year ~ normal(0,10);
 
    // Likelihood
    surv_t1 ~ bernoulli_logit(mu);
    }
    
