## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates survival kernel written in STAN with mixed effects, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
LTREB_endodemog <- 
  read.csv("~/Documents/R projects/LTREBendodemog/endo_demog_long.csv")

# View(LTREB_endodemog)
str(LTREB_endodemog)
dim(LTREB_endodemog)

POAL_data <- LTREB_endodemog %>% 
  filter(species == "POAL")
# View(POAL_data)
str(POAL_data)
dim(POAL_data)

## create new column with log(size) and recode years and plot (including fixing values entered as R where they should be plot 16)
POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, 'R'=8, '3'=1, '4'=2, '8'=3, '9'=4, '10'=5, '11'=6, '15'=7, '16'=8, '17'=9, '19'=10, '151'=11, '152'=12, '153'=13, '154'=14, '155'=15, '156'=16, '157'=17, '158'=18)))
str(POAL_data1)
dim(POAL_data1)
# View(POAL_data1)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

## here is the Stan model ##
## run this to optimize computer system settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 200
nb <- 100
nc <- 1

## data for GLMM
POAL_data1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1)) 


dim(POAL_data1)
str(POAL_data1)
# View(POAL_data1)
surv_dat1 <- POAL_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
size_dat1 <- POAL_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
year_dat <- POAL_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
endo_dat <- POAL_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
origin_dat <- POAL_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1,
                                        origin == 16 ~ 1)))
plot_dat <- POAL_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))


dim(surv_dat1)
dim(size_dat1)
dim(year_dat)
dim(endo_dat)
dim(origin_dat)
dim(plot_dat)

POAL_data_list1 <- list(surv_t1 = surv_dat1$surv_t1, logsize_t = size_dat1$logsize_t, year_t = year_dat$year_t_index, N = 3241L, Y = 11L)

str(POAL_data_list1)

## GLMM for survival vs. log(size) with year effects

sink("endodemog_surv_year3.stan")
cat("
    
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> Y;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    }
    
    parameters {
    vector[2] beta;              // fixed intercept and slope
    vector[Y] beta_year;      // random year intercept
    }
    
    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] =  beta[1] + beta_year[year_t[n]] + beta[2]*logsize_t[n];
    }
    // Priors
    beta ~ normal(0,1e6);      // prior for slope
    beta_year ~ normal(0,1e6);   // prior for year random effects
 
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
    generated quantities{
    int yrep[N];
    vector[N] mu;
    
    // for posterior predictive check
    for (n in 1:N) {
      mu[n] =  beta[1] + beta_year[year_t[n]] + beta[2]*logsize_t[n];
     
      yrep[n] = bernoulli_logit_rng(mu[n]);
    }
    }
  
      ", fill = T)
sink()

stanmodel <- stanc("endodemog_surv_year3.stan")

## Run the model by calling stan()
sm <- stan(file = "endodemog_surv_year3.stan", data = POAL_data_list1,
           iter = ni, warmup = nb, chains = nc)

print(sm)
## save the stanfit object so that it can be 
## called later without rerunning the model
saveRDS(sm, file = "endodemog_surv_year3.rds")
sm <- readRDS(file = "endodemog_surv_year3.rds")


## check convergence and posterior distributions
surv_t1 <- as.vector(surv_dat1$surv_t1)
yrep <- as.matrix(sm, pars = "yrep")
mu <- as.matrix(sm, pars = "mu")
p <- invlogit(mu)

for(i in 1:100){
  y_resid <-  abs(surv_t1 - p);
  yrep_resid <- abs(yrep - p[i,]);
}

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))


plot(x = fit_yrep, y = fit, main = "Residuals full model")
abline(a = 0, b = 1, col = "red")




ppc_dens_overlay(surv_dat1$surv_t1, yrep[1:400, ])
ppc_stat(y = surv_dat1$surv_t1, yrep = yrep, stat = "mean")

traceplot(sm, pars = c("beta[1]", "beta[2]"))

posterior <- as.data.frame(sm)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("beta[1]", "beta[2]", "beta_year[1]", "beta_year[2]", "beta_year[3]", "beta_year[4]", "beta_year[5]"),
           prob = 0.8) + plot_title
sm_summary <- summary(sm)
inits <- get_inits(sm)
shiny <- as.shinystan(sm)
launch_shinystan(shiny)




print(sm)

plot(sm)
traceplot(sm)






# GLMM for Surv~ size +Endo + Origin  with year + plot random effects-------------------------
# Data list for this model
POAL_data_list <- list(surv_t1 = surv_dat1$surv_t1, 
                       logsize_t = size_dat1$logsize_t, 
                       endo = endo_dat$endo1, 
                       endo_index = endo_dat$endo_index,
                       origin = origin_dat$origin1,
                       year_t = year_dat$year_t_index, 
                       N = 3241L, K = 4L, nyear = 11L, nEndo = 2L)
str(POAL_data_list) 

sink("endodemog_surv_full_POAL8.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    int<lower=0, upper=1> endo[N];            // endophyte status 
    int<lower=1, upper=2> endo_index[N];       // index for endophyte effect
    int<lower=0, upper=1> origin[N];            // origin status
    }
    
    parameters {
    real alpha;                  // intercept
    vector[K] beta;              // predictor parameters
    matrix[nEndo, nyear] tau_year;      // random year effect
      
    real<lower=0> sigma_0[nEndo];        //year variance intercept
    }
    

    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
       + beta[4]*logsize_t[n]*endo[n] +
       + tau_year[endo_index[n], year_t[n]];
    }
    // Priors
    alpha ~ normal(0,1e6);      // prior for fixed intercept
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    for(n in 1:nyear){
    for(e in 1:nEndo){
    tau_year[e,n] ~ normal(0,sigma_0);   // prior for year random effects
    }}
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
   generated quantities{
    int yrep[N];
    vector[N] mu;
    
    // for posterior predictive check
    for(n in 1:N){
      mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
      +beta[4]*logsize_t[n]*endo[n]
      + tau_year[endo_index[n], year_t[n]];
      
      yrep[n] = bernoulli_logit_rng(mu[n]);
    }
    
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_surv_full_POAL8.stan")

## Run the model by calling stan()
sm <- stan(file = "endodemog_surv_full_POAL8.stan", data = POAL_data_list,
           iter = ni, warmup = nb, chains = nc)

print(sm, pars = c("beta[2]"))
## save the stanfit object so that it can be 
## called later without rerunning the model
saveRDS(sm, file = "endodemog_surv_full_POAL8.rds")
sm <- readRDS(file = "endodemog_surv_full_POAL8.rds")

print(sm)
plot(sm, pars = c("tau_year[1]", "sigma_0", "sigma_y", "sigma_e[1]"))

traceplot(sm, pars = c("beta[2]", "tau_year[1,1]"))
stan_dens(sm, pars = c("beta[2]"))


# checking residuals
surv_t1 <- as.vector(surv_dat1$surv_t1)
yrep <- as.matrix(sm, pars = "yrep")
mu <- as.matrix(sm, pars = "mu")
p <- invlogit(mu)

for(i in 1:100){
y_resid <-  abs(surv_t1 - p);
yrep_resid <- abs(yrep - p[i,]);
}

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))


plot(x = fit, y = fit_yrep, main = "surv Residuals full model")
abline(a = 0, b = 1, col = "red")

ppc_dens_overlay(surv_t1, yrep[1:100, ])
ppc_stat(y = surv_t1, yrep = yrep, stat = "mean")


posterior <-  as.data.frame(sm)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")

mcmc_areas(posterior, 
           pars = c("sigma_0[2]"),
           prob = 0.8) + plot_title

sm_summary <- summary(sm)
inits <- get_inits(sm)
shiny <- as.shinystan(sm)
launch_shinystan(shiny)
