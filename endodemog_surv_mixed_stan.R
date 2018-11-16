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

#View(LTREB_endodemog)
str(LTREB_endodemog)
dim(LTREB_endodemog)

POAL_data <- LTREB_endodemog %>% 
  filter(species == "POAL")
# View(POAL_data)
str(POAL_data)
dim(POAL_data)

## create new column with log(size) and recode years
POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11)))

str(POAL_data1)
dim(POAL_data1)
# View(POAL_data1)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

## here is the Stan model ##
## run this to optimize computer system settings
Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 100
nb <- 10
nc <- 1

## data for GLMM
POAL_data1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(year_t_index))
  

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
  filter(!is.na(year_t_index))
endo_dat <- POAL_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(endo)-4) ## recoding endo to 0 or 1



dim(surv_dat1)
dim(size_dat1)
dim(year_dat)
dim(endo_dat)

POAL_data_list <- list(surv_t1 = surv_dat1$surv_t1, logsize_t = size_dat1$logsize_t, endo = endo_dat$endo1, year_t = as.integer(year_dat$year_t_index), N = nrow(POAL_data1), Y = nlevels(year_dat$year_t_index))

str(POAL_data_list)

## GLMM for survival vs. log(size) with year effects

sink("endodemog_surv_year.stan")
cat("
    data { 
    int<lower=0> N;                     // number of observations
    int<lower=0> Y;                     // number of years (used as index)
    int<lower=0> year_t[N];             // year of recruitment
    int<lower=0,upper=1> surv_t1[N];   // plant survival at time t+1 and target variable (response)
    real<lower=0> logsize_t[N];                  // log of plant size at time t (predictor)
    }
    
    parameters {
    real alpha0;              // fixed intercept
    real beta0;               // fixed slope
    vector[Y] beta_year;      // random year intercept
    }
    
    model {
    real mu[N];
    // Priors
    alpha0 ~ normal(0,100);     // prior for intercept
    beta0 ~ normal(0,100);      // prior for slope
    beta_year ~ normal(0,100);   // prior for year random effects
    // Linear Predictor
    for(n in 1:N){
       mu[n] = alpha0 + beta_year[year_t[n]] + beta0*logsize_t[n];
    // Likelihood
      surv_t1[n] ~ bernoulli_logit(mu[n]);
      }
    }
      ", fill = T)
sink()

stanmodel <- stanc("endodemog_surv_year.stan")

## Run the model by calling stan()
sm <- stan(file = "endodemog_surv_year.stan", data = POAL_data_list,
           iter = ni, warmup = nb, chains = nc)

print(sm)
## save the stanfit object so that it can be 
## called later without rerunning the model
saveRDS(sm, file = "endodemog_surv_year.rds")
sm <- readRDS(file = "endodemog_surv_year.rds")


## check convergence and posterior distributions
plot(sm)

traceplot(sm)

posterior <- as.matrix(sm)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("alpha0", "beta0", "beta_year[1]"),
           prob = 0.8) + plot_title

shiny <- as.shinystan(sm)
launch_shinystan(shiny)


## GLMM for survival vs. log(size) with year and Endo effects

sink("endodemog_surv_y_e.stan")
cat("
    data { 
    int<lower=0> N;                     // number of observations
    int<lower=0> Y;                     // number of years (used as index)
    //int<lower=0> year_t[N];             // year of recruitment
    real endo[N];          // status as (E+ or E-)
    int<lower=0,upper=1> surv_t1[N];   // plant survival at time t+1 and target variable (response)
    real<lower=0> logsize_t[N];                  // log of plant size at time t (predictor)
    }
    
    parameters {
    real alpha0;              // fixed intercept
    real beta0;               // fixed slope
    real beta_endo;           // endo slope
   // vector[Y] beta_year;      // random year intercept
 
    }
    
    model {
    vector[N] mu;
    // Priors
    alpha0 ~ normal(0,100);        // prior for intercept
    beta0 ~ normal(0,100);         // prior for slope
    beta_endo ~ beta(1,1);         // prior for endo
   // beta_year ~ normal(0,100);   // prior for year random effects
    // Linear Predictor
    for(n in 1:N)
    mu[n] = alpha0 + beta0*logsize_t[n] + beta_endo*endo[n];
    // Likelihood
    for(n in 1:N)
    surv_t1[n] ~ bernoulli_logit(mu);
    
    }
    

    ",fill=T)
sink()

stanmodel <- stanc("endodemog_surv_y_e.stan")

## Run the model by calling stan()
sm <- stan(file = "endodemog_surv_y_e.stan", data = POAL_data_list,
           iter = ni, warmup = nb, chains = nc)
saveRDS(sm, file = "endodemog_surv_y_e.rds")
sm <- readRDS(file = "endodemog_surv_y_e.rds")


print(sm)

plot(sm)
traceplot(sm)




## I'm not sure how to pull out the coefficients from the stanfit object yet
x_dummy <- seq(min((POAL_data1$size_t),na.rm=T),max((POAL_data1$logsize_t),na.rm=T),0.1)

plot(POAL_data$surv_t1 ~ log(POAL_data$size_t), xlab = "Size in year t", ylab = "Survival in year t+1")        
points(surv_bin$mean_size, surv_bin$mean_surv, pch=16, cex=2)
lines(x_dummy,invlogit(coef(posterior)[1]+coef(posterior)[2]*x_dummy),col="red",lwd=3)

