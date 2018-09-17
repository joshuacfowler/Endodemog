## Grass endophyte population model with a bayesian framework
## Survival kernel
## Author: Joshua Fowler


setwd("~/Documents/R projects")
library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)

## Load data and remove NAs
LTREB_endodemog <- read.csv("~/Documents/R projects/LTREBendodemog/LTREB_endodemog.csv")
View(LTREB_endodemog)
POAL_data <- filter(LTREB_endodemog, species == "POAL")
View(POAL_data)
surv_dat <- POAL_data %>% 
  filter(!is.na(surv_t1))
size_dat <- POAL_data %>% 
  filter(!is.na(size_t))

dim(surv_dat)
dim(size_dat)
POAL_data1 <- POAL_data %>% 
  filter(!is.na(size_t)) %>%
  mutate(size_t, logsize_t = log(size_t)) %>% 
  filter(!is.na(surv_t1))
dim(POAL_data)
dim(POAL_data1)

## Create functions for linear predictors
invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

## Visualize survival as a function of log(size) ##

with(POAL_data,
     {plot(surv_t1 ~ log(size_t), xlab = "Volume in year t", ylab = "Survival in year t+1")}
)

surv_bin <- POAL_data %>% 
  mutate(size_bin = cut_interval(log(size_t), n=8)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_size = mean(log(size_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))

plot(POAL_data$surv_t1 ~ log(POAL_data$size_t), xlab = "Size in year t", ylab = "Survival in year t+1", col="gray")        
points(surv_bin$mean_size, surv_bin$mean_surv, pch=16, cex=2)

## Here is the basic R model for Survival ##

survival_model <- glm(surv_t1 ~ log(size_t), family = "binomial", data=POAL_data)
summary(survival_model)

## Below is the stan model ##
## Recommended setup options for MCMC simulations
rstan_options(auto_write = TRUE)              
options(mc.cores = parallel::detectCores())
set.seed(120)


## MCMC settings
ni <- 1000
nb <- 200
nc <- 3

## This first part is a model of the mean to learn Stan
## Actual mean value from survival data
mean(surv_dat$surv_t1)


## stan model of the mean of survival which has a Bernoulli distribution
sink("endodemog_surv_Bmom.stan")
cat("
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
    ",fill=T)
sink()

sm <- stanc("endodemog_surv_Bmom.stan")

## Call Stan from R and assign the output to "Bmom"
surv_data_list <- list(surv_t1 = surv_dat$surv_t1, N = nrow(surv_dat))
surv_data_list
Bmom <- stan(file = "LTREB_endodemog_Bmom.stan", data = surv_data_list,
            iter = ni, warmup = nb, chains = nc)


## Summarize posteriors
print(Bmom)
plot(Bmom)
pairs(Bmom)


## Prepare data for GLM

POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t))
surv_dat1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1))
size_dat1 <- POAL_data1 %>% 
  filter(!is.na(logsize_t))

POAL_data_list <- list(surv_t1 = surv_dat1$surv_t1, logsize_t = size_dat1$logsize_t, N = nrow(size_dat1))

str(POAL_data_list)


## GLM - Survival vs log(size)
sink("LTREB_endodemog.stan")
cat("
    data { 
    int<lower=0> N;                     // number of observations
    int<lower=0,upper=1> surv_t1[N]; // plant survival at time t+1 and target variable
    row_vector [N] logsize_t;            // log of plant size at time t
    }
    
    parameters {
    real alpha; // intercept
    real beta;  // slope
    }

    model {
    row_vector[N] mu;
    for(n in 1:N)
    mu = alpha + beta*logsize_t;  // linear predictor
  //Priors
    alpha ~ normal(0,100);
    beta ~ normal(0,100);
    
  //Likelihood
    surv_t1 ~ bernoulli_logit(mu);
    }
    ",fill=T)
sink()

sm <- stanc("LTREB_endodemog.stan")

## Call Stan from R and assign the output to "out"
out <- stan(file = "LTREB_endodemog.stan", data = POAL_data_list,
            iter = ni, warmup = nb, chains = nc)

## Summarize posteriors
print(out)
plot(out)
pairs(out)
