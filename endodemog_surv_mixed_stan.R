## Grass endophyte population model with a bayesian framework
## Survival kernel with mixed effects


setwd("~/Documents/R projects")
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
str(POAL_data)
dim(POAL_data)

## create new column with log(size)
POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t)) 

str(POAL_data1)
dim(POAL_data1)
#View(POAL_data1)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

## here is the Stan model ##
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(120)

## MCMC settings
ni <- 500
nb <- 100
nc <- 1
## actual values from data
mean(surv_dat$surv_t1)


## data for GLMM
POAL_data1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t))

dim(POAL_data1)
surv_dat1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1))
size_dat1 <- POAL_data1 %>% 
  filter(!is.na(logsize_t))
year_dat <- POAL_data1 %>% 
  filter(!is.na(year_t))

POAL_data_list <- list(surv_t1 = surv_dat1$surv_t1, logsize_t = size_dat1$logsize_t, year = year_dat$year_t, N = nrow(POAL_data1))

str(POAL_data_list)

## GLMM for survival vs. log(size) with year effects

sink("endodemog_surv.stan")
cat("
    data { 
    int<lower=0> N;                   // number of observations
    int<lower=0,upper=1> surv_t1[N];  // plant survival at time t+1 and target variable
    int<lower=0> year_t[N];            // year of planting
    vector [N] logsize_t;             // log of plant size at time t
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
    ",fill=T)
sink()

stanmodel <- stanc("endodemog_surv.stan")

## Run the model by calling stan()
sm <- stan(file = "endodemog_surv.stan", data = POAL_data_list,
           iter = ni, warmup = nb, chains = nc)
print(sm)
## save the stanfit object so that it can be 
## called later without rerunning the model
saveRDS(sm, file = "endodemog_surv.rds")
sm <- readRDS(file = "endodemog_surv.rds")


## check convergence and posterior distributions
plot(sm)

traceplot(sm)

posterior <- as.matrix(sm)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("alpha", "beta"),
           prob = 0.8) + plot_title

shiny <- as.shinystan(sm)
launch_shinystan(shiny)

## I'm not sure how to pull out the coefficients from the stanfit object yet
x_dummy <- seq(min((POAL_data1$size_t),na.rm=T),max((POAL_data1$logsize_t),na.rm=T),0.1)

plot(POAL_data$surv_t1 ~ log(POAL_data$size_t), xlab = "Size in year t", ylab = "Survival in year t+1")        
points(surv_bin$mean_size, surv_bin$mean_surv, pch=16, cex=2)
lines(x_dummy,invlogit(coef(posterior)[1]+coef(posterior)[2]*x_dummy),col="red",lwd=3)

