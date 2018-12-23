## Grass endophyte population model with a bayesian framework
## Survival kernel
## Tom makes a commit.

setwd("~/Documents/R projects")
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
str(POAL_data)
dim(POAL_data)
#View(POAL_data)

## create new column with log(size)
POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t)) 

str(POAL_data1)
dim(POAL_data1)

surv_dat1 <- c(na.omit(POAL_data1$surv_t1))
size_dat1 <- c(na.omit(POAL_data1$logsize_t))
year_dat1 <- c(na.omit(POAL_data1$year_t))

# surv_dat1 <- POAL_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
# size_dat1 <- POAL_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))

dim(POAL_data1)
dim(surv_dat1)
dim(size_dat1)
#View(POAL_data1)

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

## Here is the basic model for Survival ##

survival_model <- glm(surv_t1 ~ log(size_t), family = "binomial", data=POAL_data)
summary(survival_model)

## here is the Stan model ##
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(120)


## MCMC settings
ni <- 100
nb <- 10
nc <- 1
## actual values from data
mean(surv_dat1$surv_t1)
mean(surv_dat1)


##Bernoulli model of the mean survival
sink("endodemog_surv_Bmom.stan")
cat("
    data { 
    int<lower=0> N; // number of observations
    int surv_t1[N]; // plant survival at time t+1
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



stanmodel <- stanc("endodemog_surv_Bmom.stan")
## Call Stan from R
surv_data_list <- list(surv_t1 = surv_dat1, N = 3241)
str(surv_data_list)
Bmom <- stan(file = "endodemog_surv_Bmom.stan", data = surv_data_list,
            iter = ni, warmup = nb, chains = nc)


## Summarize parameters
print(Bmom)
plot(Bmom)

## Traceplot for parameters
traceplot(Bmom)


## data for GLM


POAL_data_list <- list(surv_t1 = surv_dat1, logsize_t = size_dat1, N = as.integer(3241), K = as.integer(1))

str(POAL_data_list)
#View(POAL_data_list)

## GLM - Survival vs log(size)
sink("endodemog_surv.stan")
cat("
    data { 
    int<lower=0> N;                   // number of observations
    //int<lower=0> K;    // number of predictors
    int surv_t1[N];  // plant survival at time t+1 and target variable
    real logsize_t[N];        // log of plant size at time t
    }
    
    parameters {
    real alpha; //fixed intercept
    real beta; // fixed slope 
    }

    model {
    vector[N] mu;
  //Priors
    alpha ~ normal(0,10);
    beta ~ normal(0,10);
  //Likelihood
  for(n in 1:N){
    mu[n] = alpha + beta*logsize_t[n]; // linear predictor
    surv_t1[n] ~ bernoulli_logit(mu); // likelihood
    }
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

