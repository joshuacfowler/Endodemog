## Grass endophyte population model
## Script for data clean up
setwd("~/Documents/R projects")
library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
LTREB_endodemog <- 
  read.csv("~/Documents/R projects/LTREBendodemog/endo_demog_long.csv")

View(LTREB_endodemog)
str(LTREB_endodemog)
dim(LTREB_endodemog)
