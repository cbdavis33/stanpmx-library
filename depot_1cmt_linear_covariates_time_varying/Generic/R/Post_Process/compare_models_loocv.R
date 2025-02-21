rm(list = ls())
cat("\014")

library(gganimate)
library(loo)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(posterior)
library(tidyverse)

filename_base <- "depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits"

fit_fixed_allometric <- read_rds(file.path(filename_base, 
                                           "prop_fixed_allometric_locf.rds"))
fit_covariates <- read_rds(file.path(filename_base, 
                                     "prop_locf.rds"))

fit_fixed_allometric_loo <- fit_fixed_allometric$loo(cores = min(4, parallel::detectCores()/2))
fit_covariates_loo <- fit_covariates$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_fixed_allometric_loo, fit_covariates_loo)
