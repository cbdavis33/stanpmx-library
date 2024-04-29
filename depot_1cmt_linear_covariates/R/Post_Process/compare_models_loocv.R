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

filename_base <- "depot_1cmt_linear_covariates/Stan/Fits"

fit_no_covariates <- read_rds(file.path(filename_base, 
                                        "depot_1cmt_prop_no_covariates.rds"))
fit_fixed_allometric <- read_rds(file.path(filename_base, 
                                           "depot_1cmt_prop_covariates_fixed_allometric.rds"))
fit_covariates <- read_rds(file.path(filename_base, 
                                     "depot_1cmt_prop_covariates.rds"))

fit_no_covariates_loo <- fit_no_covariates$loo(cores = min(4, parallel::detectCores()/2))
fit_fixed_allometric_loo <- fit_fixed_allometric$loo(cores = min(4, parallel::detectCores()/2))
fit_covariates_loo <- fit_covariates$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_no_covariates_loo, fit_fixed_allometric_loo, fit_covariates_loo)
