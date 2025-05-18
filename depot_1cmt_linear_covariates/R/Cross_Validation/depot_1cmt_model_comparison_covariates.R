rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

# Leave-one-out cross-validation
filename_base <- "depot_1cmt_linear_covariates/Stan/Fits"

fit_prop <- read_rds(file.path(filename_base, "depot_1cmt_prop.rds"))
fit_prop_covariates <- read_rds(file.path(filename_base, 
                                          "depot_1cmt_prop_covariates.rds"))
fit_prop_covariates_fixed_allometric <- 
  read_rds(file.path(filename_base, "depot_1cmt_prop_covariates_fixed_allometric.rds"))

fit_prop_loo <- fit_prop$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_covariates_loo <- 
  fit_prop_covariates$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_covariates_fixed_allometric_loo <- 
  fit_prop_covariates_fixed_allometric$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_prop_loo, fit_prop_covariates_loo,
            fit_prop_covariates_fixed_allometric_loo)


# Leave-one-group-out cross-validation

stan_data <- jsonlite::read_json(
  "depot_1cmt_linear_covariates/Stan/Fits/Stan_Data/prop.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data_covariates <- jsonlite::read_json(
  "depot_1cmt_linear_covariates/Stan/Fits/Stan_Data/prop_covariates.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data_covariates_fixed_allometric <- jsonlite::read_json(
  "depot_1cmt_linear_covariates/Stan/Fits/Stan_Data/prop_covariates_fixed_allometric.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 100
stan_data_covariates$n_sim <- 100
stan_data_covariates_fixed_allometric$n_sim <- 100


model_prop <- cmdstan_model(
  "depot_1cmt_linear/Stan/Cross_Validation/depot_1cmt_prop_logo_cv_sim.stan")

model_prop_covariates <- cmdstan_model(
  "depot_1cmt_linear_covariates/Stan/Cross_Validation/depot_1cmt_prop_logo_cv_sim_covariates.stan")

model_prop_covariates_fixed_allometric <- cmdstan_model(
  "depot_1cmt_linear_covariates/Stan/Cross_Validation/depot_1cmt_prop_logo_cv_sim_covariates_fixed_allometric.stan")

logo_prop <- model_prop$generate_quantities(fit_prop,
                                            data = stan_data,
                                            parallel_chains = 4,
                                            seed = 12345)

logo_prop_covariates <- 
  model_prop_covariates$generate_quantities(fit_prop_covariates,
                                            data = stan_data_covariates,
                                            parallel_chains = 4,
                                            seed = 12345)

logo_prop_covariates_fixed_allometric <- 
  model_prop_covariates_fixed_allometric$generate_quantities(
    fit_prop_covariates_fixed_allometric,
    data = stan_data_covariates_fixed_allometric,
    parallel_chains = 4,
    seed = 12345)

fit_prop_logo <- loo(logo_prop$draws("log_lik_subj"), 
                     cores = min(4, parallel::detectCores()/2))

fit_prop_logo_covariates <- loo(logo_prop_covariates$draws("log_lik_subj"), 
                                cores = min(4, parallel::detectCores()/2))

fit_prop_logo_covariates_fixed_allometric <- 
  loo(logo_prop_covariates_fixed_allometric$draws("log_lik_subj"), 
      cores = min(4, parallel::detectCores()/2))

loo_compare(fit_prop_logo, fit_prop_logo_covariates, 
            fit_prop_logo_covariates_fixed_allometric)


