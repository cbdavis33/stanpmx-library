rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

# Leave-one-out cross-validation
filename_base <- "transit_fixed_ntr_1cmt_linear/Stan/Fits"

fit_ppa <- read_rds(file.path(filename_base, "transit_fixed_ntr_1cmt_ppa.rds"))
fit_ppa_no_delay <- read_rds(file.path(filename_base, "ppa_no_delay.rds"))
fit_ppa_zero_into_depot <- read_rds(file.path(filename_base, 
                                              "ppa_zero_into_depot.rds"))

fit_ppa_loo <- fit_ppa$loo(cores = min(4, parallel::detectCores()/2))
fit_ppa_no_delay_loo <- fit_ppa_no_delay$loo(cores = min(4, parallel::detectCores()/2))
fit_ppa_zero_into_depot_loo <- 
  fit_ppa_zero_into_depot$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_ppa_loo, fit_ppa_no_delay_loo, fit_ppa_zero_into_depot_loo)


# Leave-one-group-out cross-validation

stan_data_ppa <- jsonlite::read_json(
  "transit_fixed_ntr_1cmt_linear/Stan/Fits/Stan_Data/ppa.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data_ppa_no_delay <- jsonlite::read_json(
  "transit_fixed_ntr_1cmt_linear/Stan/Fits/Stan_Data/ppa_no_delay.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data_ppa_zero_into_depot <- jsonlite::read_json(
  "transit_fixed_ntr_1cmt_linear/Stan/Fits/Stan_Data/ppa_zero_into_depot.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data_ppa$n_sim <- stan_data_ppa_no_delay$n_sim <- stan_data_ppa_zero_into_depot$n_sim <- 1000

model_ppa <- cmdstan_model(
  "transit_fixed_ntr_1cmt_linear/Stan/Cross_Validation/transit_fixed_ntr_1cmt_ppa_logo_cv_sim.stan")

model_ppa_no_delay <- cmdstan_model(
  "depot_1cmt_linear/Stan/Cross_Validation/depot_1cmt_ppa_logo_cv_sim.stan")

model_ppa_zero_into_depot <- cmdstan_model(
  "zero_into_depot_1cmt_linear/Stan/Cross_Validation/zero_into_depot_1cmt_ppa_logo_cv_sim.stan")

logo_ppa <- model_ppa$generate_quantities(fit_ppa,
                                          data = stan_data_ppa,
                                          parallel_chains = 4,
                                          seed = 12345)

logo_ppa_zero_into_depot <- 
  model_ppa_zero_into_depot$generate_quantities(fit_ppa_zero_into_depot,
                                                data = stan_data_ppa_zero_into_depot,
                                                parallel_chains = 4,
                                                seed = 12345)

logo_ppa_no_delay <- model_ppa_no_delay$generate_quantities(fit_ppa_no_delay,
                                                            data = stan_data_ppa_no_delay,
                                                            parallel_chains = 4,
                                                            seed = 12345)

fit_ppa_logo <- loo(logo_ppa$draws("log_lik_subj"), 
                    cores = min(4, parallel::detectCores()/2))

fit_ppa_no_delay_logo <- loo(logo_ppa_no_delay$draws("log_lik_subj"), 
                             cores = min(4, parallel::detectCores()/2))

fit_ppa_zero_into_depot_logo <- loo(logo_ppa_zero_into_depot$draws("log_lik_subj"), 
                                    cores = min(4, parallel::detectCores()/2))

fit_ppa_logo
fit_ppa_no_delay_logo
fit_ppa_zero_into_depot_logo

loo_compare(fit_ppa_logo, fit_ppa_no_delay_logo, fit_ppa_zero_into_depot_logo)


