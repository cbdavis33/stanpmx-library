rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

# Leave-one-out cross-validation
filename_base <- "zero_into_depot_2cmt_linear/Stan/Fits"

fit_ppa <- read_rds(file.path(filename_base, "zero_into_depot_2cmt_ppa.rds"))
fit_prop <- read_rds(file.path(filename_base, "ppa_fit_with_prop.rds"))
fit_exp <- read_rds(file.path(filename_base, "ppa_fit_with_exp.rds"))
fit_ppa_no_delay <- read_rds(file.path(filename_base, "ppa_no_delay.rds"))

fit_ppa_loo <- fit_ppa$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_loo <- fit_prop$loo(cores = min(4, parallel::detectCores()/2))
fit_exp_loo <- fit_exp$loo(cores = min(4, parallel::detectCores()/2))
fit_ppa_no_delay_loo <- fit_ppa_no_delay$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_ppa_loo, fit_prop_loo, fit_exp_loo, fit_ppa_no_delay_loo)



# Leave-one-group-out cross-validation

stan_data <- jsonlite::read_json(
  "zero_into_depot_2cmt_linear/Stan/Fits/Stan_Data/ppa.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 1000

model_ppa <- cmdstan_model(
  "zero_into_depot_2cmt_linear/Stan/Cross_Validation/zero_into_depot_2cmt_ppa_logo_cv_sim.stan")

model_prop <- cmdstan_model(
  "zero_into_depot_2cmt_linear/Stan/Cross_Validation/zero_into_depot_2cmt_prop_logo_cv_sim.stan")

model_exp <- cmdstan_model(
  "zero_into_depot_2cmt_linear/Stan/Cross_Validation/zero_into_depot_2cmt_exp_logo_cv_sim.stan")

model_ppa_no_delay <- cmdstan_model(
  "depot_2cmt_linear/Stan/Cross_Validation/depot_2cmt_ppa_logo_cv_sim.stan")

logo_ppa <- model_ppa$generate_quantities(fit_ppa,
                                          data = stan_data,
                                          parallel_chains = 4,
                                          seed = 12345)

logo_prop <- model_prop$generate_quantities(fit_prop,
                                            data = stan_data,
                                            parallel_chains = 4,
                                            seed = 12345)

logo_exp <- model_exp$generate_quantities(fit_exp,
                                          data = stan_data,
                                          parallel_chains = 4,
                                          seed = 12345)

stan_data$rate <- rep(0, times = stan_data$n_total)

logo_ppa_no_delay <- model_ppa_no_delay$generate_quantities(fit_ppa_no_delay,
                                                            data = stan_data,
                                                            parallel_chains = 4,
                                                            seed = 12345)

fit_ppa_logo <- loo(logo_ppa$draws("log_lik_subj"), 
                    cores = min(4, parallel::detectCores()/2))

fit_prop_logo <- loo(logo_prop$draws("log_lik_subj"), 
                     cores = min(4, parallel::detectCores()/2))

fit_exp_logo <- loo(logo_exp$draws("log_lik_subj"), 
                    cores = min(4, parallel::detectCores()/2))

fit_ppa_no_delay_logo <- loo(logo_ppa_no_delay$draws("log_lik_subj"), 
                             cores = min(4, parallel::detectCores()/2))

fit_ppa_logo
fit_prop_logo
fit_exp_logo
fit_ppa_no_delay_logo

loo_compare(fit_ppa_logo, fit_prop_logo, fit_exp_logo, fit_ppa_no_delay_logo)


