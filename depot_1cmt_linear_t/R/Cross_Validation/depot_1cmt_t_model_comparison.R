rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

# Leave-one-out cross-validation
filename_base <- "depot_1cmt_linear_t/Stan/Fits"

fit_ppa <- read_rds(file.path(filename_base, "depot_1cmt_ppa_t.rds"))
fit_prop <- read_rds(file.path(filename_base, "ppa_fit_with_prop.rds"))
fit_ppa_norm <- read_rds(file.path(filename_base, "ppa_fit_with_ppa_norm.rds"))
fit_prop_norm <- read_rds(file.path(filename_base, "ppa_fit_with_prop_norm.rds"))

fit_ppa_loo <- fit_ppa$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_loo <- fit_prop$loo(cores = min(4, parallel::detectCores()/2))
fit_ppa_norm_loo <- fit_ppa_norm$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_norm_loo <- fit_prop_norm$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_ppa_loo, fit_prop_loo, fit_ppa_norm_loo, fit_prop_norm_loo)



# Leave-one-group-out cross-validation

stan_data <- jsonlite::read_json(
  "depot_1cmt_linear_t/Stan/Fits/Stan_Data/ppa.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 200

model_ppa <- cmdstan_model(
  "depot_1cmt_linear_t/Stan/Cross_Validation/depot_1cmt_ppa_t_logo_cv_sim.stan")

model_prop <- cmdstan_model(
  "depot_1cmt_linear_t/Stan/Cross_Validation/depot_1cmt_prop_t_logo_cv_sim.stan")

model_ppa_norm <- cmdstan_model(
  "depot_1cmt_linear/Stan/Cross_Validation/depot_1cmt_ppa_logo_cv_sim.stan")

model_prop_norm <- cmdstan_model(
  "depot_1cmt_linear/Stan/Cross_Validation/depot_1cmt_prop_logo_cv_sim.stan")

logo_ppa <- model_ppa$generate_quantities(fit_ppa,
                                          data = stan_data,
                                          parallel_chains = 4,
                                          seed = 12345)

logo_prop <- model_prop$generate_quantities(fit_prop,
                                            data = stan_data,
                                            parallel_chains = 4,
                                            seed = 12345)

logo_ppa_norm <- model_ppa_norm$generate_quantities(fit_ppa_norm,
                                                    data = stan_data,
                                                    parallel_chains = 4,
                                                    seed = 12345)

logo_prop_norm <- model_prop_norm$generate_quantities(fit_prop_norm,
                                                      data = stan_data,
                                                      parallel_chains = 4,
                                                      seed = 12345)

fit_ppa_logo <- loo(logo_ppa$draws("log_lik_subj"), 
                    cores = min(4, parallel::detectCores()/2))

fit_prop_logo <- loo(logo_prop$draws("log_lik_subj"), 
                     cores = min(4, parallel::detectCores()/2))

fit_ppa_norm_logo <- loo(logo_ppa_norm$draws("log_lik_subj"), 
                         cores = min(4, parallel::detectCores()/2))

fit_prop_norm_logo <- loo(logo_prop_norm$draws("log_lik_subj"), 
                          cores = min(4, parallel::detectCores()/2))

fit_ppa_logo
fit_prop_logo
fit_ppa_norm_logo
fit_prop_norm_logo

loo_compare(fit_ppa_logo, fit_prop_logo, fit_ppa_norm_logo, fit_prop_norm_logo)


