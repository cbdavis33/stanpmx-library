rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_prop_covariates.rds")

stan_data <- jsonlite::read_json(
  "iv_2cmt_linear_covariates/Stan/Fits/Stan_Data/prop_covariates.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 100

model <- cmdstan_model(
  "iv_2cmt_linear_covariates/Stan/Cross_Validation/iv_2cmt_prop_logo_cv_sim_covariates.stan")

logo_cv <- model$generate_quantities(fit,
                                     data = stan_data,
                                     parallel_chains = 4,
                                     seed = 12345)

logo <- loo(logo_cv$draws("log_lik_subj"), 
            cores = min(4, parallel::detectCores()/2))

plot(logo, label_points = TRUE)
