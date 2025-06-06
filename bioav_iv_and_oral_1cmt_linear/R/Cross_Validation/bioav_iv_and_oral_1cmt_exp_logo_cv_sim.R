rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("bioav_iv_and_oral_1cmt_linear/Stan/Fits/bioav_iv_and_oral_1cmt_exp.rds")

stan_data <- jsonlite::read_json(
  "bioav_iv_and_oral_1cmt_linear/Stan/Fits/Stan_Data/exp.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 1000

model <- cmdstan_model(
  "bioav_iv_and_oral_1cmt_linear/Stan/Cross_Validation/bioav_iv_and_oral_1cmt_exp_logo_cv_sim.stan")

logo_cv <- model$generate_quantities(fit,
                                     data = stan_data,
                                     parallel_chains = 4,
                                     seed = 12345)

logo <- loo(logo_cv$draws("log_lik_subj"), 
            cores = min(4, parallel::detectCores()/2))
logo
plot(logo, label_points = TRUE)
