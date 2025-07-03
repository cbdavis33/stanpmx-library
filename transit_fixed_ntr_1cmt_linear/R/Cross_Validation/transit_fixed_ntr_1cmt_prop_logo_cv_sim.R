rm(list = ls())
cat("\014")

library(loo)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds(
  "transit_fixed_ntr_1cmt_linear/Stan/Fits/transit_fixed_ntr_1cmt_prop.rds")

stan_data <- jsonlite::read_json(
  "transit_fixed_ntr_1cmt_linear/Stan/Fits/Stan_Data/prop.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$n_sim <- 1000

model <- cmdstan_model(
  "transit_fixed_ntr_1cmt_linear/Stan/Cross_Validation/transit_fixed_ntr_1cmt_prop_logo_cv_sim.stan")

logo_cv <- model$generate_quantities(fit,
                                     data = stan_data,
                                     parallel_chains = 4,
                                     seed = 12345)

logo <- loo(logo_cv$draws("log_lik_subj"), 
            cores = min(4, parallel::detectCores()/2))
logo
plot(logo, label_points = TRUE)
