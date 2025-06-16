rm(list = ls())
cat("\014")

library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_1cmt_linear_ir1/Data/depot_1cmt_exp_ir1_exp.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID, cmt) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

fit <- read_rds("depot_1cmt_linear_ir1/Stan/Fits/depot_1cmt_exp_ir1_exp.rds")

## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka",
                                               "kin", "kout", "ic50")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_pk", "sigma_pd")

summary <- summarize_draws(fit$draws(parameters_to_summarize), 
                           mean, median, sd, mcse_mean,
                           ~quantile2(.x, probs = c(0.025, 0.975)), rhat,
                           ess_bulk, ess_tail)

summary %>%
  mutate(rse = sd/mean*100) %>% 
  select(variable, mean, sd, mcse_mean, q2.5, median, q97.5, rse, rhat, 
         starts_with("ess")) %>% 
  print(n = nrow(.))

mcmc_dens(fit$draws(parameters_to_summarize %>% str_subset("TV")))
mcmc_dens(fit$draws(parameters_to_summarize %>% str_subset("omega")))
mcmc_dens(fit$draws(parameters_to_summarize %>% str_subset("sigma")))
mcmc_dens(fit$draws(parameters_to_summarize %>% str_subset("cor")))

draws_df <- fit$draws(format = "draws_df")

# Individual point estimates (posterior mean)
est_ind <- draws_df %>%
  spread_draws(c(CL, VC, KA, KIN, KOUT, IC50, 
                 eta_cl, eta_vc, eta_ka, eta_kin, eta_kout, eta_ic50)[ID]) %>%
  mean_qi() %>% 
  select(ID, CL, VC, KA, KIN, KOUT, IC50, 
         eta_cl, eta_vc, eta_ka, eta_kin, eta_kout, eta_ic50) %>% 
  mutate(ID = factor(ID))

## Standardized Random Effects (posterior mean)
eta_std <- est_ind %>% 
  select(ID, starts_with("eta")) %>% 
  pivot_longer(cols = starts_with("eta"), 
               names_to = "parameter", 
               values_to = "eta",
               names_prefix = "eta_") %>% 
  group_by(parameter) %>% 
  mutate(eta_std = (eta - mean(eta))/sd(eta)) %>% 
  ungroup() %>% 
  mutate(parameter = toupper(parameter))

(est_ind %>% 
    ggplot(aes(x = eta_cl, y = eta_vc)) + 
    geom_point() +
    theme_bw() + 
    geom_smooth(method = "loess", span = 0.9) +
    ggtitle("Point Estimates")) +
  (draws_df %>%
     sample_draws(10) %>%
     spread_draws(`eta_.*`[ID], regex = TRUE) %>% 
     ungroup()  %>% 
     mutate(ID = factor(ID)) %>% 
     ggplot(aes(x = eta_cl, y = eta_vc)) + 
     geom_point() +
     theme_bw() + 
     geom_smooth(method = "loess", span = 0.9) +
     ggtitle("Samples"))

(est_ind %>% 
    ggplot(aes(x = eta_kin, y = eta_kout)) + 
    geom_point() +
    theme_bw() + 
    geom_smooth(method = "loess", span = 0.9) +
    ggtitle("Point Estimates")) +
  (draws_df %>%
     sample_draws(10) %>%
     spread_draws(`eta_.*`[ID], regex = TRUE) %>% 
     ungroup()  %>% 
     mutate(ID = factor(ID)) %>% 
     ggplot(aes(x = eta_kin, y = eta_kout)) + 
     geom_point() +
     theme_bw() + 
     geom_smooth(method = "loess", span = 0.9) +
     ggtitle("Samples"))
