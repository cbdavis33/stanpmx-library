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

nonmem_data <- read_csv("depot_1cmt_linear_t/Data/depot_1cmt_prop_t.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)


## Read in fit
fit <- read_rds("depot_1cmt_linear_t/Stan/Fits/depot_1cmt_prop_t.rds")

## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_p", "nu")

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
mcmc_dens(fit$draws(parameters_to_summarize %>% 
                      map(c("sigma", "nu"), str_subset, string = .) %>%
                      reduce(union)))
mcmc_dens(fit$draws(parameters_to_summarize %>% str_subset("cor")))

draws_df <- fit$draws(format = "draws_df")

# Individual point estimates (posterior mean)
est_ind <- draws_df %>%
  spread_draws(c(CL, VC, KA, eta_cl, eta_vc, eta_ka)[ID]) %>%
  mean_qi() %>% 
  select(ID, CL, VC, KA,
         eta_cl, eta_vc, eta_ka) %>% 
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


