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

nonmem_data <- read_csv(
  "zero_into_depot_2cmt_linear/Data/zero_into_depot_2cmt_exp.csv",
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
fit <- read_rds(
  "zero_into_depot_2cmt_linear/Stan/Fits/zero_into_depot_2cmt_exp.rds")


## Most of the below can be done interactively with shinystan
# as.shinystan(fit) %>% 
#   launch_shinystan()


## Overall summary
# https://mc-stan.org/docs/cmdstan-guide/diagnose_utility.html#building-the-diagnose-command
# https://mc-stan.org/misc/warnings.html
# treedepth, divergences, E-BFMI, ESS, R-hat
fit$cmdstan_diagnose()

## Divergences, max_treedepth, E-BFMI values in a list
fit$diagnostic_summary()

# fileConn <- file("zero_into_depot_2cmt_linear/output.txt")
# writeLines(fit$cmdstan_summary()$stdout, fileConn)
# close(fileConn)

## Plots of divergences, treedepth, and E-BFMI
## See https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html
##   for a more in-depth examination of divergences
mcmc_nuts_divergence(nuts_params(fit), log_posterior(fit))
mcmc_nuts_treedepth(nuts_params(fit), log_posterior(fit))
mcmc_nuts_energy(nuts_params(fit), log_posterior(fit))

# Plots
parameters_to_examine <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                           str_c("omega_", c("cl", "vc", "ka", "q", "vp", "dur")),
                           str_subset(fit$metadata()$stan_variables, "cor_"),
                           "sigma")

## R-hat
rhats <- bayesplot::rhat(fit, pars = parameters_to_examine) %>% 
  enframe(name = "variable", value = "r_hat") 

rhats %>% 
  arrange(desc(r_hat))

rhats %>% 
  deframe() %>% 
  mcmc_rhat() +
  yaxis_text(hjust = 1)


## Effective sample size

# neffs <- bayesplot::neff_ratio(fit, pars = parameters_to_examine) %>% 
#   enframe(name = "variable", value = "neff_ratio")

neffs <- summarize_draws(fit$draws(parameters_to_examine), 
                         ess_bulk, ess_tail) 

# Want all effective sample sizes > 400
neffs
neffs %>% 
  summarize(across(starts_with("ess"), min))

neff_ratios <- neffs %>% 
  mutate(across(starts_with("ess"), 
                \(x) x/(fit$metadata()$iter_sampling*fit$num_chains())))

(neff_ratios %>% 
    select(variable, ess_bulk) %>% 
    deframe() %>% 
    mcmc_neff() +
    yaxis_text(hjust = 1) +
    ggtitle("Bulk")) +
  (neff_ratios %>% 
     select(variable, ess_tail) %>% 
     deframe() %>% 
     mcmc_neff() +
     yaxis_text(hjust = 1) +
     ggtitle("Tail")) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Autocorrelation (not hugely important once you've already looked at N_eff)
mcmc_acf(fit$draws(parameters_to_examine))


## Other MCMC diagnostics
# dens_overlay shows marginal posterior distributions by chain. We want to see 
#   overlapping densities
# traceplots show how chains mix over time. We want to see a "fuzzy caterpillar"
# In the violin plot, we want to see similar violins across chains with matching
#   quantile lines
mcmc_combo(fit$draws(parameters_to_examine %>% str_subset("TV")),
           combo = c("dens_overlay", "trace", "violin"))
mcmc_combo(fit$draws(parameters_to_examine %>% str_subset("omega")),
           combo = c("dens_overlay", "trace", "violin"))
mcmc_combo(fit$draws(parameters_to_examine %>% str_subset("sigma")),
           combo = c("dens_overlay", "trace", "violin"))
mcmc_combo(fit$draws(parameters_to_examine %>% str_subset("cor")),
           combo = c("dens_overlay", "trace", "violin"))

# rank_hist plots look how MCMC samples mix between chains. We want to see 
# uniform distributions
mcmc_rank_hist(fit$draws(parameters_to_examine %>% str_subset("TV")))
mcmc_rank_hist(fit$draws(parameters_to_examine %>% str_subset("omega")))
mcmc_rank_hist(fit$draws(parameters_to_examine %>% str_subset("sigma")))
mcmc_rank_hist(fit$draws(parameters_to_examine %>% str_subset("cor")))

# Alternatively, visualize the ranks from mcmc_rank_hist() with overlaid lines.
# We want to see the lines mixed together with none systematically above or
# below others
mcmc_rank_overlay(fit$draws(parameters_to_examine %>% str_subset("TV")))
mcmc_rank_overlay(fit$draws(parameters_to_examine %>% str_subset("omega")))
mcmc_rank_overlay(fit$draws(parameters_to_examine %>% str_subset("sigma")))
mcmc_rank_overlay(fit$draws(parameters_to_examine %>% str_subset("cor")))


## rank_ecdf plots similarly look how MCMC samples mix between chains. We want 
## lines to be within the confidence bands. I find it easier to see with 
## `plot_diff = TRUE`

# mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("TV")))
# mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("omega")))
# mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("sigma")))
# mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("cor")))

mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("TV")),
               plot_diff = TRUE)
mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("omega")),
               plot_diff = TRUE)
mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("sigma")),
               plot_diff = TRUE)
mcmc_rank_ecdf(fit$draws(parameters_to_examine %>% str_subset("cor")),
               plot_diff = TRUE)


bayesplot::color_scheme_set('viridis')
mcmc_pairs(fit$draws(map(c("TV", "omega"), str_subset, 
                         string = parameters_to_examine) %>% 
                       reduce(union)),
           diag_fun = "dens",
           off_diag_fun = "hex")

mcmc_pairs(fit$draws(c("omega_cl", str_c("CL[", 1:9, "]"))),
           transform = list(omega_cl = "log"),
           diag_fun = "dens",
           off_diag_fun = "hex")


mcmc_pairs(fit$draws(c("omega_cl", "CL[1]", "omega_vc", "VC[1]",
                       "omega_ka", "KA[1]", "omega_dur", "DUR[1]")),
           transform = list(omega_cl = "log",
                            omega_vc = "log",
                            omega_ka = "log",
                            omega_dur = "log"),
           diag_fun = "dens",
           off_diag_fun = "hex")

