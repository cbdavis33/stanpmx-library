rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "iv_2cmt_linear_covariates/Stan/Fits/Stan_Data/ppa_covariates.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model(
  "iv_2cmt_linear_covariates/Stan/Fit/iv_2cmt_ppa_covariates.stan",
  cpp_options = list(stan_threads = TRUE))

priors <- model$sample(data = stan_data,
                       seed = 235813,
                       chains = 4,
                       parallel_chains = 4,
                       threads_per_chain = 1,
                       iter_warmup = 500,
                       iter_sampling = 1000,
                       adapt_delta = 0.8,
                       refresh = 500,
                       max_treedepth = 10,
                       init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                                              TVVC = rlnorm(1, log(3), 0.3),
                                              TVQ = rlnorm(1, log(1), 0.3),
                                              TVVP = rlnorm(1, log(4), 0.3),
                                              theta_cl_wt = rnorm(1), 
                                              theta_vc_wt = rnorm(1), 
                                              theta_q_wt = rnorm(1), 
                                              theta_vp_wt = rnorm(1), 
                                              theta_vc_race_asian = rnorm(1), 
                                              theta_cl_egfr = rnorm(1),
                                              omega = rlnorm(4, log(0.3), 0.3),
                                              sigma_p = rlnorm(1, log(0.2), 0.3)))


fit <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_ppa_covariates.rds")

draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_subset(fit$metadata()$stan_variables, "theta"),
                             str_c("omega_", c("cl", "vc", "q", "vp")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_p", "sigma_a")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma_a, names_to = "variable", values_to = "value")

(target_comparison_tv <- draws_all_df %>% 
    filter(str_detect(variable, "TV")) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("CL", "VC", "Q", "VP")))) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() +
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1))

(target_comparison_covariates <- draws_all_df %>% 
    filter(str_detect(variable, "theta_")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("theta_cl_wt", "theta_vc_wt",
                                        "theta_q_wt", "theta_vp_wt",
                                        "theta_vc_race_asian", "theta_cl_egfr")),
           variable = fct_recode(variable, 
                                 "theta[CL[WT]]" = "theta_cl_wt",
                                 "theta[VC[WT]]" = "theta_vc_wt",
                                 "theta[Q[WT]]" = "theta_q_wt",
                                 "theta[VP[WT]]" = "theta_vp_wt",
                                 "theta[VC[ASIAN]]" = "theta_vc_race_asian",
                                 "theta[CL[eGFR]]" = "theta_cl_egfr")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_omega <- draws_all_df %>% 
    filter(str_detect(variable, "omega_")) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("cl", "vc", "q", "vp"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[Q]" = "omega_q",
                                 "omega[VP]" = "omega_vp")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor <- draws_all_df %>% 
    filter(variable %in% c("cor_cl_vc", "cor_cl_q", "cor_cl_vp",
                           "cor_vc_q", "cor_vc_vp", "cor_q_vp")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("cor_cl_vc", "cor_cl_q", "cor_cl_vp",
                                        "cor_vc_q", "cor_vc_vp", "cor_q_vp")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', Q)]" = "cor_cl_q",
                                 "rho[paste(CL, ', ', VP)]" = "cor_cl_vp",
                                 "rho[paste(VC, ', ', Q)]" = "cor_vc_q",
                                 "rho[paste(VC, ', ', VP)]" = "cor_vc_vp",
                                 "rho[paste(Q, ', ', VP)]" = "cor_q_vp")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))


(target_comparison_error <- draws_all_df %>% 
    filter(variable %in% c("sigma_p", "sigma_a", "cor_p_a")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("sigma_p", "sigma_a", "cor_p_a")),
           variable = fct_recode(variable, "sigma[prop]" = "sigma_p",
                                 "sigma[add]" = "sigma_a",
                                 "rho[paste(p, ', ', a)]" = "cor_p_a")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))


layout <- c(
  area(t = 1, l = 1, b = 1.5, r = 6),
  area(t = 2, l = 1, b = 2.5, r = 6),
  area(t = 3, l = 1, b = 3.5, r = 6),
  area(t = 4, l = 1, b = 4.5, r = 6),
  area(t = 5, l = 1, b = 5.5, r = 6)
)

target_comparison_tv /
  target_comparison_covariates /
  target_comparison_omega /
  target_comparison_cor /
  target_comparison_error +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom")

