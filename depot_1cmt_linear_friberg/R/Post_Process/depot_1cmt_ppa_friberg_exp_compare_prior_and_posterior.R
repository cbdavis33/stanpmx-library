rm(list = ls())
cat("\014")

# library(trelliscopejs)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "depot_1cmt_linear_friberg/Stan/Fits/Stan_Data/ppa_exp.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model(
  "depot_1cmt_linear_friberg/Stan/Fit/depot_1cmt_ppa_friberg_exp.stan",
  cpp_options = list(stan_threads = TRUE))

priors <- model$sample(
  data = stan_data,
  seed = 235813,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 24,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  init = function() list(TVCL = rlnorm(1, log(0.6), 0.3),
                         TVVC = rlnorm(1, log(18), 0.3),
                         TVKA = rlnorm(1, log(1), 0.3),
                         omega = rlnorm(3, log(0.3), 0.3),
                         sigma = c(rlnorm(1, log(0.2), 0.3), 
                                   rlnorm(1, log(0.8), 0.3)),
                         TVMTT = rlnorm(1, log(120), 0.3),
                         TVCIRC0 = rlnorm(1, log(5), 0.3),
                         TVGAMMA = rlnorm(1, log(0.2), 0.3),
                         TVALPHA = rlnorm(1, log(3e-4), 0.3),
                         omega_pd = rlnorm(4, log(0.35), 0.3),
                         sigma_pd = rlnorm(1, log(0.2), 0.3)))

fit <- read_rds("depot_1cmt_linear_friberg/Stan/Fits/depot_1cmt_ppa_friberg_exp.rds")

draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka",
                                               "mtt", "circ0", "gamma", 
                                               "alpha")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             str_c("sigma_", c("p", "a", "pd")))

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma_pd, names_to = "variable", 
               values_to = "value")

(target_comparison_tv_pk <- draws_all_df %>% 
    filter(variable %in% str_c("TV", c("CL", "VC", "KA"))) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("CL", "VC", "KA")))) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() +
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1))

(target_comparison_tv_pd <- draws_all_df %>% 
    filter(variable %in% str_c("TV", c("MTT", "CIRC0", "GAMMA", "ALPHA"))) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("MTT", "CIRC0", "GAMMA", 
                                                    "ALPHA")))) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() +
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1))

(target_comparison_omega_pk <- draws_all_df %>% 
    filter(variable %in% str_c("omega_", c("cl", "vc", "ka"))) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("cl", "vc", "ka"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[KA]" = "omega_ka")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_omega_pd <- draws_all_df %>% 
    filter(variable %in% str_c("omega_", c("mtt", "circ0", "gamma", 
                                           "alpha"))) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("mtt", "circ0", "gamma",
                                                        "alpha"))),
           variable = fct_recode(variable, "omega[MTT]" = "omega_mtt",
                                 "omega[Circ[0]]" = "omega_circ0",
                                 "omega[gamma]" = "omega_gamma",
                                 "omega[alpha]" = "omega_alpha")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor_pk <- draws_all_df %>% 
    filter(variable %in% c("cor_cl_vc", "cor_cl_ka", "cor_vc_ka")) %>% 
    mutate(variable = 
             factor(variable, 
                    levels = c("cor_cl_vc", "cor_cl_ka", "cor_vc_ka")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 2, labeller = label_parsed))

(target_comparison_cor_pd <- draws_all_df %>% 
    filter(variable %in% c("cor_mtt_circ0", "cor_mtt_gamma", "cor_mtt_alpha",
                           "cor_circ0_gamma", "cor_circ0_alpha",
                           "cor_gamma_alpha")) %>% 
    mutate(variable = 
             factor(variable, 
                    levels = c("cor_mtt_circ0", "cor_mtt_gamma", "cor_mtt_alpha",
                               "cor_circ0_gamma", "cor_circ0_alpha",
                               "cor_gamma_alpha")),
           variable = fct_recode(variable, 
                                 "rho[paste(MTT, ', ', Circ[0])]" = 
                                   "cor_mtt_circ0",
                                 "rho[paste(MTT, ', ', gamma)]" = 
                                   "cor_mtt_gamma",
                                 "rho[paste(MTT, ', ', alpha)]" = 
                                   "cor_mtt_alpha",
                                 "rho[paste(Circ[0], ', ', gamma)]" = 
                                   "cor_circ0_gamma",
                                 "rho[paste(Circ[0], ', ', alpha)]" = 
                                   "cor_circ0_alpha",
                                 "rho[paste(gamma, ', ', alpha)]" = 
                                   "cor_gamma_alpha")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))


(target_comparison_error_pk <- draws_all_df %>% 
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

(target_comparison_error_pd <- draws_all_df %>% 
    filter(variable %in% c("sigma_pd")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("sigma_pd")),
           variable = fct_recode(variable, "sigma[PD]" = "sigma_pd")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

layout_pk <- c(
  area(t = 1, l = 1, b = 1.5, r = 6),
  area(t = 2, l = 1, b = 2.5, r = 6),
  area(t = 3, l = 1, b = 4.5, r = 6),
  area(t = 5, l = 1, b = 5.5, r = 6)
)

target_comparison_tv_pk /
  target_comparison_omega_pk /
  target_comparison_cor_pk /
  target_comparison_error_pk +
  plot_layout(guides = 'collect', 
              design = layout_pk) &
  theme(legend.position = "bottom")

layout_pd <- c(
  area(t = 1, l = 1, b = 1.5, r = 6),
  area(t = 2, l = 1, b = 2.5, r = 6),
  area(t = 3, l = 1, b = 4.5, r = 6),
  area(t = 5, l = 3, b = 5.5, r = 4)
)

target_comparison_tv_pd /
  target_comparison_omega_pd /
  target_comparison_cor_pd /
  target_comparison_error_pd +
  plot_layout(guides = 'collect', 
              design = layout_pd) &
  theme(legend.position = "bottom")

