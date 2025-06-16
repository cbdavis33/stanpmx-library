rm(list = ls())
cat("\014")

# library(trelliscopejs)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "depot_1cmt_linear_ir1/Stan/Fits/Stan_Data/exp_exp.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model(
  "depot_1cmt_linear_ir1/Stan/Fit/depot_1cmt_exp_ir1_exp.stan",
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
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_pk = abs(rnorm(1, 0, scale_sigma_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVIC50 = rlnorm(1, log(location_tvic50), scale_tvic50),
              omega_pd = abs(rnorm(3, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_ic50))),
              sigma_pd = abs(rnorm(1, 0, scale_sigma_pd)))))


fit <- read_rds("depot_1cmt_linear_ir1/Stan/Fits/depot_1cmt_exp_ir1_exp.rds")
draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka",
                                               "kin", "kout", "ic50")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             c("sigma_pk", "sigma_pd"))

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
                             levels = str_c("TV", c("CL", "VC", "Q", 
                                                    "VP", "KA")))) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() +
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1))

(target_comparison_tv_pd <- draws_all_df %>% 
    filter(variable %in% str_c("TV", c("KIN", "KOUT", 
                                       "IC50"))) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("KIN", "KOUT", 
                                                    "IC50")))) %>% 
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
                             levels = str_c("omega_", c("cl", "vc", "q", 
                                                        "vp", "ka"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[Q]" = "omega_q",
                                 "omega[VP]" = "omega_vp",
                                 "omega[KA]" = "omega_ka")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_omega_pd <- draws_all_df %>% 
    filter(variable %in% str_c("omega_", c("kin", "kout", "ic50"))) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("kin", "kout", 
                                                        "ic50"))),
           variable = fct_recode(variable, "omega[K[`in`]]" = "omega_kin",
                                 "omega[K[out]]" = "omega_kout",
                                 "omega[IC[50]]" = "omega_ic50")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor_pk <- draws_all_df %>% 
    filter(variable %in% c("cor_cl_vc", "cor_cl_q", "cor_cl_vp", "cor_cl_ka",
                           "cor_vc_q", "cor_vc_vp", "cor_vc_ka",
                           "cor_q_vp", "cor_q_ka",
                           "cor_vp_ka")) %>% 
    mutate(variable = 
             factor(variable, 
                    levels = c("cor_cl_vc", "cor_cl_q", "cor_cl_vp", "cor_cl_ka",
                               "cor_vc_q", "cor_vc_vp", "cor_vc_ka",
                               "cor_q_vp", "cor_q_ka",
                               "cor_vp_ka")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', Q)]" = "cor_cl_q",
                                 "rho[paste(CL, ', ', VP)]" = "cor_cl_vp",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(VC, ', ', Q)]" = "cor_vc_q",
                                 "rho[paste(VC, ', ', VP)]" = "cor_vc_vp",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka",
                                 "rho[paste(Q, ', ', VP)]" = "cor_q_vp",
                                 "rho[paste(Q, ', ', KA)]" = "cor_q_ka",
                                 "rho[paste(VP, ', ', KA)]" = "cor_vp_ka")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor_pd <- draws_all_df %>% 
    filter(variable %in% c("cor_kin_kout", "cor_kin_ic50", "cor_kout_ic50")) %>% 
    mutate(variable = 
             factor(variable, 
                    levels = c("cor_kin_kout", "cor_kin_ic50", "cor_kout_ic50")),
           variable = fct_recode(variable, 
                                 "rho[paste(K[`in`], ', ', K[out])]" = 
                                   "cor_kin_kout",
                                 "rho[paste(K[`in`], ', ', IC[50])]" = 
                                   "cor_kin_ic50",
                                 "rho[paste(K[out], ', ', IC[50])]" = 
                                   "cor_kout_ic50")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))


(target_comparison_error_pk <- draws_all_df %>% 
    filter(variable %in% c("sigma_pk")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("sigma_pk")),
           variable = fct_recode(variable, "sigma[PK]" = "sigma_pk")) %>% 
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
  area(t = 5, l = 3, b = 5.5, r = 4)
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
