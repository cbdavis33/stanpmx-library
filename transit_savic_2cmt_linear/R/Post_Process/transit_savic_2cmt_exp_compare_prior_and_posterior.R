rm(list = ls())
cat("\014")

library(cmdstanr)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "transit_savic_2cmt_linear/Stan/Fits/Stan_Data/exp.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model(
  "transit_savic_2cmt_linear/Stan/Fit/transit_savic_2cmt_exp.stan",
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
                       init = function() list(TVCL = rlnorm(1, log(0.4), 0.3),
                                              TVVC = rlnorm(1, log(15), 0.3),
                                              TVQ = rlnorm(1, log(6), 0.3),
                                              TVVP = rlnorm(1, log(35), 0.3),
                                              TVKA = rlnorm(1, log(4), 0.3),
                                              TVNTR = rlnorm(1, log(5), 0.3),
                                              TVMTT = rlnorm(1, log(1), 0.3),
                                              omega = rlnorm(7, log(0.3), 0.3),
                                              sigma = rlnorm(1, log(0.2), 0.3)))


fit <- read_rds("transit_savic_2cmt_linear/Stan/Fits/transit_savic_2cmt_exp.rds")
draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_c("TV", c("CL", "VC", "Q", "VP", "KA",
                                           "NTR", "MTT")),
                             str_c("omega_", c("cl", "vc", "q", "vp", "ka",
                                               "ntr", "mtt")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma, names_to = "variable", values_to = "value")

(target_comparison_tv <- draws_all_df %>% 
    filter(str_detect(variable, "TV")) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("CL", "VC", "Q", "VP", "KA", 
                                                    "NTR", "MTT")))) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() +
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1))

(target_comparison_omega <- draws_all_df %>% 
    filter(str_detect(variable, "omega_")) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("cl", "vc", "q", "vp", 
                                                        "ka", "ntr", "mtt"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[Q]" = "omega_q",
                                 "omega[VP]" = "omega_vp",
                                 "omega[KA]" = "omega_ka",
                                 "omega[NTR]" = "omega_ntr",
                                 "omega[MTT]" = "omega_mtt")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor <- draws_all_df %>% 
    filter(str_detect(variable, "cor_")) %>% 
    mutate(variable = 
             factor(variable, 
                    levels = c("cor_cl_vc", "cor_cl_q", "cor_cl_vp", 
                               "cor_cl_ka", "cor_cl_ntr", "cor_cl_mtt",
                               "cor_vc_q", "cor_vc_vp", "cor_vc_ka", 
                               "cor_vc_ntr", "cor_vc_mtt",
                               "cor_q_vp", "cor_q_ka", "cor_q_ntr", "cor_q_mtt",
                               "cor_vp_ka", "cor_vp_ntr", "cor_vp_mtt",
                               "cor_ka_ntr", "cor_ka_mtt",
                               "cor_ntr_mtt")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', Q)]" = "cor_cl_q",
                                 "rho[paste(CL, ', ', VP)]" = "cor_cl_vp",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(CL, ', ', NTR)]" = "cor_cl_ntr",
                                 "rho[paste(CL, ', ', MTT)]" = "cor_cl_mtt",
                                 "rho[paste(VC, ', ', Q)]" = "cor_vc_q",
                                 "rho[paste(VC, ', ', VP)]" = "cor_vc_vp",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka",
                                 "rho[paste(VC, ', ', NTR)]" = "cor_vc_ntr",
                                 "rho[paste(VC, ', ', MTT)]" = "cor_vc_mtt",
                                 "rho[paste(Q, ', ', VP)]" = "cor_q_vp",
                                 "rho[paste(Q, ', ', KA)]" = "cor_q_ka",
                                 "rho[paste(Q, ', ', NTR)]" = "cor_q_ntr",
                                 "rho[paste(Q, ', ', MTT)]" = "cor_q_mtt",
                                 "rho[paste(VP, ', ', KA)]" = "cor_vp_ka",
                                 "rho[paste(VP, ', ', NTR)]" = "cor_vp_ntr",
                                 "rho[paste(VP, ', ', MTT)]" = "cor_vp_mtt",
                                 "rho[paste(KA, ', ', NTR)]" = "cor_ka_ntr",
                                 "rho[paste(KA, ', ', MTT)]" = "cor_ka_mtt",
                                 "rho[paste(NTR, ', ', MTT)]" = "cor_ntr_mtt")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 3, labeller = label_parsed))


(target_comparison_sigma <- draws_all_df %>% 
    filter(str_detect(variable, "sigma")) %>% 
    mutate(variable = factor(variable, 
                             levels = "sigma"),
           variable = fct_recode(variable, "sigma" = "sigma")) %>% 
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
  area(t = 3, l = 1, b = 4.5, r = 6),
  area(t = 5, l = 3, b = 5.5, r = 4)
)

target_comparison_tv /
  target_comparison_omega /
  target_comparison_cor /
  target_comparison_sigma +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom")

