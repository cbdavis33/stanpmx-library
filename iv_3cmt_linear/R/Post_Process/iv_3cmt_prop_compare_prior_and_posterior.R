rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "iv_3cmt_linear/Stan/Fits/Stan_Data/prop.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model("iv_3cmt_linear/Stan/Fit/iv_3cmt_prop.stan",
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
                       init = function() list(TVCL = rlnorm(1, log(2), 0.3),
                                              TVVC = rlnorm(1, log(5), 0.3),
                                              TVQ1 = rlnorm(1, log(2), 0.3),
                                              TVVP1 = rlnorm(1, log(3.5), 0.3),
                                              TVQ2 = rlnorm(1, log(0.8), 0.3),
                                              TVVP2 = rlnorm(1, log(4), 0.3),
                                              omega = rlnorm(6, log(0.3), 0.3),
                                              sigma_p = rlnorm(1, log(0.2), 0.3)))


fit <- read_rds("iv_3cmt_linear/Stan/Fits/iv_3cmt_prop.rds")
draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "q1", "vp1",
                                               "q2", "vp2")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_p")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma_p, names_to = "variable", values_to = "value")

(target_comparison_tv <- draws_all_df %>% 
    filter(str_detect(variable, "TV")) %>%
    mutate(variable = factor(variable, 
                             levels = str_c("TV", c("CL", "VC", "Q1", "VP1",
                                                    "Q2", "VP2")))) %>% 
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
                             levels = str_c("omega_", c("cl", "vc", "q1", "vp1",
                                                        "q2", "vp2"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[Q[1]]" = "omega_q1",
                                 "omega[VP[1]]" = "omega_vp1",
                                 "omega[Q[2]]" = "omega_q2",
                                 "omega[VP[2]]" = "omega_vp2")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

(target_comparison_cor <- draws_all_df %>% 
    filter(str_detect(variable, "cor_")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("cor_cl_vc", "cor_cl_q1", "cor_cl_vp1",
                                        "cor_cl_q2", "cor_cl_vp2",
                                        "cor_vc_q1", "cor_vc_vp1", "cor_vc_q2", 
                                        "cor_vc_vp2",
                                        "cor_q1_vp1", "cor_q1_q2", "cor_q1_vp2",
                                        "cor_vp1_q2", "cor_vp1_vp2",
                                        "cor_q2_vp2")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', Q1)]" = "cor_cl_q1",
                                 "rho[paste(CL, ', ', VP1)]" = "cor_cl_vp1",
                                 "rho[paste(CL, ', ', Q2)]" = "cor_cl_q2",
                                 "rho[paste(CL, ', ', VP2)]" = "cor_cl_vp2",
                                 "rho[paste(VC, ', ', Q1)]" = "cor_vc_q1",
                                 "rho[paste(VC, ', ', VP1)]" = "cor_vc_vp1",
                                 "rho[paste(VC, ', ', Q2)]" = "cor_vc_q2",
                                 "rho[paste(VC, ', ', VP2)]" = "cor_vc_vp2",
                                 "rho[paste(Q1, ', ', VP1)]" = "cor_q1_vp1",
                                 "rho[paste(Q1, ', ', Q2)]" = "cor_q1_q2",
                                 "rho[paste(Q1, ', ', VP2)]" = "cor_q1_vp2",
                                 "rho[paste(VP1, ', ', Q2)]" = "cor_vp1_q2",
                                 "rho[paste(VP1, ', ', VP2)]" = "cor_vp1_vp2",
                                 "rho[paste(Q2, ', ', VP2)]" = "cor_q2_vp2")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 3, labeller = label_parsed))


(target_comparison_sigma <- draws_all_df %>% 
    filter(str_detect(variable, "sigma_")) %>% 
    mutate(variable = factor(variable, 
                             levels = "sigma_p"),
           variable = fct_recode(variable, "sigma[prop]" = "sigma_p")) %>% 
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
  area(t = 4, l = 3, b = 4.5, r = 4)
)

target_comparison_tv /
  target_comparison_omega /
  target_comparison_cor /
  target_comparison_sigma +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom")

