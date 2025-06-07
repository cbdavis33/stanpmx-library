rm(list = ls())
cat("\014")

library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "bioav_iv_and_oral_2cmt_linear/Stan/Fits/Stan_Data/ppa.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model("bioav_iv_and_oral_2cmt_linear/Stan/Fit/bioav_iv_and_oral_2cmt_ppa.stan",
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
                       init = function() 
                         with(stan_data,
                              list(TVCL = rlnorm(1, log(0.6), 0.3),
                                   TVVC = rlnorm(1, log(18), 0.3),
                                   TVQ = rlnorm(1, log(2), 0.3),
                                   TVVP = rlnorm(1, log(40), 0.3),
                                   TVKA = rlnorm(1, log(1), 0.3),
                                   TVBIOAV = rbeta(1, location_tvbioav*scale_tvbioav, 
                                                   (1 - location_tvbioav)*scale_tvbioav),
                                   omega = rlnorm(6, log(0.3), 0.3),
                                   sigma = rlnorm(2, log(0.4), 0.3))))

fit <- read_rds("bioav_iv_and_oral_2cmt_linear/Stan/Fits/bioav_iv_and_oral_2cmt_ppa.rds")

draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "q", "vp", "ka", "bioav")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_p", "sigma_a")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma_a, names_to = "variable", values_to = "value")

draws_for_tv <- draws_all_df %>% 
  filter(str_detect(variable, "TV")) %>%
  mutate(variable = factor(variable, 
                           levels = str_c("TV", c("CL", "VC", "Q", "VP", "KA", 
                                                  "BIOAV"))))

p_cl_vc_ka <- draws_for_tv %>%
  filter(variable %in% c("TVCL", "TVVC", "TVQ", "TVVP", "TVKA")) %>% 
  ggplot() +
  geom_density(aes(x = value, fill = target), alpha = 0.25) +
  theme_bw() +
  scale_fill_manual(name = "Distribution",
                    values = c("prior" = "blue", "posterior" = "red")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  facet_wrap(~ variable, scales = "free", nrow = 1)

p_bioav <- draws_for_tv %>%
  filter(variable %in% c("TVBIOAV")) %>% 
  ggplot() +
  geom_density(aes(x = value, fill = target), alpha = 0.25) +
  theme_bw() +
  scale_fill_manual(name = "Distribution",
                    values = c("prior" = "blue", "posterior" = "red")) +
  theme(legend.position = "bottom") +
  facet_wrap(~ variable, scales = "free", nrow = 1)

layout_tv <- c(
  area(t = 1, l = 1, b = 2, r = 4.9),
  area(t = 1, l = 5, b = 2, r = 6))

(target_comparison_tv <-  p_cl_vc_ka +
    p_bioav +
    plot_layout(guides = 'collect', 
                design = layout_tv) &
    theme(legend.position = "bottom"))

(target_comparison_omega <- draws_all_df %>% 
    filter(str_detect(variable, "omega_")) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", c("cl", "vc", "q", "vp", 
                                                        "ka", "bioav"))),
           variable = fct_recode(variable, 
                                 "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[Q]" = "omega_q",
                                 "omega[VP]" = "omega_vp",
                                 "omega[KA]" = "omega_ka",
                                 "omega[BIOAV]" = "omega_bioav")) %>% 
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
                    levels = c("cor_cl_vc", "cor_cl_q", "cor_cl_vp", "cor_cl_ka", "cor_cl_bioav",
                               "cor_vc_q", "cor_vc_vp", "cor_vc_ka", "cor_vc_bioav",
                               "cor_q_vp", "cor_q_ka", "cor_q_bioav",
                               "cor_vp_ka", "cor_vp_bioav",
                               "cor_ka_bioav")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', Q)]" = "cor_cl_q",
                                 "rho[paste(CL, ', ', VP)]" = "cor_cl_vp",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(CL, ', ', BIOAV)]" = "cor_cl_bioav",
                                 "rho[paste(VC, ', ', Q)]" = "cor_vc_q",
                                 "rho[paste(VC, ', ', VP)]" = "cor_vc_vp",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka",
                                 "rho[paste(VC, ', ', BIOAV)]" = "cor_vc_bioav",
                                 "rho[paste(Q, ', ', VP)]" = "cor_q_vp",
                                 "rho[paste(Q, ', ', KA)]" = "cor_q_ka",
                                 "rho[paste(Q, ', ', BIOAV)]" = "cor_q_bioav",
                                 "rho[paste(VP, ', ', KA)]" = "cor_vp_ka",
                                 "rho[paste(VP, ', ', BIOAV)]" = "cor_vp_bioav",
                                 "rho[paste(KA, ', ', BIOAV)]" = "cor_ka_bioav")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 2, labeller = label_parsed))


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
  area(t = 4, l = 1, b = 4.5, r = 6)
)

target_comparison_tv /
  target_comparison_omega /
  target_comparison_cor /
  target_comparison_error +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom") 
