rm(list = ls())
cat("\014")

# library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

stan_data <- jsonlite::read_json(
  "depot_1cmt_linear_plus_mm/Stan/Fits/Stan_Data/exp.json") %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

model <- cmdstan_model(
  "depot_1cmt_linear_plus_mm/Stan/Fit/depot_1cmt_linear_plus_mm_exp.stan",
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
                       init = function() list(TVCL = rlnorm(1, log(0.6), 0.3),
                                              TVVC = rlnorm(1, log(10), 0.3),
                                              TVVMAX = rlnorm(1, log(42), 0.3),
                                              TVKM = rlnorm(1, log(4), 0.3),
                                              TVKA = rlnorm(1, log(1), 0.3),
                                              omega = rlnorm(5, log(0.3), 0.3),
                                              sigma = rlnorm(1, log(0.2), 0.3)))

fit <- read_rds(
  "depot_1cmt_linear_plus_mm/Stan/Fits/depot_1cmt_linear_plus_mm_exp.rds")

draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "vmax", "km", "ka")),
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
                             levels = str_c("TV", c("CL", "VC", "VMAX", 
                                                    "KM", "KA")))) %>% 
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
                             levels = str_c("omega_", c("cl", "vc", "vmax", 
                                                        "km", "ka"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[VMAX]" = "omega_vmax",
                                 "omega[KM]" = "omega_km",
                                 "omega[KA]" = "omega_ka")) %>% 
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
                             levels = c("cor_cl_vc", "cor_cl_vmax", 
                                        "cor_cl_km", "cor_cl_ka",
                                        "cor_vc_vmax", "cor_vc_km", "cor_vc_ka",
                                        "cor_vmax_km", "cor_vmax_ka", "cor_km_ka")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', VMAX)]" = "cor_cl_vmax",
                                 "rho[paste(CL, ', ', KM)]" = "cor_cl_km",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(VC, ', ', VMAX)]" = "cor_vc_vmax",
                                 "rho[paste(VC, ', ', KM)]" = "cor_vc_km",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka",
                                 "rho[paste(VMAX, ', ', KM)]" = "cor_vmax_km",
                                 "rho[paste(VMAX, ', ', KA)]" = "cor_vmax_ka",
                                 "rho[paste(KM, ', ', KA)]" = "cor_km_ka")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))


(target_comparison_error <- draws_all_df %>% 
    filter(variable %in% c("sigma")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("sigma")),
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
  area(t = 3, l = 1, b = 3.5, r = 6),
  area(t = 4, l = 3, b = 4.5, r = 4)
)

target_comparison_tv /
  target_comparison_omega /
  target_comparison_cor /
  target_comparison_error +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom")
