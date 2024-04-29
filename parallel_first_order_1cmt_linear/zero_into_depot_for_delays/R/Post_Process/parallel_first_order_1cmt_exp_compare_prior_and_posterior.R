rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

n_depots <- 3
n_depots_with_delay <- n_depots # "n_depots - 1" means no delay on the 
                                    # fastest depot. "n_depots" means a delay on
                                    # fastest depot

with_or_without_delay_string <- if_else(n_depots_with_delay == n_depots, 
                                        "with_initial_delay", 
                                        "no_initial_delay")

prior_string <- "prior_kan_durn"

stan_data <- jsonlite::read_json(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fits",
            "Stan_Data",
            str_c("parallel_", n_depots, "_first_order_1cmt_",
                  with_or_without_delay_string, "_", 
                  prior_string, "_exp.json"))) %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data$prior_only <- stan_data$no_gq_predictions <- 1

rsimplex <- function(n_depots){
  
  x <- 1 + runif(n_depots, -0.5, 0.5)
  
  return(x/sum(x))
  
}

model <- cmdstan_model(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fit",
            "parallel_first_order_1cmt_exp.stan"),
  cpp_options = list(stan_threads = TRUE))

priors <- model$sample(
  data = stan_data,
  seed = 235813,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  init = function() list(TVKA = sort(rlnorm(n_depots, log(0.01), 1)),
                         TVCL = rlnorm(1, log(1), 0.3),
                         TVVC = rlnorm(1, log(8), 0.3),
                         TVDUR_backward = sort(rlnorm(n_depots_with_delay, 
                                                      log(500), 0.7), 
                                               decreasing = FALSE),
                         TVFRAC = rsimplex(n_depots),
                         omega = rlnorm(n_depots + n_depots + 2, 
                                        log(0.3), 0.3),
                         sigma = rlnorm(1, log(0.2), 0.3)))

fit <- read_rds(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fits",
            str_c("parallel_", n_depots, "_first_order_1cmt_",
                  with_or_without_delay_string, "_", prior_string,
                  "_exp.rds")))

draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "back", 
                                        negate = TRUE) %>% 
                               str_subset("TV"),
                             str_c("omega_", c("cl", "vc", "ka", "dur")),
                             "sigma")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(contains(parameters_to_summarize)), target) %>%
  select(!contains("back")) %>%
  pivot_longer(cols = `TVKA[1]`:sigma, names_to = "variable", 
               values_to = "value")

draws_for_tv <- draws_all_df %>% 
  filter(str_detect(variable, "TV")) %>%
  mutate(variable = 
           factor(variable, 
                  levels = str_c("TV", 
                                 c("CL", "VC", 
                                   str_c("KA[", c(1:n_depots), "]"),
                                   str_c("DUR[", c(1:n_depots_with_delay), "]"),
                                   str_c("FRAC[", c(1:n_depots), "]")))))

p_cl_vc_ka_dur <- draws_for_tv %>% 
  filter(variable %in% c("TVCL", "TVVC") | str_detect(variable, "TVKA") |
           str_detect(variable, "TVDUR")) %>% 
  ggplot() +
  geom_density(aes(x = value, fill = target), alpha = 0.25) +
  theme_bw() +
  scale_fill_manual(name = "Distribution",
                    values = c("prior" = "blue", "posterior" = "red")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  facet_wrap(~ variable, scales = "free", nrow = 2)

p_frac <- draws_for_tv %>%
  filter(str_detect(variable, "TVFRAC")) %>% 
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

(target_comparison_tv <-  p_cl_vc_ka_dur +
    p_frac +
    plot_layout(guides = 'collect', 
                design = layout_tv) &
    theme(legend.position = "bottom"))


(target_comparison_omega <- draws_all_df %>% 
    filter(str_detect(variable, "omega_")) %>% 
    mutate(variable = factor(variable, 
                             levels = str_c("omega_", 
                                            c("cl", "vc", 
                                              str_c("ka[", c(1:n_depots), "]"),
                                              str_c("dur[", c(1:n_depots_with_delay), "]"))))) %>% 
    # variable = fct_recode(variable, "omega[CL]" = "omega_cl",
    #                       "omega[VC]" = "omega_vc",
    #                       "omega[KA]" = "omega_ka",
    #                       "omega[DUR]" = "omega_dur")) %>% 
    ggplot() +
    geom_density(aes(x = value, fill = target), alpha = 0.25) +
    theme_bw() + 
    scale_fill_manual(name = "Distribution",
                      values = c("prior" = "blue", "posterior" = "red")) +
    theme(legend.position = "bottom") +
    facet_wrap(~ variable, scales = "free", nrow = 1, labeller = label_parsed))

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
  area(t = 1, l = 1, b = 2.5, r = 6),
  # area(t = 2, l = 1, b = 2.5, r = 6),
  area(t = 3, l = 1, b = 3.5, r = 6),
  area(t = 4, l = 3, b = 4.5, r = 4)
)

target_comparison_tv /
  target_comparison_omega /
  target_comparison_sigma +
  plot_layout(guides = 'collect', 
              design = layout) &
  theme(legend.position = "bottom")
