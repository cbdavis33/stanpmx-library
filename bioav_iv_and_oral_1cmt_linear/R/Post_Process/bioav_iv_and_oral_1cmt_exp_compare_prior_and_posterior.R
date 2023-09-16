rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv(
  "bioav_iv_and_oral_1cmt_linear/Data/bioav_iv_and_oral_1cmt_exp.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

n_subjects <- nonmem_data %>%  # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data)   # total number of records

i_obs <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  filter(evid == 0) %>%
  select(row_num) %>%
  deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(row_num) %>%
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  rate = nonmem_data$rate,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  location_tvcl = 1,
                  location_tvvc = 8,
                  location_tvka = 0.8,
                  location_tvbioav = 0.5,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_tvbioav = 3.5,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  scale_omega_bioav = 0.3,
                  lkj_df_omega = 2,
                  scale_sigma = 0.5,
                  prior_only = 1)

model <- cmdstan_model("bioav_iv_and_oral_1cmt_linear/Stan/Fit/bioav_iv_and_oral_1cmt_exp.stan",
                       cpp_options = list(stan_threads = TRUE))

priors <- model$sample(data = stan_data,
                       seed = 11235,
                       chains = 4,
                       parallel_chains = 4,
                       threads_per_chain = parallel::detectCores()/4,
                       iter_warmup = 500,
                       iter_sampling = 1000,
                       adapt_delta = 0.8,
                       refresh = 500,
                       max_treedepth = 10,
                       init = function() list(TVCL = rlnorm(1, log(1), 0.3),
                                              TVVC = rlnorm(1, log(8), 0.3),
                                              TVKA = rlnorm(1, log(0.8), 0.3),
                                              TVBIOAV = rbeta(1, 3, 3),
                                              omega = rlnorm(4, log(0.3), 0.3),
                                              sigma = rlnorm(1, log(0.4), 0.3)))


fit <- read_rds("bioav_iv_and_oral_1cmt_linear/Stan/Fits/bioav_iv_and_oral_1cmt_exp.rds")
draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka", "bioav")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma")

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  drop_na() %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma, names_to = "variable", values_to = "value")


draws_for_tv <- draws_all_df %>% 
  filter(str_detect(variable, "TV")) %>%
  mutate(variable = factor(variable, 
                           levels = str_c("TV", c("CL", "VC", "KA", 
                                                  "BIOAV"))))

p_cl_vc_ka <- draws_for_tv %>%
  filter(variable %in% c("TVCL", "TVVC", "TVKA")) %>% 
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
                             levels = str_c("omega_", c("cl", "vc", "ka",
                                                        "bioav"))),
           variable = fct_recode(variable, "omega[CL]" = "omega_cl",
                                 "omega[VC]" = "omega_vc",
                                 "omega[KA]" = "omega_ka",
                                 "omega[F]" = "omega_bioav")) %>% 
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
                             levels = c("cor_cl_vc", "cor_cl_ka", "cor_cl_bioav",
                                        "cor_vc_ka", "cor_vc_bioav",
                                        "cor_ka_bioav")),
           variable = fct_recode(variable, 
                                 "rho[paste(CL, ', ', VC)]" = "cor_cl_vc",
                                 "rho[paste(CL, ', ', KA)]" = "cor_cl_ka",
                                 "rho[paste(CL, ', ', F)]" = "cor_cl_bioav",
                                 "rho[paste(VC, ', ', KA)]" = "cor_vc_ka",
                                 "rho[paste(VC, ', ', F)]" = "cor_vc_bioav",
                                 "rho[paste(KA, ', ', F)]" = "cor_ka_bioav")) %>% 
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
