rm(list = ls())
cat("\014")

# library(trelliscopejs)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_2cmt_linear_ir1/Data/depot_2cmt_ppa_ir1_ppa.csv",
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
                  location_tvcl = 0.75,
                  location_tvvc = 18,
                  location_tvq = 3,
                  location_tvvp = 35,
                  location_tvka = 0.8,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  scale_sigma_a = 1,
                  lkj_df_sigma = 2,
                  location_tvkin = 5,
                  location_tvkout = 0.2,
                  location_tvic50 = 10,
                  scale_tvkin = 1,
                  scale_tvkout = 1,
                  scale_tvic50 = 1,
                  scale_omega_kin = 0.4,
                  scale_omega_kout = 0.4,
                  scale_omega_ic50 = 0.4,
                  lkj_df_omega_pd = 2,
                  scale_sigma_p_pd = 0.5,
                  scale_sigma_a_pd = 1,
                  lkj_df_sigma_pd = 2,
                  prior_only = 1)

model <- cmdstan_model(
  "depot_2cmt_linear_ir1/Stan/Fit/depot_2cmt_ppa_ir1_ppa.stan",
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
                         TVQ = rlnorm(1, log(2), 0.3),
                         TVVP = rlnorm(1, log(40), 0.3),
                         TVKA = rlnorm(1, log(1), 0.3),
                         omega = rlnorm(5, log(0.3), 0.3),
                         sigma = c(rlnorm(1, log(0.2), 0.3), 
                                   rlnorm(1, log(0.8), 0.3)),
                         TVKIN = rlnorm(1, log(4), 0.3),
                         TVKOUT = rlnorm(1, log(0.3), 0.3),
                         TVIC50 = rlnorm(1, log(15), 0.3),
                         omega_pd = rlnorm(3, log(0.35), 0.3),
                         sigma_pd = c(rlnorm(1, log(0.2), 0.3), 
                                      rlnorm(1, log(0.8), 0.3))))


fit <- read_rds("depot_2cmt_linear_ir1/Stan/Fits/depot_2cmt_ppa_ir1_ppa.rds")
draws_df <- fit$draws(format = "draws_df")

parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "q", "vp", "ka",
                                               "kin", "kout", "ic50")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             str_c("sigma_", c("p", "a", "p_pd", "a_pd")))

draws_all_df <- priors$draws(format = "draws_df") %>% 
  mutate(target = "prior") %>% 
  drop_na() %>% 
  bind_rows(draws_df %>% 
              mutate(target = "posterior")) %>% 
  select(all_of(parameters_to_summarize), target) %>% 
  pivot_longer(cols = TVCL:sigma_a_pd, names_to = "variable", 
               values_to = "value")

(target_comparison_tv_pk <- draws_all_df %>% 
    filter(variable %in% str_c("TV", c("CL", "VC", "Q", 
                                       "VP", "KA"))) %>%
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
    filter(variable %in% str_c("omega_", c("cl", "vc", "q", 
                                           "vp", "ka"))) %>% 
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
    facet_wrap(~ variable, scales = "free", nrow = 2, labeller = label_parsed))

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
    filter(variable %in% c("sigma_p_pd", "sigma_a_pd", "cor_p_a_pd")) %>% 
    mutate(variable = factor(variable, 
                             levels = c("sigma_p_pd", "sigma_a_pd", "cor_p_a_pd")),
           variable = fct_recode(variable, "sigma[prop[PD]]" = "sigma_p_pd",
                                 "sigma[add[PD]]" = "sigma_a_pd",
                                 "rho[paste(p, ', ', a[PD])]" = "cor_p_a_pd")) %>% 
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
  area(t = 5, l = 1, b = 5.5, r = 6)
)

target_comparison_tv_pd /
  target_comparison_omega_pd /
  target_comparison_cor_pd /
  target_comparison_error_pd +
  plot_layout(guides = 'collect', 
              design = layout_pd) &
  theme(legend.position = "bottom")
