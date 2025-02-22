rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

mod <- mread(
  "depot_1cmt_linear_covariates_time_varying/Generic/mrgsolve/CPP/depot_1cmt_covariates_time_varying.cpp")

TVCL <- 0.5
TVVC <- 12
TVKA <- 1

theta_cl_wt <- 0.75
theta_vc_wt <- 1
theta_ka_cmppi <- -1.5 # -0.8
theta_cl_egfr <- 0 # 0.5

omega_cl <- 0 # 0.3
omega_vc <- 0 # 0.3
omega_ka <- 0 # 0.3

R <- diag(rep(1, times = 3))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 1e-9 # 0.2
sigma_a <- 0    # 0.5

cor_p_a <- 0

n_subjects_per_dose <- 2

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 1, amt = c(400), ss = 0, tinf = 0, 
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 


# dense_grid <- seq(0, 24*7, by = 0.5)
dense_grid <- sort(c(seq(0, 24*7, by = 0.5), 47.99))

sampling_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times <- c(sampling_times, 72, 144, 144 + sampling_times)

times_to_simulate <- dense_grid
# times_to_simulate <- realistic_times

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         RATE = 0,
         TIME = times_to_simulate) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  group_by(ID) %>% 
  # mutate(SEXF = rbinom(1, 1, 0.5),
  #        AGE = sample(18:80, 1),
  #        RACE = sample(1:4, 1),
  #        CMPPI = if_else(TIME < 72, 0, 1),
  #        EGFR_baseline = 90, # runif(1, 15, 120),
  #        EGFR = if_else(TIME == 0, EGFR_baseline, rnorm(n(), EGFR_baseline, 2)),
  #        WT_baseline = rlnorm(1, log(70), 0.2),
  #        WT_beta = rnorm(1, 0, 0.004),) %>% 
  mutate(SEXF = rbinom(1, 1, 0.5),
         AGE = sample(18:80, 1),
         RACE = sample(1:4, 1),
         CMPPI = if_else(TIME < 72, 0, 1),
         EGFR = if_else(TIME < 24, 80, 90),
         WT = if_else(TIME < 72, 70, 80)) %>% 
  ungroup() %>% 
  # select(-WT_baseline, -WT_beta, -EGFR_baseline) %>% 
  arrange(ID, TIME, AMT) %>% 
  filter(!(ID == 1 & between(TIME, 36, 47.99)))

n_subjects <- nonmem_data_simulate %>%  # number of individuals to simulate
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_total <- nrow(nonmem_data_simulate) # total number of time points at which to predict

subj_start <- nonmem_data_simulate %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total) 

wt <- nonmem_data_simulate %>% 
  pull(WT)

cmppi <- nonmem_data_simulate %>% 
  pull(CMPPI)

egfr <- nonmem_data_simulate %>% 
  pull(EGFR)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  time = nonmem_data_simulate$TIME,
                  amt = nonmem_data_simulate$AMT,
                  cmt = nonmem_data_simulate$CMT,
                  evid = nonmem_data_simulate$EVID,
                  rate = nonmem_data_simulate$RATE,
                  ii = nonmem_data_simulate$II,
                  addl = nonmem_data_simulate$ADDL,
                  ss = nonmem_data_simulate$SS,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  wt = wt,
                  cmppi = cmppi,
                  egfr = egfr,
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVKA = TVKA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  theta_cl_wt = theta_cl_wt,
                  theta_vc_wt = theta_vc_wt,
                  theta_ka_cmppi = theta_ka_cmppi,
                  theta_cl_egfr = theta_cl_egfr,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  solver = 3) # analytical = 1, mat exp = 2, rk45 = 3

nonmem_data_simulate_locf <- nonmem_data_simulate %>%
  group_by(ID) %>% 
  mutate(WT_locf = lag(WT, default = first(WT)),
         EGFR_locf = lag(EGFR, default = first(EGFR)),
         CMPPI_locf = lag(CMPPI, default = first(CMPPI))) %>% 
  ungroup()

n_subjects <- nonmem_data_simulate_locf %>%  # number of individuals to simulate
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_total <- nrow(nonmem_data_simulate_locf) # total number of time points at which to predict

subj_start <- nonmem_data_simulate_locf %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total) 

wt <- nonmem_data_simulate_locf %>% 
  pull(WT_locf)

cmppi <- nonmem_data_simulate_locf %>% 
  pull(CMPPI_locf)

egfr <- nonmem_data_simulate_locf %>% 
  pull(EGFR)

stan_data_locf <- list(n_subjects = n_subjects,
                       n_total = n_total,
                       time = nonmem_data_simulate_locf$TIME,
                       amt = nonmem_data_simulate_locf$AMT,
                       cmt = nonmem_data_simulate_locf$CMT,
                       evid = nonmem_data_simulate_locf$EVID,
                       rate = nonmem_data_simulate_locf$RATE,
                       ii = nonmem_data_simulate_locf$II,
                       addl = nonmem_data_simulate_locf$ADDL,
                       ss = nonmem_data_simulate_locf$SS,
                       subj_start = subj_start,
                       subj_end = subj_end,
                       wt = wt,
                       cmppi = cmppi,
                       egfr = egfr,
                       TVCL = TVCL,
                       TVVC = TVVC,
                       TVKA = TVKA,
                       omega_cl = omega_cl,
                       omega_vc = omega_vc,
                       omega_ka = omega_ka,
                       theta_cl_wt = theta_cl_wt,
                       theta_vc_wt = theta_vc_wt,
                       theta_ka_cmppi = theta_ka_cmppi,
                       theta_cl_egfr = theta_cl_egfr,
                       R = R,
                       sigma_p = sigma_p,
                       sigma_a = sigma_a,
                       cor_p_a = cor_p_a,
                       solver = 3) # analytical = 1, mat exp = 2, rk45 = 3

model <- cmdstan_model(
  "depot_1cmt_linear_covariates_time_varying/Generic/Stan/Simulate/depot_1cmt_ppa_linear_covariates_time_varying_generic.stan") 

simulated_data_nocb <- model$sample(data = stan_data,
                                    fixed_param = TRUE,
                                    seed = 11235,
                                    iter_warmup = 0,
                                    iter_sampling = 1,
                                    chains = 1,
                                    parallel_chains = 1)

simulated_data_locf <- model$sample(data = stan_data_locf,
                                    fixed_param = TRUE,
                                    seed = 11235,
                                    iter_warmup = 0,
                                    iter_sampling = 1,
                                    chains = 1,
                                    parallel_chains = 1)


params_ind_nocb <- nonmem_data_simulate %>%
  mutate(i = row_number()) %>% 
  inner_join(simulated_data_nocb$draws(c("CL", "VC", "KA")) %>%
               spread_draws(c(CL, VC, KA)[i]) ,
             by = "i") %>% 
  select(ID, TIME, SEXF, AGE, RACE, WT, CMPPI, EGFR, CL, VC, KA)

params_ind_locf <- nonmem_data_simulate_locf %>%
  mutate(i = row_number()) %>% 
  inner_join(simulated_data_locf$draws(c("CL", "VC", "KA")) %>%
               spread_draws(c(CL, VC, KA)[i]) ,
             by = "i") %>% 
  select(ID, TIME, SEXF, AGE, RACE, WT, CMPPI, EGFR, CL, VC, KA)


data_stan_nocb <- simulated_data_nocb$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         SEXF, AGE, RACE, WT, CMPPI, EGFR,
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0, # 1,
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(software = "stan",
         method = "nocb") %>% 
  inner_join(params_ind_nocb) %>% 
  select(ID, TIME, EVID, CMT, AMT, CL, VC, KA, 
         WT, EGFR, CMPPI, IPRED, software, method)

data_stan_locf <- simulated_data_locf$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate_locf %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         SEXF, AGE, RACE, WT, CMPPI, EGFR,
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0, # 1,
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(software = "stan",
         method = "locf") %>% 
  inner_join(params_ind_locf) %>% 
  select(ID, TIME, EVID, CMT, AMT, CL, VC, KA, 
         WT, EGFR, CMPPI, IPRED, software, method)

data_mrgsolve_nocb <- mod %>% 
  param(THETA_CL_WT = theta_cl_wt,
        THETA_VC_WT = theta_vc_wt,
        THETA_KA_CMPPI = theta_ka_cmppi,
        THETA_CL_EGFR = theta_cl_egfr) %>% 
  data_set(nonmem_data_simulate) %>% 
  mrgsim_df(carry_out = c("AMT", "EVID", "CMT", "WT", "CMPPI", "EGFR"), 
            nocb = TRUE) %>% 
  as_tibble() %>% 
  mutate(software = "mrgsolve",
         method = "nocb") %>% 
  select(ID, TIME, EVID, CMT, AMT, CL = CLi, VC = VCi, KA = KAi, 
         WT, EGFR, CMPPI, IPRED, software, method)

data_mrgsolve_locf <- mod %>% 
  param(THETA_CL_WT = theta_cl_wt,
        THETA_VC_WT = theta_vc_wt,
        THETA_KA_CMPPI = theta_ka_cmppi,
        THETA_CL_EGFR = theta_cl_egfr) %>% 
  data_set(nonmem_data_simulate) %>% 
  mrgsim_df(carry_out = c("AMT", "EVID", "CMT", "WT", "CMPPI", "EGFR"), 
            nocb = FALSE) %>% 
  as_tibble() %>% 
  mutate(software = "mrgsolve",
         method = "locf") %>% 
  select(ID, TIME, EVID, CMT, AMT, CL = CLi, VC = VCi, KA = KAi, 
         WT, EGFR, CMPPI, IPRED, software, method)

data <- bind_rows(data_stan_nocb, data_stan_locf,
                  data_mrgsolve_nocb, data_mrgsolve_locf) %>% 
  distinct()

# data_1 <- mget(str_subset(ls(), "^data_[a-z]")) %>%
#   bind_rows()


(p_1 <- ggplot(data %>% 
                 filter(EVID == 0)) +
    # geom_point(mapping = aes(x = TIME, y = IPRED, 
    #                          group = interaction(ID, software, method), 
    #                          color = interaction(ID, software, method),
    #                          linetype = software)) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = IPRED, 
                            group = interaction(ID, software, method), 
                            color = interaction(ID, software, method),
                            linetype = software)) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID), color = "magenta") +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))))

(p_2 <- ggplot(data %>% 
                 filter(EVID == 0)) +
    # geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = IPRED, 
                            group = interaction(ID, software, method), 
                            color = interaction(ID, software, method))) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID), color = "magenta") +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))) +
    facet_wrap(~software + method))


(p_3 <- ggplot(data %>% 
                 filter(EVID == 0)) +
    # geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = TIME, y = KA, 
                             group = interaction(software, method, ID), 
                             color = interaction(software, method, ID))) +
    # geom_line(mapping = aes(x = TIME, y = DV, group = ID), color = "magenta") +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))))


data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - stan_nocb_2)) +
  geom_point() +
  theme_bw() 

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - stan_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - stan_locf_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - mrgsolve_nocb_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - mrgsolve_nocb_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_1 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()

##########################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - stan_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - stan_locf_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - mrgsolve_nocb_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - mrgsolve_nocb_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_nocb_2 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()


#######################################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_1 - stan_locf_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_1 - mrgsolve_nocb_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_1 - mrgsolve_nocb_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_1 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_1 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()


#######################################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_2 - mrgsolve_nocb_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_2 - mrgsolve_nocb_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_2 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = stan_locf_2 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()

#######################################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_nocb_1 - mrgsolve_nocb_2)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_nocb_1 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_nocb_1 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()

#######################################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_nocb_2 - mrgsolve_locf_1)) +
  geom_point() +
  theme_bw()

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_nocb_2 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()

#######################################

data %>% 
  mutate(IPRED = round(IPRED, 4)) %>%
  filter(TIME > 0) %>% 
  select(ID, TIME, IPRED, software, method) %>% 
  pivot_wider(values_from = IPRED, names_from = c(software, method, ID)) %>% 
  ggplot(aes(x = TIME, y = mrgsolve_locf_1 - mrgsolve_locf_2)) +
  geom_point() +
  theme_bw()

