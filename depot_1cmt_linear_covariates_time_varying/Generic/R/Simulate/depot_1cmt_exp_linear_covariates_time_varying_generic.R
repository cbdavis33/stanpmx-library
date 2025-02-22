rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 12
TVKA <- 1

theta_cl_wt <- 0.75
theta_vc_wt <- 1
theta_ka_cmppi <- -0.8
theta_cl_egfr <- 0.5

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3

R <- diag(rep(1, times = 3))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma <- 0.2

locf <- TRUE # last-one-carried-forward (locf, like Monolix) or next-observation-carried-backward (nocb, like NONMEM) 

n_subjects_per_dose <- 25

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 1, amt = c(50, 100, 200, 400), ss = 0, tinf = 0, 
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

dense_grid <- seq(0, 24*7, by = 0.5)

sampling_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times <- c(sampling_times, 72, 144, 144 + sampling_times)

# times_to_simulate <- dense_grid
times_to_simulate <- realistic_times

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
  mutate(SEXF = rbinom(1, 1, 0.5),
         AGE = sample(18:80, 1),
         RACE = sample(1:4, 1),
         CMPPI = if_else(TIME < 72, 0, 1),
         EGFR_baseline = runif(1, 15, 120),
         EGFR = if_else(TIME == 0, EGFR_baseline, rnorm(n(), EGFR_baseline, 2)),
         WT_baseline = rlnorm(1, log(70), 0.2),
         WT_beta = rnorm(1, 0, 0.004),
         WT = WT_baseline + WT_beta*TIME) %>% 
  ungroup() %>% 
  select(-WT_baseline, -WT_beta, -EGFR_baseline) %>% 
  arrange(ID, TIME, AMT) 

if(isTRUE(locf)){
  nonmem_data_simulate <- nonmem_data_simulate %>%
    group_by(ID) %>%
    mutate(WT = lag(WT, default = first(WT)),
           EGFR = lag(EGFR, default = first(EGFR)),
           CMPPI = lag(CMPPI, default = first(CMPPI))) %>%
    ungroup()
}

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
                  sigma = sigma,
                  solver = 3) # analytical = 1, mat exp = 2, rk45 = 3

model <- cmdstan_model(
  "depot_1cmt_linear_covariates_time_varying/Generic/Stan/Simulate/depot_1cmt_exp_linear_covariates_time_varying_generic.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- nonmem_data_simulate %>%
  mutate(i = row_number()) %>% 
  inner_join(simulated_data$draws(c("CL", "VC", "KA")) %>%
               spread_draws(c(CL, VC, KA)[i]) ,
             by = "i") %>% 
  select(ID, TIME, SEXF, AGE, RACE, WT, CMPPI, EGFR, CL, VC, KA)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         SEXF, AGE, RACE, WT, CMPPI, EGFR,
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 1,
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) 

(p_1 <- ggplot(data %>% 
                 group_by(ID) %>% 
                 mutate(Dose = factor(max(AMT))) %>% 
                 ungroup() %>% 
                 filter(!is.na(DV))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))))
# p_1 +
#   facet_trelliscope(~ID, nrow = 2, ncol = 2)


((params_ind %>% 
    ggplot(aes(x = WT, y = CL)) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    theme_bw()) |
    (params_ind %>% 
       ggplot(aes(x = WT, y = VC)) + 
       geom_point() + 
       geom_smooth(method = "lm") +
       theme_bw())) /
  ((params_ind %>%
      mutate(CMPPI = factor(CMPPI)) %>% 
      ggplot(aes(x = CMPPI, y = KA, group = CMPPI)) + 
      geom_boxplot() + 
      theme_bw()) +
     (params_ind %>% 
        ggplot(aes(x = EGFR, y = CL)) + 
        geom_point() + 
        geom_smooth(method = "lm") +
        theme_bw()))

data %>%
  select(-IPRED) %>% 
  # write_csv("depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_exp_covariates_time_varying_generic_nocb.csv",
  #           na = ".")
  write_csv("depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_exp_covariates_time_varying_generic_locf.csv",
            na = ".")

params_ind %>%
  # write_csv("depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_exp_params_ind_covariates_time_varying_generic_nocb.csv")
  write_csv("depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_exp_params_ind_covariates_time_varying_generic_locf.csv")


