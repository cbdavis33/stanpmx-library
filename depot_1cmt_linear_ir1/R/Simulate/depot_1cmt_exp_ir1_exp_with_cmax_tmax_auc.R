rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 12
TVKA <- 1

TVKIN <- 5
TVKOUT <- 0.2
TVIC50 <- 10

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3

omega_kin <- 0.4
omega_kout <- 0.4
omega_ic50 <- 0.4

R_pk <- diag(rep(1, times = 3))
R_pk[1, 2] <- R_pk[2, 1] <- 0.4 # Put in some correlation between CL and VC

R_pd <- diag(rep(1, times = 3))
R_pd[1, 2] <- R_pd[2, 1] <- 0.7 # Put in some correlation between KIN and KOUT

sigma_pk <- 0.2

sigma_pd <- 0.15

n_subjects_per_dose <- 6

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 1, amt = c(100, 200, 400, 800), ss = 0, tinf = 0, 
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 


dense_grid <- seq(0, 24*7, by = 2)

sampling_times_pk <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times_pk <- c(sampling_times_pk, 72, 144, 144 + sampling_times_pk)

sampling_times_pd <- c(0, 0.25, 1, 2, 6, 12, 24)
realistic_times_pd <- c(sampling_times_pd, 72, 144, 144 + sampling_times_pd)

# times_to_simulate_pk <- dense_grid
# times_to_simulate_pd <- dense_grid

times_to_simulate_pk <- realistic_times_pk
times_to_simulate_pd <- realistic_times_pd

nonmem_data_simulate_pk <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_pk))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_pk) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT)

nonmem_data_simulate_pd <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_pd))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 3,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_pd) %>% 
  ungroup() %>%
  arrange(ID, TIME, AMT)

nonmem_data_simulate <- nonmem_data_simulate_pk %>% 
  bind_rows(nonmem_data_simulate_pd) %>% 
  arrange(ID, TIME, EVID, CMT)

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
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVKA = TVKA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  TVKIN = TVKIN,
                  TVKOUT = TVKOUT,
                  TVIC50 = TVIC50,
                  omega_kin = omega_kin,
                  omega_kout = omega_kout,
                  omega_ic50 = omega_ic50,
                  R_pk = R_pk,
                  sigma_pk = sigma_pk,
                  R_pd = R_pd,
                  sigma_pd = sigma_pd,
                  t_1 = 144,
                  t_2 = 168,
                  coupled = 0,
                  want_auc_cmax = 1,
                  solver = 1) # 1 = rk45, 2 = bdf, 3 = adams

model <- cmdstan_model(
  "depot_1cmt_linear_ir1/Stan/Simulate/depot_1cmt_exp_ir1_exp_with_cmax_tmax_auc.stan")

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               # seed = 112358,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- if(stan_data$want_auc_cmax){
  
  simulated_data$draws(c("CL", "VC", "KA",
                         "KIN", "KOUT", "IC50",
                         "auc_t1_t2", "c_max", "t_max", 
                         "t_half", "r_min", "t_min")) %>% 
    spread_draws(c(CL, VC, KA, KIN, KOUT, IC50,
                   auc_t1_t2, c_max, t_max,
                   t_half, r_min, t_min)[i]) %>% 
    inner_join(dosing_data %>% 
                 mutate(i = 1:n()),
               by = "i") %>% 
    ungroup() %>%
    select(ID, CL, VC, KA, KIN, KOUT, IC50,
           auc_t1_t2, c_max, t_max,
           t_half, r_min, t_min)
  
}else{
  
  simulated_data$draws(c("CL", "VC", "KA",
                         "KIN", "KOUT", "IC50")) %>% 
    spread_draws(c(CL, VC, KA, KIN, KOUT, IC50)[i]) %>% 
    inner_join(dosing_data %>% 
                 mutate(i = 1:n()),
               by = "i") %>% 
    ungroup() %>%
    select(ID, CL, VC, KA, KIN, KOUT, IC50)
  
}

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(c(dv, ipred)[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = if_else(CMT == 2, 1, 2), 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV)

(p_pk <- ggplot(data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(AMT))) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), CMT == 2) %>% 
                  mutate(cmt = factor(case_when(CMT == 2 ~ "PK",
                                                CMT == 3 ~ "PD",
                                                # .default = "Dosing"),
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))) +
    facet_wrap(~ cmt, scales = "free"))

(p_pd <- ggplot(data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(AMT))) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), CMT == 3) %>% 
                  mutate(cmt = factor(case_when(CMT == 2 ~ "PK",
                                                CMT == 3 ~ "PD",
                                                # .default = "Dosing"),
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Response (unit)"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))) +
    facet_wrap(~ cmt, scales = "free"))

p_pk + 
  p_pd +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")

data %>%
  select(-IPRED) %>%
  write_csv(str_c("depot_1cmt_linear_ir1/Data/depot_1cmt_exp_ir1_exp.csv"),
            na = ".")

params_ind %>%
  write_csv(str_c("depot_1cmt_linear_ir1/Data/depot_1cmt_exp_ir1_exp_params_ind.csv"),
            na = ".")
