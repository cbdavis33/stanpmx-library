rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 12

omega_cl <- 0
omega_vc <- 0

n_depots <- 3
delay_on_fastest_depot <- 0 # 0/1

n_subjects_per_dose <- 1

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 0, ii = 24, 
                         cmt = 1, amt = 800, ss = 0, tinf = 0,
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 



TVKA <- c(0.001, 0.005, 0.05)
TVFRAC <- c(0.5, 0.4, 0.1)
TVDUR <- c(1344, 336, 2)

omega_ka <- c(0, 0, 0)
omega_dur <- c(0, 0, 0)

dosing_data <- dosing_data %>% 
  bind_rows(dosing_data %>%
              mutate(CMT = 2),
            dosing_data %>%
              mutate(CMT = 3))

if(delay_on_fastest_depot == 0){
  
  omega_dur[n_depots] <- 0
  TVDUR[n_depots] <- 0
  
}

R <- diag(rep(1, times = n_depots*2 + 2))
R[n_depots + 1, n_depots + 2] <- R[n_depots + 2, n_depots + 1] <- 0.4 # Put in some correlation between CL and VC

sigma <- 0.2

dense_grid <- c(seq(0, 167, by = 1), 
                seq(168, 502, by = 2), 
                seq(504, 2000, by = 4))

realistic_times <- c(1, 2, 4, 6, 8, 12, 24, 48, 72, 96, 120, 144, 168, 216, 312, 
                     384, 480, 552, 648, 816, 984, 1152, 1320, 1488, 1656, 1824,
                     2160, 2496, 2832)

# times_to_simulate <- dense_grid
times_to_simulate <- realistic_times

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = n_depots + 1,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, CMT, AMT)

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
                  ii = nonmem_data_simulate$II,
                  addl = nonmem_data_simulate$ADDL,
                  ss = nonmem_data_simulate$SS,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  n_depots = n_depots,
                  delay_on_fastest_depot = delay_on_fastest_depot,
                  TVKA = TVKA,
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVDUR = TVDUR,
                  TVFRAC = TVFRAC,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  omega_dur = omega_dur,
                  R = R,
                  sigma = sigma,
                  solver = 1) # mat exp = 1, rk45 = 2

model <- cmdstan_model(file.path("parallel_first_order_1cmt_linear",
                                 "zero_into_depot_for_delays",
                                 "Stan", "Simulate",
                                 "parallel_3_first_order_1cmt_exp.stan"))

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

data_1 <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0.01, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(solver = "mat_exp",
         type = "fixed_3")


stan_data$solver <- 2
simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

data_2 <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0.01, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(solver = "ode",
         type = "fixed_3")



stan_data$solver <- 1
model <- cmdstan_model(file.path("parallel_first_order_1cmt_linear",
                                 "zero_into_depot_for_delays",
                                 "Stan", "Simulate",
                                 "parallel_first_order_1cmt_exp.stan"))

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

data_3 <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0.01, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(solver = "mat_exp",
         type = "adjustable_3")


stan_data$solver <- 2
simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

data_4 <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0.01, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  mutate(solver = "ode",
         type = "adjustable_3")

mget(str_subset(ls(), "data_[0-9]")) %>% 
  bind_rows() %>% 
  filter(EVID == 0) %>% 
  mutate(DV = if_else(BLOQ == 1, LLOQ/2, DV)) %>% 
  ggplot(aes(x = TIME, y = IPRED, group = type, color = solver)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~solver)

mget(str_subset(ls(), "data_[0-9]")) %>% 
  bind_rows() %>% 
  filter(EVID == 0) %>% 
  mutate(DV = if_else(BLOQ == 1, LLOQ/2, DV)) %>% 
  ggplot(aes(x = TIME, y = IPRED, group = solver, color = type)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~type)

mget(str_subset(ls(), "data_[0-9]")) %>% 
  bind_rows() %>% 
  filter(EVID == 0) %>% 
  mutate(DV = if_else(BLOQ == 1, LLOQ/2, DV)) %>% 
  ggplot(aes(x = TIME, y = IPRED, group = interaction(solver, type), 
             color = interaction(solver, type))) +
  geom_line() +
  theme_bw() 
