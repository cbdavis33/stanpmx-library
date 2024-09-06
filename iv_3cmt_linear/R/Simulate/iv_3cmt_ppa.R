rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 2.5    # L/m
TVVC <- 4.5    # L
TVQ1 <- 2      # L/m
TVVP1 <- 4     # L
TVQ2 <- 0.1    # L/m
TVVP2 <- 5     # L

omega_cl <- 0.3
omega_vc <- 0.3
omega_q1 <- 0.3
omega_vp1 <- 0.3
omega_q2 <- 0.3
omega_vp2 <- 0.3

R <- diag(rep(1, times = 6))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 0.2
sigma_a <- 0.05 # LLOQ = 0.1

cor_p_a <- 0

n_subjects_per_dose <- 6

solver <- 1 # mat exp = 1, rk45 = 2

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 0, ii = 0, 
                         cmt = 1, amt = c(750, 1000, 1500, 2500, 4000), ss = 0, 
                         tinf = 20, evid = 1) %>%
  as_tibble() %>% 
  mutate(cmt = 1) %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

dense_grid <- seq(0, 180, by = 1)

sampling_times <- c(5, 10, 20, 25, 30, seq(40, 120, by = 10), 140, 160, 180)
realistic_times <- sampling_times

# times_to_simulate <- dense_grid
times_to_simulate <- realistic_times

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 1,
         EVID = 0,
         RATE = 0,
         TIME = times_to_simulate) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT)

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
                  TVQ1 = TVQ1,
                  TVVP1 = TVVP1,
                  TVQ2 = TVQ2,
                  TVVP2 = TVVP2,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_q1 = omega_q1,
                  omega_vp1 = omega_vp1,
                  omega_q2 = omega_q2,
                  omega_vp2 = omega_vp2,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  solver = solver) # mat exp = 1, rk45 = 2

model <- cmdstan_model("iv_3cmt_linear/Stan/Simulate/iv_3cmt_ppa.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 112358,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "Q1", "VP1", "Q2", "VP2")) %>% 
  spread_draws(CL[i], VC[i], Q1[i], VP1[i], Q2[i], VP2[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>%
  select(ID, CL, VC, Q1, VP1, Q2, VP2)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 0.1, 
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
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(ng/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (m)",
                       breaks = seq(0, max(data$TIME), by = 15),
                       labels = seq(0, max(data$TIME), by = 15),
                       limits = c(0, max(data$TIME))))
p_1 +
  facet_trelliscope(~ID, nrow = 2, ncol = 2)

data %>%
  select(-IPRED) %>%
  # write_csv("iv_3cmt_linear/Data/iv_3cmt_prop.csv", na = ".")
  write_csv("iv_3cmt_linear/Data/iv_3cmt_ppa.csv", na = ".")

params_ind %>%
  # write_csv("iv_3cmt_linear/Data/iv_3cmt_prop_params_ind.csv")
  write_csv("iv_3cmt_linear/Data/iv_3cmt_ppa_params_ind.csv")


