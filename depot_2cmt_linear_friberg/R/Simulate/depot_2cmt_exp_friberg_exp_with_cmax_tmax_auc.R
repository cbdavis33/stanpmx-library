rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 20
TVQ <- 4
TVVP <- 40
TVKA <- 1

TVMTT <- 125
TVCIRC0 <- 5
TVGAMMA <- 0.17
TVALPHA <- 3e-4

omega_cl <- 0.3
omega_vc <- 0.3
omega_q <- 0.3
omega_vp <- 0.3
omega_ka <- 0.3

omega_mtt <- 0.4
omega_circ0 <- 0.4
omega_gamma <- 0.4
omega_alpha <- 0.4

R <- diag(rep(1, times = 5))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

R_pd <- diag(rep(1, times = 4))
R_pd[1, 2] <- R_pd[2, 1] <- 0.7 # Put in some correlation between MTT and CIRC0

sigma <- 0.2

sigma_pd <- 0.15

n_subjects_per_dose <- 6

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24,
                         cmt = 1, amt = c(10, 20, 40, 80)*1000, tinf = 0,
                         ss = 0, evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

realistic_times_pk <- c(24, 48, 72, 72.25, 72.5, 72.75, 73, 73.5, 74, 74.5, 
                        75, 76, 78, 80, 84, 88, 92, 96, 120, 144, 144.25, 
                        144.5, 144.75, 145, 145.5, 146, 146.5, 147, 148, 150, 
                        152, 156, 160, 164, 168, 180, 192, 204, 216)

realistic_times_pd <- seq(0, 672, by = 24)

dense_grid_pk <- seq(0, 24*7, by = 2)
dense_grid_pd <- seq(0, 672, by = 2)

# times_to_simulate_pk <- dense_grid_pk
# times_to_simulate_pd <- dense_grid_pd

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
  arrange(ID, TIME, AMT, CMT)

nonmem_data_simulate_pd <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_pd))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 4,
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
                  TVQ = TVQ,
                  TVVP = TVVP,
                  TVKA = TVKA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_q = omega_q,
                  omega_vp = omega_vp,
                  omega_ka = omega_ka,
                  TVMTT = TVMTT,
                  TVCIRC0 = TVCIRC0,
                  TVGAMMA = TVGAMMA,
                  TVALPHA = TVALPHA,
                  omega_mtt = omega_mtt,
                  omega_circ0 = omega_circ0,
                  omega_gamma = omega_gamma,
                  omega_alpha = omega_alpha,
                  R = R,
                  sigma = sigma,
                  R_pd = R_pd,
                  sigma_pd = sigma_pd,
                  t_1 = 144,
                  t_2 = 168)

model <- cmdstan_model(
  "depot_2cmt_linear_friberg/Stan/Simulate/depot_2cmt_exp_friberg_exp_with_cmax_tmax_auc.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               # seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "Q", "VP", "KA",
                                     "MTT", "CIRC0", "GAMMA", "ALPHA",
                                     "auc_t1_t2", "c_max", "t_max", 
                                     "t_half_alpha", "t_half_terminal",
                                     "r_min")) %>% 
  spread_draws(c(CL, VC, Q, VP, KA,
                 MTT, CIRC0, GAMMA, ALPHA,
                 auc_t1_t2, c_max, t_max,
                 t_half_alpha, t_half_terminal,
                 r_min)[i]) %>% 
  inner_join(dosing_data %>% 
               mutate(i = 1:n()),
             by = "i") %>% 
  ungroup() %>%
  select(ID, CL, VC, Q, VP, KA, MTT, CIRC0, GAMMA, ALPHA,
         auc_t1_t2, c_max, t_max,
         t_half_alpha, t_half_terminal, r_min)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(c(dv, ipred)[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = if_else(CMT == 2, 10, 1), 
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
                                                CMT == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(ng/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME[data$CMT == 2]), by = 24),
                       labels = seq(0, max(data$TIME[data$CMT == 2]), by = 24),
                       limits = c(0, max(data$TIME[data$CMT == 2]))) +
    facet_wrap(~ cmt, scales = "free"))

(p_pd <- ggplot(data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(AMT))) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), CMT == 4) %>% 
                  mutate(cmt = factor(case_when(CMT == 2 ~ "PK",
                                                CMT == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Neutrophils $(\\times 10^9/L)"),
                       trans = "identity") + 
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(data$TIME[data$CMT == 4]), by = 168),
                       labels = seq(0, max(data$TIME[data$CMT == 4])/168, 
                                    by = 168/168),
                       limits = c(0, max(data$TIME[data$CMT == 4]))) +
    facet_wrap(~ cmt, scales = "free"))

p_pk + 
  p_pd +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")

data %>%
  select(-IPRED) %>% 
  write_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_exp_friberg_exp.csv", na = ".")

params_ind %>%
  write_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_exp_friberg_exp_params_ind.csv")


