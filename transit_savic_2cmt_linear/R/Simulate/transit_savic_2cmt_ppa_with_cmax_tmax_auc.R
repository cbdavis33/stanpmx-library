rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 1
TVVC <- 20
TVQ <- 4
TVVP <- 40
TVKA <- 2
TVNTR <- 6
TVMTT <- 1

omega_cl <- 0.3
omega_vc <- 0.3
omega_q <- 0.3
omega_vp <- 0.3
omega_ka <- 0.3
omega_ntr <- 0.3
omega_mtt <- 0.3

R <- diag(rep(1, times = 7))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 0.2
sigma_a <- 0.5

cor_p_a <- 0

n_subjects_per_dose <- 4

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 1, ii = 24, 
                         cmt = 1, amt = c(100, 200, 400, 800), ss = 0, tinf = 0,
                         evid = 1) %>%
  realize_addl() %>% 
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything())

dense_grid <- c(seq(0.1, 1.9, by = 0.1), seq(2, 24, by = 1),
                seq(24.1, 25.9, by = 0.1), seq(26, 48, by = 1))


sampling_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times <- c(sampling_times, 24 + sampling_times)

times_to_simulate <- realistic_times # dense_grid

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         TIME = times_to_simulate) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT) %>% 
  select(-TINF)

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

nonmem_data_dose <- nonmem_data_simulate %>% 
  filter(EVID == 1)

n_dose <- nonmem_data_dose %>% 
  nrow()

subj_start_dose <- nonmem_data_dose %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end_dose <- c(subj_start_dose[-1] - 1, n_dose) 

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
                  TVNTR = TVNTR,
                  TVMTT = TVMTT,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_q = omega_q,
                  omega_vp = omega_vp,
                  omega_ka = omega_ka,
                  omega_ntr = omega_ntr,
                  omega_mtt = omega_mtt,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  n_dose = n_dose,
                  dosetime = nonmem_data_dose$TIME,
                  doseamt = nonmem_data_dose$AMT,
                  subj_start_dose = subj_start_dose,
                  subj_end_dose = subj_end_dose,
                  t_1 = 0,
                  t_2 = 24,
                  solver = 4)

model <- cmdstan_model(
  "transit_savic_2cmt_linear/Stan/Simulate/transit_savic_2cmt_ppa_with_cmax_tmax_auc.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235813,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "Q", "VP", "KA", "NTR", "MTT", 
                                     "auc_t1_t2", "c_max", "t_max",
                                     "t_half_alpha", "t_half_terminal")) %>% 
  spread_draws(CL[ID], VC[ID], Q[ID], VP[ID], KA[ID], NTR[ID], MTT[ID],
               auc_t1_t2[ID], c_max[ID], t_max[ID], 
               t_half_alpha[ID], t_half_terminal[ID]) %>% 
  inner_join(dosing_data %>%
               group_by(ID) %>% 
               slice_head(n = 1) %>%
               ungroup(),
             by = "ID") %>% 
  ungroup() %>%
  select(ID, CL, VC, Q, VP, KA, NTR, MTT, auc_t1_t2, c_max, t_max, 
         t_half_alpha, t_half_terminal)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
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
p_1 +
  facet_trelliscope(~ID, nrow = 2, ncol = 2)

data %>%
  select(-IPRED) %>%
  write_csv("transit_savic_2cmt_linear/Data/transit_savic_2cmt_ppa.csv", na = ".")
# write_csv("transit_savic_2cmt_linear/Data/transit_savic_2cmt_prop.csv", na = ".")

params_ind %>%
  write_csv("transit_savic_2cmt_linear/Data/transit_savic_2cmt_ppa_params_ind.csv")
# write_csv("transit_savic_2cmt_linear/Data/transit_savic_2cmt_prop_params_ind.csv")

