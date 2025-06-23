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
TVKA <- 1
TVDUR <- 1

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3
omega_dur <- 0.3

R <- diag(rep(1, times = 4))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 0.2
sigma_a <- 0.5

cor_p_a <- 0

n_subjects_per_dose <- 6

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 1, amt = c(50, 100, 200, 400), ss = 0, 
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

dense_grid <- seq(0, 24*7, by = 0.5)

sampling_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times <- c(sampling_times, 72, 144, 144 + sampling_times)

times_to_simulate <- realistic_times # dense_grid

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
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
                  ii = nonmem_data_simulate$II,
                  addl = nonmem_data_simulate$ADDL,
                  ss = nonmem_data_simulate$SS,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVKA = TVKA,
                  TVDUR = TVDUR,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  omega_dur = omega_dur,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  solver = 4,
                  t_1 = 144,
                  t_2 = 168)

model <- cmdstan_model(
  "zero_into_depot_1cmt_linear/Stan/Simulate/zero_into_depot_1cmt_ppa_with_cmax_tmax_auc.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- if(stan_data$solver == 4){
  
  simulated_data$draws(c("CL", "VC", "KA", "DUR",
                         "auc_t1_t2", "c_max", "t_max", "t_half")) %>% 
    spread_draws(c(CL, VC, KA, DUR,
                   auc_t1_t2, c_max, t_max, t_half)[ID]) %>%
    ungroup() %>% 
    left_join(dosing_data %>% 
                filter(EVID == 1) %>% 
                distinct(ID, AMT),
              by = "ID") %>% 
    ungroup() %>%
    select(ID, AMT, CL, VC, KA, DUR, auc_t1_t2, c_max, t_max, t_half)
  
}else{
  
  simulated_data$draws(c("CL", "VC", "KA", "DUR", "t_half")) %>% 
    spread_draws(c(CL, VC, KA, DUR, t_half)[ID]) %>%
    ungroup() %>% 
    left_join(dosing_data %>% 
                filter(EVID == 1) %>% 
                distinct(ID, AMT),
              by = "ID") %>% 
    ungroup() %>%
    select(ID, AMT, CL, VC, KA, DUR, t_half)
  
}

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, CMT, EVID, SS, TIME, 
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

prop_or_ppa <- if_else(sigma_a == 0, "prop", "ppa")

data %>%
  select(-IPRED) %>%
  write_csv(str_c("zero_into_depot_1cmt_linear/Data/zero_into_depot_1cmt_",
                  prop_or_ppa, ".csv"),
            na = ".")

params_ind %>%
  write_csv(str_c("zero_into_depot_1cmt_linear/Data/zero_into_depot_1cmt_",
                  prop_or_ppa, "_params_ind.csv"),
            na = ".")


