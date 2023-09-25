rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 20
TVQ <- 4
TVVP <- 40
TVKA <- 1

FMAX <- 0.8 # Maximum Bioavailability Inhibition. Bioav at infinity is 1 - FMAX
F50 <- 300  # Amount of dose above reference dose that has half of full inhibition
            #   e.g. at reference_dose + F50, bioav = 1 - FMAX/2
HILL <- 1   # Adjust steepness of inhibition

omega_cl <- 0.3
omega_vc <- 0.3
omega_q <- 0.3
omega_vp <- 0.3
omega_ka <- 0.3

R <- diag(rep(1, times = 5))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 0.2
sigma_a <- 0.5

cor_p_a <- 0

n_subjects_per_dose <- 10

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                         cmt = 1, amt = c(100, 200, 400, 800), 
                         ss = 0, tinf = 0, evid = 1) %>%
  as_tibble() %>% 
  rename_with(toupper, .cols = everything()) %>%
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
         # RATE = 0,
         TIME = times_to_simulate) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT) %>% 
  rename_with(toupper, .cols = everything()) %>% 
  group_by(ID) %>% 
  mutate(DOSE = max(AMT)) %>% 
  ungroup()

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

dose <- nonmem_data_simulate %>% 
  group_by(ID) %>%
  distinct(DOSE) %>% 
  ungroup() %>% 
  pull(DOSE)

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
                  FMAX = FMAX,
                  F50 = F50,
                  HILL = HILL,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_q = omega_q,
                  omega_vp = omega_vp,
                  omega_ka = omega_ka,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  reference_dose = nonmem_data_simulate %>% 
                    filter(EVID == 1) %>%
                    pull(AMT) %>% 
                    min(),
                  dose = dose,
                  solver = 3) # analytical = 1, mat exp = 2, rk45 = 3

model <- cmdstan_model(
  "bioav_nonlinear_dose_2cmt_linear/Stan/Simulate/bioav_nonlinear_dose_2cmt_ppa.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235813,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "Q", "VP", "KA", "BIOAV")) %>% 
  spread_draws(CL[ID], VC[ID], Q[ID], VP[ID], KA[ID], BIOAV[ID]) %>% 
  inner_join(dosing_data %>%
               select(ID, DOSE = AMT) %>% 
               distinct(),
             by = "ID") %>%
  ungroup() %>%
  select(ID, CL, VC, Q, VP, KA, BIOAV, DOSE)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, DOSE,
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 1, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0),
         DOSE = factor(DOSE)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  rename_with(toupper, .cols = everything())

(p_1 <- ggplot(data %>% 
                 filter(!is.na(DV))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = DOSE)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = DOSE)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(data$TIME), by = 24),
                       labels = seq(0, max(data$TIME), by = 24),
                       limits = c(0, max(data$TIME))))
# p_1 +
#   facet_trelliscope(~ID, nrow = 2, ncol = 2)
# 

if(sigma_a == 0){

  data %>%
    select(-IPRED) %>%
    write_csv(file.path("bioav_nonlinear_dose_2cmt_linear",
                        "Data", 
                        "bioav_nonlinear_dose_2cmt_prop.csv"),
              na = ".")
  
  
  params_ind %>%
    write_csv(file.path("bioav_nonlinear_dose_2cmt_linear",
                        "Data", 
                        "bioav_nonlinear_dose_2cmt_prop_params_ind.csv"),
              na = ".")
  
}else{
  
  data %>%
    select(-IPRED) %>%
    write_csv(file.path("bioav_nonlinear_dose_2cmt_linear",
                        "Data", 
                        "bioav_nonlinear_dose_2cmt_ppa.csv"),
              na = ".")
  
  
  params_ind %>%
    write_csv(file.path("bioav_nonlinear_dose_2cmt_linear",
                        "Data", 
                        "bioav_nonlinear_dose_2cmt_ppa_params_ind.csv"),
              na = ".")
  
}
