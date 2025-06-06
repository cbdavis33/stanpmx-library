rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 12
TVKA <- 1
TVBIOAV <- 0.4

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3
omega_bioav <- 0.2

R <- diag(rep(1, times = 4))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC
R[2, 4] <- R[4, 2] <- 0.2 # Put in some correlation between VC and F

sigma_p <- 0.2
sigma_a <- 0.5

cor_p_a <- 0

n_subjects_per_dose <- 10

dosing_data_oral <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                              cmt = 1, amt = c(100, 300), ss = 0, tinf = 0, 
                              evid = 1) %>%
  as_tibble() %>% 
  rename_with(toupper, .cols = everything()) %>%
  select(ID, TIME, everything()) %>% 
  mutate(type = "Oral",
         cohort = if_else(ID <= n_subjects_per_dose, 1, 2))

dosing_data_iv <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                            cmt = 2, amt = c(100, 300), ss = 0, tinf = 2, 
                            evid = 1) %>%
  as_tibble() %>% 
  rename_with(toupper, .cols = everything()) %>%
  select(ID, TIME, everything()) %>% 
  mutate(cohort = if_else(ID <= n_subjects_per_dose, 3, 4),
         ID = max(dosing_data_oral$ID) + ID,
         type = "IV")

dosing_data_iv_and_oral_iv <- expand.ev(ID = 1:n_subjects_per_dose, addl = 1, 
                                        ii = 24, cmt = 2, amt = c(100, 300), 
                                        ss = 0, tinf = 2, evid = 1) %>% 
  as.ev()

dosing_data_iv_and_oral_oral <- expand.ev(ID = 1:(2*n_subjects_per_dose), 
                                          addl = 4, ii = 24, cmt = 1, amt = 300, 
                                          ss = 0, tinf = 0, evid = 1) %>% 
  as.ev()

dosing_data_iv_and_oral <- ev_seq(dosing_data_iv_and_oral_iv,
                                  dosing_data_iv_and_oral_oral) %>% 
  as_tibble() %>% 
  rename_with(toupper, .cols = everything()) %>%
  arrange(ID, TIME) %>% 
  mutate(cohort = if_else(ID <= n_subjects_per_dose, 5, 6),
         ID = max(dosing_data_iv$ID) + ID, 
         type = "IV and Oral")

dosing_data <- dosing_data_oral %>% 
  bind_rows(dosing_data_iv,
            dosing_data_iv_and_oral) %>% 
  mutate(cohort = factor(cohort)) %>% 
  rename_with(toupper, .cols = everything()) # %>% 
  # mutate(TYPE = "want dense grid")

dense_grid <- seq(0, 24*7, by = 0.5)

sampling_times_oral <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times_oral <- c(sampling_times_oral, 72, 144, 
                          144 + sampling_times_oral)

sampling_times_iv <- c(2, 3, 4, 6, 8, 12, 24)
realistic_times_iv <- c(sampling_times_iv, 72, 144, 
                        144 + sampling_times_iv)

realistic_times_iv_and_oral <- c(sampling_times_iv, 72, 144, 
                                 144 + sampling_times_oral)

nonmem_data_simulate <- dosing_data %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>% 
  mutate(TIME = case_when(TYPE == "IV" ~ list(realistic_times_iv),
                          TYPE == "Oral" ~ list(realistic_times_oral),
                          TYPE == "IV and Oral" ~ list(realistic_times_iv_and_oral),
                          .default = list(dense_grid))) %>% 
  ungroup() %>%
  unnest(TIME) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         RATE = 0,
         TINF = 0,
         EVID = 0) %>% 
  ungroup() %>%
  bind_rows(dosing_data) %>% 
  arrange(ID, TIME, AMT) %>% 
  rename_with(toupper, .cols = everything())

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
                  TVBIOAV = TVBIOAV,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  omega_bioav = omega_bioav,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  t_1 = 144,
                  t_2 = 168,
                  solver = 4)

model <- cmdstan_model(
  "bioav_iv_and_oral_1cmt_linear/Stan/Simulate/bioav_iv_and_oral_1cmt_ppa_with_cmax_tmax_auc.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 54321,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)


params_ind <- if(stan_data$solver == 4){
  
  simulated_data$draws(c("CL", "VC", "KA", "BIOAV",
                         "auc_t1_t2", "c_max", "t_max", "t_half")) %>% 
    spread_draws(c(CL, VC, KA, BIOAV,
                   auc_t1_t2, c_max, t_max, t_half)[ID]) %>%
    ungroup() %>% 
    left_join(dosing_data %>% 
                distinct(ID, TYPE, COHORT),
              by = "ID") %>% 
    ungroup() %>%
    select(ID, CL, VC, KA, BIOAV, auc_t1_t2, c_max, t_max, t_half, 
           TYPE, COHORT)
  
}else{
  
  simulated_data$draws(c("CL", "VC", "KA", "BIOAV")) %>% 
    spread_draws(c(CL, VC, KA, BIOAV)[ID]) %>%
    ungroup() %>% 
    left_join(dosing_data %>% 
                distinct(ID, TYPE, COHORT),
              by = "ID") %>% 
    ungroup() %>%
    select(ID, CL, VC, KA, BIOAV, TYPE, COHORT)
  
}


data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(dv[i], ipred[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, TYPE, COHORT,
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 1, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV) %>% 
  rename_with(toupper, .cols = everything())

(p_1 <- ggplot(data %>% 
                 filter(!is.na(DV))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = COHORT)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = COHORT)) +
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
  write_csv(str_c("bioav_iv_and_oral_1cmt_linear/Data/bioav_iv_and_oral_1cmt_",
                  prop_or_ppa, ".csv"),
            na = ".")

params_ind %>%
  write_csv(str_c("bioav_iv_and_oral_1cmt_linear/Data/bioav_iv_and_oral_1cmt_",
                  prop_or_ppa, "_params_ind.csv"),
            na = ".")


