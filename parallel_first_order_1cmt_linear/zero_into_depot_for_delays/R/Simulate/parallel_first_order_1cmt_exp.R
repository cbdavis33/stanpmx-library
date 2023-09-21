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

omega_cl <- 0.3
omega_vc <- 0.3

n_depots <- 3
delay_on_fastest_depot <- 1 # 0/1

n_subjects_per_dose <- 6

dosing_data <- expand.ev(ID = 1:n_subjects_per_dose, addl = 0, ii = 24, 
                         cmt = 1, amt = c(100, 200, 400, 800), ss = 0, tinf = 0,
                         evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

if(n_depots == 2){
  
  TVKA <- c(0.01, 0.05)
  TVFRAC <- c(0.5, 0.5)
  TVDUR <- c(168, 2)
  
  omega_ka <- c(0.3, 0.3)
  omega_dur <- c(0.4, 0.3)
  
  dosing_data <- dosing_data %>% 
    bind_rows(dosing_data %>%
                mutate(CMT = 2))
  
}else if(n_depots == 3){
  
  TVKA <- c(0.001, 0.005, 0.05)
  TVFRAC <- c(0.5, 0.4, 0.1)
  TVDUR <- c(1344, 336, 2)
  
  omega_ka <- c(0.3, 0.3, 0.3)
  omega_dur <- c(0.3, 0.3, 0.3)
  
  dosing_data <- dosing_data %>% 
    bind_rows(dosing_data %>%
                mutate(CMT = 2),
              dosing_data %>%
                mutate(CMT = 3))
  
}else if(n_depots == 4){
  
  TVKA <- c(0.0001, 0.001, 0.01, 0.1)
  TVFRAC <- c(0.5, 0.3, 0.15, 0.05)
  TVDUR <- c(504, 168, 48, 2)
  
  omega_ka <- c(0.3, 0.3, 0.3, 0.3)
  omega_dur <- c(0.3, 0.3, 0.3, 0.3)
  
  dosing_data <- dosing_data %>% 
    bind_rows(dosing_data %>%
                mutate(CMT = 2),
              dosing_data %>%
                mutate(CMT = 3), 
              dosing_data %>%
                mutate(CMT = 4))
  
}

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
                                 "parallel_first_order_1cmt_exp.stan"))

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "KA", "DUR")) %>% 
  spread_draws(CL[ID], VC[ID], KA[ID, Depot], DUR[ID, Depot]) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Depot", values_from = c("KA", "DUR")) %>% 
  select(ID, CL, VC, starts_with("KA"), starts_with("DUR"))

data <- simulated_data$draws(c("dv", "ipred")) %>% 
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
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(data$TIME), by = 24*7),
                       labels = seq(0, max(data$TIME)/(24*7), by = 1),
                       limits = c(0, max(data$TIME))))

plotly::ggplotly(p_1 +
                   scale_y_continuous(name = "Drug Conc. (ug/mL)",
                                      trans = "log10"))


# p_1 +
#   facet_trelliscope(~ID, nrow = 2, ncol = 2)

delay_string <- if_else(delay_on_fastest_depot == 0, "no_initial_delay",
                        "initial_delay")

data %>%
  select(-IPRED) %>%
  write_csv(file.path("parallel_first_order_1cmt_linear",
                      "zero_into_depot_for_delays",
                      "Data", 
                      str_c("parallel_", n_depots, "_first_order_",
                            delay_string, "_exp.csv")),
            na = ".")

params_ind %>%
  write_csv(file.path("parallel_first_order_1cmt_linear",
                      "zero_into_depot_for_delays",
                      "Data", 
                      str_c("parallel_", n_depots, "_first_order_",
                            delay_string, "_exp_params_ind.csv")),
            na = ".")




