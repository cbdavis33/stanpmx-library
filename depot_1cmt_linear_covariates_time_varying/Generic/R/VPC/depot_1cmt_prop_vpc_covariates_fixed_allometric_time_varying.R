rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(vpc)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

locf <- FALSE

file_path <- if_else(isTRUE(locf), 
                     "depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_prop_covariates_time_varying_generic_locf.csv",
                     "depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_prop_covariates_time_varying_generic_nocb.csv")


fit <- if(isTRUE(locf)){
  read_rds("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/prop_fixed_allometric_locf.rds")
}else{
  read_rds("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/prop_fixed_allometric_nocb.rds")
}

stan_data <- jsonlite::read_json(
  str_c("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data/", 
        "prop_fixed_allometric_nocb.json")) %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

stan_data <- jsonlite::read_json(
  str_c("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data/", 
        "prop_fixed_allometric_", if_else(locf, "locf", "nocb"), ".json")) %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

nonmem_data <- read_csv(
  file_path,
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

new_data <- nonmem_data %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, time, .keep_all = TRUE) %>% 
  select(-DV)

n_subjects <- new_data %>%  # number of individuals
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_time_new <- nrow(new_data) # total number of time points at which to predict

subj_start <- new_data %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_time_new)  

wt <- nonmem_data %>% 
  pull(wt)

cmppi <- nonmem_data %>% 
  pull(cmppi)

egfr <- nonmem_data %>% 
  pull(egfr)

stan_data <- list(n_subjects = n_subjects,
                  n_subjects_new = n_subjects,
                  n_time_new = n_time_new,
                  time = new_data$time,
                  amt = new_data$amt,
                  cmt = new_data$cmt,
                  evid = new_data$evid,
                  rate = new_data$rate,
                  ii = new_data$ii,
                  addl = new_data$addl,
                  ss = new_data$ss,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  t_1 = 144,
                  t_2 = 168,
                  wt = wt,
                  cmppi = cmppi,
                  egfr = egfr,
                  want_auc_cmax = 0,
                  theta_cl_wt = stan_data$theta_cl_wt,
                  theta_vc_wt = stan_data$theta_vc_wt)

model <- cmdstan_model(
  "depot_1cmt_linear_covariates_time_varying/Generic/Stan/Predict/depot_1cmt_prop_predict_new_subjects_covariates_fixed_allometric_time_varying.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

preds_col <- preds_df %>% 
  spread_draws(pred[i]) %>% 
  summarize(pred_mean = mean(pred)) %>% 
  ungroup() %>% 
  bind_cols(new_data %>% 
              select(ID, time, evid, cmppi)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, cmppi, pred_mean, -i)

sim <- preds_df %>%
  spread_draws(ipred[i], dv[i]) %>% 
  ungroup() %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         mdv = new_data$mdv[i],
         evid = new_data$evid[i],
         bloq = new_data$bloq[i]) %>% 
  filter(evid == 0) %>%
  mutate(dv = if_else(bloq == 1, 0, dv)) %>% 
  select(ID, sim = ".draw", time, ipred, dv) %>% 
  arrange(ID, sim, time) %>% 
  left_join(preds_col, by = c("ID", "time"))

obs <- nonmem_data %>% 
  select(ID, time, dv = "DV", amt, evid, mdv, bloq) %>% 
  filter(evid == 0) %>% 
  mutate(dv = if_else(bloq == 1, 0, dv)) %>% 
  left_join(preds_col, by = c("ID", "time"))

(p_vpc <- vpc(sim = sim, obs = obs,
              sim_cols = list(idv = "time",
                              dv = "dv",
                              id = "ID",
                              pred = "pred",
                              sim = "sim"),
              obs_cols = list(dv = "dv",
                              idv = "time",
                              id = "ID",
                              pred = "pred"),
              show = list(obs_dv = TRUE, 
                          obs_ci = TRUE,
                          pi = TRUE,
                          pi_as_area = FALSE,
                          pi_ci = TRUE,
                          obs_median = TRUE,
                          sim_median = TRUE,
                          sim_median_ci = TRUE),
              log_y = TRUE,
              lloq = 1) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "identity") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

(p_pcvpc <- vpc(sim = sim, obs = obs,
                sim_cols = list(idv = "time",
                                dv = "dv",
                                id = "ID",
                                pred = "pred_mean",
                                sim = "sim"),
                obs_cols = list(dv = "dv",
                                idv = "time",
                                id = "ID",
                                pred = "pred_mean"),
                pred_corr = TRUE,
                show = list(obs_dv = TRUE, 
                            obs_ci = TRUE,
                            pi = TRUE,
                            pi_as_area = FALSE,
                            pi_ci = TRUE,
                            obs_median = TRUE,
                            sim_median = TRUE,
                            sim_median_ci = TRUE),
                log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

p_vpc + 
  p_pcvpc

(p_vpc_cmppi <- vpc(sim = sim, obs = obs,
                    sim_cols = list(idv = "time",
                                    dv = "dv",
                                    id = "ID",
                                    pred = "pred",
                                    sim = "sim"),
                    obs_cols = list(dv = "dv",
                                    idv = "time",
                                    id = "ID",
                                    pred = "pred"),
                    show = list(obs_dv = TRUE, 
                                obs_ci = TRUE,
                                pi = TRUE,
                                pi_as_area = FALSE,
                                pi_ci = TRUE,
                                obs_median = TRUE,
                                sim_median = TRUE,
                                sim_median_ci = TRUE),
                    log_y = TRUE,
                    stratify = "cmppi",
                    lloq = 1) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "identity") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

(p_pcvpc_cmppi <- vpc(sim = sim, obs = obs,
                      sim_cols = list(idv = "time",
                                      dv = "dv",
                                      id = "ID",
                                      pred = "pred_mean",
                                      sim = "sim"),
                      obs_cols = list(dv = "dv",
                                      idv = "time",
                                      id = "ID",
                                      pred = "pred_mean"),
                      pred_corr = TRUE,
                      show = list(obs_dv = TRUE, 
                                  obs_ci = TRUE,
                                  pi = TRUE,
                                  pi_as_area = FALSE,
                                  pi_ci = TRUE,
                                  obs_median = TRUE,
                                  sim_median = TRUE,
                                  sim_median_ci = TRUE),
                      stratify = "cmppi",
                      log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

p_vpc_cmppi + 
  p_pcvpc_cmppi
