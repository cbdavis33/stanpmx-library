rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(vpc)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("depot_1cmt_linear_ir4/Stan/Fits/depot_1cmt_exp_ir4_exp.rds")

nonmem_data <- read_csv(
  "depot_1cmt_linear_ir4/Data/depot_1cmt_exp_ir4_exp.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# VPCs can be done two ways, and they're the same up to Monte Carlo error
#   1) Use the _predict_new_subjects.stan file to simulate new datasets with
#        the same design. 
#   2) Use the fitted object that has already done this in generated quantities

new_data <- nonmem_data %>% 
  arrange(ID, time, evid, cmt) %>% 
  distinct(ID, time, cmt, .keep_all = TRUE) %>% 
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
                  want_auc_cmax = 0)

model <- cmdstan_model(
  "depot_1cmt_linear_ir4/Stan/Predict/depot_1cmt_exp_ir4_exp_predict_new_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

preds_col <- preds_df %>% 
  spread_draws(epred[i]) %>% 
  summarize(epred_mean = mean(epred)) %>% 
  ungroup() %>% 
  bind_cols(new_data %>% 
              select(ID, time, evid, cmt)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, cmt, epred_mean, -i) %>% 
  mutate(cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing")))

sim <- preds_df %>%
  spread_draws(c(epred_stan, epred)[i]) %>% 
  ungroup() %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         mdv = new_data$mdv[i],
         cmt = new_data$cmt[i],
         evid = new_data$evid[i],
         bloq = new_data$bloq[i]) %>% 
  filter(evid == 0) %>%
  # mutate(epred = if_else(bloq == 1, 0, epred)) %>% 
  select(ID, sim = ".draw", time, cmt, epred) %>% 
  mutate(cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing"))) %>% 
  arrange(ID, sim, time, cmt) %>% 
  left_join(preds_col, by = c("ID", "time", "cmt", "cmtn")) 

obs <- nonmem_data %>% 
  select(ID, time, dv = "DV", cmt, amt, evid, mdv, bloq, lloq) %>% 
  filter(evid == 0) %>%
  mutate(dv = if_else(bloq == 1, lloq, dv),
         cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing"))) %>% 
  left_join(preds_col, by = c("ID", "time", "cmt", "cmtn")) 

(p_vpc_pk <- vpc(sim = sim %>% 
                   filter(cmtn == "PK"), 
                 obs = obs %>% 
                   filter(cmtn == "PK"),
                 sim_cols = list(idv = "time",
                                 dv = "epred",
                                 id = "ID",
                                 pred = "epred_mean",
                                 sim = "sim"),
                 obs_cols = list(dv = "dv",
                                 idv = "time",
                                 id = "ID",
                                 pred = "epred_mean"),
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

(p_pcvpc_pk <- vpc(sim = sim %>% 
                     filter(cmtn == "PK"), 
                   obs = obs %>% 
                     filter(cmtn == "PK"),
                   sim_cols = list(idv = "time",
                                   dv = "epred",
                                   id = "ID",
                                   pred = "epred_mean",
                                   sim = "sim"),
                   obs_cols = list(dv = "dv",
                                   idv = "time",
                                   id = "ID",
                                   pred = "epred_mean"),
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
    # scale_y_continuous(name = "Response (units)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

p_vpc_pk + 
  p_pcvpc_pk

(p_vpc_pd <- vpc(sim = sim %>% 
                   filter(cmtn == "PD"), 
                 obs = obs %>% 
                   filter(cmtn == "PD"),
                 sim_cols = list(idv = "time",
                                 dv = "epred",
                                 id = "ID",
                                 pred = "epred_mean",
                                 sim = "sim"),
                 obs_cols = list(dv = "dv",
                                 idv = "time",
                                 id = "ID",
                                 pred = "epred_mean"),
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

(p_pcvpc_pd <- vpc(sim = sim %>% 
                     filter(cmtn == "PD"), 
                   obs = obs %>% 
                     filter(cmtn == "PD"),
                   sim_cols = list(idv = "time",
                                   dv = "epred",
                                   id = "ID",
                                   pred = "epred_mean",
                                   sim = "sim"),
                   obs_cols = list(dv = "dv",
                                   idv = "time",
                                   id = "ID",
                                   pred = "epred_mean"),
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

p_vpc_pd + 
  p_pcvpc_pd

p_pcvpc_pk +
  p_pcvpc_pd

(p_pcvpc <- vpc(sim = sim, obs = obs,
                sim_cols = list(idv = "time",
                                dv = "epred",
                                id = "ID",
                                pred = "epred_mean",
                                sim = "sim"),
                obs_cols = list(dv = "dv",
                                idv = "time",
                                id = "ID",
                                pred = "epred_mean"),
                pred_corr = TRUE,
                show = list(obs_dv = TRUE, 
                            obs_ci = TRUE,
                            pi = TRUE,
                            pi_as_area = FALSE,
                            pi_ci = TRUE,
                            obs_median = TRUE,
                            sim_median = TRUE,
                            sim_median_ci = TRUE),
                stratify = "cmtn",
                log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())


## Alternatively, use the output in the fitted object 

draws_df <- fit$draws(format = "draws_df")

preds_col <- draws_df %>% 
  spread_draws(epred[i]) %>% 
  summarize(epred_mean = mean(epred)) %>% 
  ungroup() %>% 
  bind_cols(nonmem_data %>% 
              filter(evid == 0) %>% 
              select(ID, time, evid, cmt)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, cmt, epred_mean, -i) %>% 
  mutate(cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing")))

sim <- draws_df %>%
  spread_draws(c(epred_stan, epred)[i]) %>% 
  ungroup() %>% 
  mutate(ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i],
         mdv = nonmem_data$mdv[nonmem_data$evid == 0][i],
         cmt = nonmem_data$cmt[nonmem_data$evid == 0][i],
         evid = nonmem_data$evid[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i]) %>% 
  filter(evid == 0) %>%
  # mutate(epred = if_else(bloq == 1, 0, epred)) %>% 
  select(ID, sim = ".draw", time, cmt, epred) %>% 
  mutate(cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing"))) %>% 
  arrange(ID, sim, time, cmt) %>% 
  left_join(preds_col, by = c("ID", "time", "cmt", "cmtn")) 

obs <- nonmem_data %>% 
  select(ID, time, dv = "DV", cmt, amt, evid, mdv, bloq, lloq) %>% 
  filter(evid == 0) %>%
  mutate(dv = if_else(bloq == 1, lloq, dv),
         cmtn = factor(case_when(cmt == 2 ~ "PK",
                                 cmt == 3 ~ "PD",
                                 TRUE ~ "Dosing"),
                       levels = c("PK", "PD", "Dosing"))) %>% 
  left_join(preds_col, by = c("ID", "time", "cmt", "cmtn")) 

(p_vpc_pk <- vpc(sim = sim %>% 
                   filter(cmtn == "PK"), 
                 obs = obs %>% 
                   filter(cmtn == "PK"),
                 sim_cols = list(idv = "time",
                                 dv = "epred",
                                 id = "ID",
                                 pred = "epred_mean",
                                 sim = "sim"),
                 obs_cols = list(dv = "dv",
                                 idv = "time",
                                 id = "ID",
                                 pred = "epred_mean"),
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


(p_pcvpc_pk <- vpc(sim = sim %>% 
                     filter(cmtn == "PK"), 
                   obs = obs %>% 
                     filter(cmtn == "PK"),
                   sim_cols = list(idv = "time",
                                   dv = "epred",
                                   id = "ID",
                                   pred = "epred_mean",
                                   sim = "sim"),
                   obs_cols = list(dv = "dv",
                                   idv = "time",
                                   id = "ID",
                                   pred = "epred_mean"),
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
    # scale_y_continuous(name = "Response (units)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

p_vpc_pk + 
  p_pcvpc_pk

(p_vpc_pd <- vpc(sim = sim %>% 
                   filter(cmtn == "PD"), 
                 obs = obs %>% 
                   filter(cmtn == "PD"),
                 sim_cols = list(idv = "time",
                                 dv = "epred",
                                 id = "ID",
                                 pred = "epred_mean",
                                 sim = "sim"),
                 obs_cols = list(dv = "dv",
                                 idv = "time",
                                 id = "ID",
                                 pred = "epred_mean"),
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

(p_pcvpc_pd <- vpc(sim = sim %>% 
                     filter(cmtn == "PD"), 
                   obs = obs %>% 
                     filter(cmtn == "PD"),
                   sim_cols = list(idv = "time",
                                   dv = "epred",
                                   id = "ID",
                                   pred = "epred_mean",
                                   sim = "sim"),
                   obs_cols = list(dv = "dv",
                                   idv = "time",
                                   id = "ID",
                                   pred = "epred_mean"),
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

p_vpc_pd + 
  p_pcvpc_pd

p_pcvpc_pk +
  p_pcvpc_pd

(p_pcvpc <- vpc(sim = sim, obs = obs,
                sim_cols = list(idv = "time",
                                dv = "epred",
                                id = "ID",
                                pred = "epred_mean",
                                sim = "sim"),
                obs_cols = list(dv = "dv",
                                idv = "time",
                                id = "ID",
                                pred = "epred_mean"),
                pred_corr = TRUE,
                show = list(obs_dv = TRUE, 
                            obs_ci = TRUE,
                            pi = TRUE,
                            pi_as_area = FALSE,
                            pi_ci = TRUE,
                            obs_median = TRUE,
                            sim_median = TRUE,
                            sim_median_ci = TRUE),
                stratify = "cmtn",
                log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())
