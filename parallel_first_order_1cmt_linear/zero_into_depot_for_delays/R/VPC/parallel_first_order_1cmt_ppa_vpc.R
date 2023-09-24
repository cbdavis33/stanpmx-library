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

n_depots <- 3
n_depots_with_delay <- n_depots - 1 # "n_depots - 1" means no delay on the 
                                # fastest depot. "n_depots" means a delay on
                                # fastest depot

with_or_without_delay_string <- if_else(n_depots_with_delay == n_depots, 
                                        "with_initial_delay", 
                                        "no_initial_delay")
prior_string <- "prior_kan_dur1"

nonmem_data <- read_csv(file.path("parallel_first_order_1cmt_linear",
                                  "zero_into_depot_for_delays",
                                  "Data", 
                                  str_c("parallel_", n_depots, "_first_order_",
                                        with_or_without_delay_string, 
                                        "_ppa.csv")),
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

fit <- read_rds(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fits",
            str_c("parallel_", n_depots, "_first_order_1cmt_",
                  with_or_without_delay_string, "_", prior_string,
                  "_ppa.rds")))

new_data <- nonmem_data %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, cmt, time, .keep_all = TRUE) %>% 
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
                  ii = new_data$ii,
                  addl = new_data$addl,
                  ss = new_data$ss,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  n_depots = n_depots,
                  n_depots_with_delay = n_depots_with_delay,
                  t_1 = 0,
                  t_2 = 1500)

model <- cmdstan_model(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Predict",
            "parallel_first_order_1cmt_ppa_predict_new_subjects.stan"))

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
              select(ID, time, evid)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, pred_mean, -i)

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
              lloq = 0.01) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "identity") +
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(nonmem_data$time), by = 24*7),
                       labels = seq(0, max(nonmem_data$time)/(24*7), by = 1),
                       limits = c(0, max(nonmem_data$time))) +
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
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(nonmem_data$time), by = 24*7),
                       labels = seq(0, max(nonmem_data$time)/(24*7), by = 1),
                       limits = c(0, max(nonmem_data$time))) +
    theme_bw())

p_vpc + 
  p_pcvpc
