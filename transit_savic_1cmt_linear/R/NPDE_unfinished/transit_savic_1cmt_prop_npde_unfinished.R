rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(vpc)
library(npde)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("transit_savic_1cmt_linear/Stan/Fits/transit_savic_1cmt_prop.rds")

nonmem_data <- read_csv("transit_savic_1cmt_linear/Data/transit_savic_1cmt_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

new_data <- nonmem_data %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, time, evid, .keep_all = TRUE) %>% 
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

new_data_dose <- new_data %>% 
  filter(evid == 1)

n_dose <- new_data_dose %>% 
  nrow()

subj_start_dose <- new_data_dose %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end_dose <- c(subj_start_dose[-1] - 1, n_dose) 

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
                  n_dose = n_dose,
                  dosetime = new_data_dose$time,
                  doseamt = new_data_dose$amt,
                  subj_start_dose = subj_start_dose,
                  subj_end_dose = subj_end_dose,
                  t_1 = 0,
                  t_2 = 24)

model <- cmdstan_model(
  "transit_savic_1cmt_linear/Stan/Predict/transit_savic_1cmt_prop_predict_new_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 123456) 

preds_df <- preds$draws(format = "draws_df")

ipreds_col <- preds_df %>% 
  spread_draws(ipred[i]) %>% 
  summarize(ipred_mean = mean(ipred)) %>% 
  ungroup() %>% 
  bind_cols(new_data %>% 
              select(ID, time, evid)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, ipred_mean, -i)

obs <- nonmem_data %>% 
  select(ID, time, dv = "DV", amt, evid, bloq, mdv) %>% 
  filter(evid == 0) %>% 
  mutate(dv = if_else(bloq == 1, 1, dv)) %>% 
  left_join(ipreds_col, by = c("ID", "time"))

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
  arrange(.draw, ID, time) %>% 
  select(ID, time, dv, bloq) %>% 
  left_join(ipreds_col, by = c("ID", "time"))

my_npde <- autonpde(namobs = obs %>% 
                      select(ID, time, dv, bloq, ipred_mean), 
                    namsim = sim, 
                    iid = "ID", ix = "time", iy = "dv", icens = "bloq",
                    iipred = "ipred_mean",
                    cens.method = "omit")

my_npde@results@res %>% 
  ggplot() + 
  geom_density(aes(x = npde), fill = "gray") + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  theme_bw()

plot(my_npde, plot.type = "data", 
     xlab = "Time (hr)", ylab = "Drug Conc. (ug/mL)", 
     line.loq = TRUE, ylim = c(0, NA), 
     main = "LOQ imputed to individual prediction")

plot(my_npde, plot.type = "vpc", 
     xlab = "Time (hr)", ylab = "Drug Conc. (ug/mL)", 
     line.loq = TRUE, ylim = c(0, NA), 
     main = "VPC")

