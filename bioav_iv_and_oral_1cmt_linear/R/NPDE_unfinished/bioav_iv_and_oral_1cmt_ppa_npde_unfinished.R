rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(npde)
library(patchwork)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds(
  "bioav_iv_and_oral_1cmt_linear/Stan/Fits/bioav_iv_and_oral_1cmt_ppa.rds")

nonmem_data <- read_csv(
  "bioav_iv_and_oral_1cmt_linear/Data/bioav_iv_and_oral_1cmt_ppa.csv",
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
                  t_2 = 168)

model <- cmdstan_model(
  "bioav_iv_and_oral_1cmt_linear/Stan/Predict/bioav_iv_and_oral_1cmt_ppa_predict_new_subjects.stan")

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

