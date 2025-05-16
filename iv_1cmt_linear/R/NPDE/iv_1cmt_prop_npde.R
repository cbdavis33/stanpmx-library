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

fit <- read_rds("iv_1cmt_linear/Stan/Fits/iv_1cmt_prop.rds")

nonmem_data <- read_csv("iv_1cmt_linear/Data/iv_1cmt_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), # This value can be anything except NA. It'll be indexed away 
         cmt = 1)

# NPDEs can be done two ways, and they're the same up to Monte Carlo error
#   1) Use the _predict_new_subjects.stan file 
#   2) Use the fitted object that already has already done everything in
#        generated quantities

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
                  t_1 = 70,
                  t_2 = 84,
                  want_auc_cmax = 0)

model <- cmdstan_model(
  "iv_1cmt_linear/Stan/Predict/iv_1cmt_prop_predict_new_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

data_obs <- nonmem_data %>% 
  filter(evid == 0) %>% 
  mutate(DV = if_else(bloq == 1, lloq, DV)) 

simdata <- preds_df %>%
  spread_draws(epred[i]) %>% 
  ungroup() %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         mdv = new_data$mdv[i],
         evid = new_data$evid[i],
         bloq = new_data$bloq[i]) %>% 
  filter(evid == 0) %>%
  arrange(.draw, ID, time) 

my_npde <- autonpde(data_obs %>% 
                      filter(evid == 0) %>% 
                      select(ID, amt, time, DV), 
                    simdata %>% 
                      select(ID, time, epred), 
                    1, 3, 4, boolsave = FALSE, cens.method = "cdf")

results_npde <- data_obs %>% 
  select(ID, lloq, bloq, time, DV) %>% 
  bind_cols(my_npde@results@res %>% 
              tibble() %>% 
              rename(ewres = ydobs) %>% 
              select(ypred, ycomp, ewres, npde))

# NPDE density plot
results_npde %>% 
  ggplot() + 
  geom_density(aes(x = npde), fill = "gray") + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  scale_x_continuous(name = "NPDE") +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

# NPDE vs. Population Prediction
results_npde %>% 
  ggplot(aes(x = ypred, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Population Prediction (EPRED)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()

# NPDE vs. Time
results_npde %>% 
  ggplot(aes(x = time, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Time (d)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()


## Alternatively, use the output in the fitted object

draws_df <- fit$draws(format = "draws_df")

data_obs <- nonmem_data %>% 
  filter(evid == 0) %>% 
  mutate(DV = if_else(bloq == 1, lloq, DV)) 

simdata <- draws_df %>%
  spread_draws(epred[i]) %>% 
  ungroup() %>% 
  mutate(ID = data_obs$ID[data_obs$evid == 0][i],
         time = data_obs$time[data_obs$evid == 0][i]) %>% 
  left_join(data_obs %>% 
              filter(evid == 0) %>%
              select(ID, time, DV), 
            by = c("ID", "time")) %>%
  select(ID, time, everything()) %>% 
  arrange(.draw, ID, time)

my_npde <- autonpde(data_obs %>% 
                      filter(evid == 0) %>% 
                      select(ID, amt, time, DV), 
                    simdata %>% 
                      select(ID, time, epred), 
                    1, 3, 4, boolsave = FALSE, cens.method = "cdf")

results_npde <- data_obs %>% 
  select(ID, lloq, bloq, time, DV) %>% 
  bind_cols(my_npde@results@res %>% 
              tibble() %>% 
              rename(ewres = ydobs) %>% 
              select(ypred, ycomp, ewres, npde))

# NPDE density plot
results_npde %>% 
  ggplot() + 
  geom_density(aes(x = npde), fill = "gray") + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  scale_x_continuous(name = "NPDE") +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

# NPDE vs. Population Prediction
results_npde %>% 
  ggplot(aes(x = ypred, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Population Prediction (EPRED)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()

# NPDE vs. Time
results_npde %>% 
  ggplot(aes(x = time, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Time (h)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()

