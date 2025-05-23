rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(mrgsolve)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_prop_covariates.rds")

nonmem_data <- read_csv("iv_2cmt_linear_covariates/Data/iv_2cmt_prop_covariates.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), # This value can be anything except NA. It'll be indexed away 
         cmt = 1)

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

end_inf_data <- nonmem_data %>% 
  filter(evid == 1) %>% 
  realize_addl() %>% 
  mutate(tinf = amt/rate,
         time = time + tinf) %>% 
  select(-tinf, -sexf, -age, -race, -race_asian, -wt, -egfr)

new_data_to_simulate <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  # expand(time = seq(time[1], time[2], by = 0.5)) %>%
  expand(time = seq(time[1], time[2], by = 12/24)) %>%
  ungroup() %>% 
  bind_rows(end_inf_data) %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 1,
         mdv = 1, 
         ss = 0,
         lloq = 1,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything()) %>% 
  left_join(nonmem_data %>% 
              select(ID, sexf, age, race, race_asian, wt, egfr) %>% 
              distinct(),
            by = "ID")

observed_data_to_simulate <- nonmem_data %>% 
  mutate(evid = if_else(evid == 1, 1, 2))

new_data <- new_data_to_simulate %>% 
  bind_rows(observed_data_to_simulate) %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, time, .keep_all = TRUE) %>% 
  select(-DV, -lloq, -bloq)

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

wt <- new_data %>% 
  group_by(ID) %>% 
  distinct(wt) %>% 
  ungroup() %>% 
  pull(wt)

race_asian <- new_data %>% 
  group_by(ID) %>% 
  distinct(race_asian) %>% 
  ungroup() %>% 
  pull(race_asian)

egfr <- new_data %>% 
  group_by(ID) %>% 
  distinct(egfr) %>% 
  ungroup() %>% 
  pull(egfr)

stan_data <- list(n_subjects = n_subjects,
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
                  wt = wt,
                  race_asian = race_asian,
                  egfr = egfr,
                  t_1 = 70,
                  t_2 = 84)

model <- cmdstan_model(
  "iv_2cmt_linear_covariates/Stan/Predict/iv_2cmt_prop_predict_observed_subjects_covariates.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

# preds <- model$generate_quantities(fit$draws() %>%
#                                      thin_draws(100),
#                                    data = stan_data,
#                                    parallel_chains = 4,
#                                    seed = 1234)

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))

post_preds_summary <- preds_df %>%
  spread_draws(pred[i], ipred[i], dv[i]) %>%
  mean_qi(pred, ipred, dv) %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>% 
  left_join(nonmem_data %>% 
              filter(mdv == 0) %>%
              select(ID, time, DV), 
            by = c("ID", "time")) %>% 
  select(ID, time, everything(), -i)

ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = pred), linetype = 2, size = 1.05) +
  geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
             color = "red") +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "log10",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (d)",
                     breaks = seq(0, max(nonmem_data$time), by = 14),
                     labels = seq(0, max(nonmem_data$time), by = 14),
                     limits = c(0, max(nonmem_data$time))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_trelliscope(~ ID, nrow = 2, ncol = 2, scales = "free")


tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 2, ncol = 2,
                               page = 1, scales = "free")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(post_preds_summary, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
          geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          scale_y_log10(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                        limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, max(nonmem_data$time), by = 14),
                             labels = seq(0, max(nonmem_data$time), by = 14),
                             limits = c(0, max(nonmem_data$time))) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 2, ncol = 2,
                                       page = i, scales = "free"))
  
}

## Individual estimates (posterior mean)
est_ind <- preds_df %>%
  spread_draws(c(CL, VC, Q, VP,
                 auc_ss, t_half_alpha, t_half_terminal)[ID]) %>% 
  mean_qi() %>% 
  select(ID, CL, VC, Q, VP, 
         auc_ss, t_half_alpha, t_half_terminal) %>% 
  inner_join(post_preds_summary %>% 
               filter(time == 84) %>% 
               select(ID, c_trough = "ipred") %>% 
               distinct(),
             by = "ID") %>% 
  inner_join(post_preds_summary %>% 
               filter(time == 70 + 1/24) %>% 
               select(ID, c_max = "ipred") %>% 
               distinct(),
             by = "ID")

est_ind
