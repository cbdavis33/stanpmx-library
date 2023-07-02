rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_prop_covariates.rds")

nonmem_data <- read_csv(
  "depot_1cmt_linear_covariates/Data/depot_1cmt_prop_covariates.csv", 
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

new_data_to_simulate <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  expand(time = seq(time[1], time[2], by = 0.5)) %>%
  # expand(time = seq(time[1], time[2], by = 2)) %>%
  ungroup() %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 2,
         mdv = 1, 
         ss = 0,
         lloq = 1,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything()) %>% 
  left_join(nonmem_data %>% 
              select(ID, sexf, age, race, wt, cmppi, egfr) %>% 
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

cmppi <- new_data %>% 
  group_by(ID) %>% 
  distinct(cmppi) %>% 
  ungroup() %>% 
  pull(cmppi)

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
                  cmppi = cmppi,
                  egfr = egfr,
                  t_1 = 144,
                  t_2 = 168)

model <- cmdstan_model(
  "depot_1cmt_linear_covariates/Stan/Predict/depot_1cmt_prop_predict_observed_subjects_covariates.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))

post_preds_summary <- preds_df %>%
  spread_draws(pred[i], ipred[i], dv[i]) %>%
  mean_qi(pred, ipred, dv) %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         wt = new_data$wt[i],
         cmppi = new_data$cmppi[i],
         egfr = new_data$egfr[i]) %>% 
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
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, NA)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_trelliscope(~ ID, nrow = 3, ncol = 4, scales = "free")


tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 2, ncol = 2,
                               page = 1, scales = "free")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(post_preds_summary, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
          geom_line(aes(y = pred), linetype = 2, size = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          scale_y_log10(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                        limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, 168, by = 24),
                             labels = seq(0, 168, by = 24),
                             limits = c(0, 168)) +
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
  spread_draws(CL[ID], VC[ID], KA[ID], 
               auc_ss[ID], c_max[ID], t_max[ID], t_half[ID]) %>% 
  mean_qi() %>% 
  select(ID, CL, VC, KA, 
         auc_ss, c_max, t_max, t_half) %>% 
  inner_join(post_preds_summary %>% 
               filter(time == 168) %>% 
               select(ID, c_trough = "ipred") %>% 
               distinct(),
             by = "ID")

est_ind
